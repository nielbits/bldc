/*
	Copyright 2016 - 2022 Benjamin Vedder	benjamin@vedder.se

	This file is part of the VESC firmware.

	The VESC firmware is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    The VESC firmware is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include "foc_math.h"
#include "utils_math.h"
#include "encoder/encoder.h"
#include <math.h>

// See http://cas.ensmp.fr/~praly/Telechargement/Journaux/2010-IEEE_TPEL-Lee-Hong-Nam-Ortega-Praly-Astolfi.pdf
void foc_observer_update(float v_alpha, float v_beta, float i_alpha, float i_beta,
		float dt, observer_state *state, float *phase, motor_all_state_t *motor) {

	mc_configuration *conf_now = motor->m_conf;

	float R = conf_now->foc_motor_r;
	float L = conf_now->foc_motor_l;
	float lambda = conf_now->foc_motor_flux_linkage;

	// Saturation compensation
	switch(conf_now->foc_sat_comp_mode) {
	case SAT_COMP_LAMBDA:
		// Here we assume that the inductance drops by the same amount as the flux linkage. I have
		// no idea if this is a valid or even a reasonable assumption.
		if (conf_now->foc_observer_type >= FOC_OBSERVER_ORTEGA_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXLEMMING_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			L = L * (state->lambda_est / lambda);
		}
		break;

	case SAT_COMP_FACTOR: {
		const float comp_fact = conf_now->foc_sat_comp * (motor->m_motor_state.i_abs_filter / conf_now->l_current_max);
		L -= L * comp_fact;
		lambda -= lambda * comp_fact;
	} break;

	case SAT_COMP_LAMBDA_AND_FACTOR: {
		if (conf_now->foc_observer_type >= FOC_OBSERVER_ORTEGA_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXLEMMING_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type >= FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			L = L * (state->lambda_est / lambda);
		}
		const float comp_fact = conf_now->foc_sat_comp * (motor->m_motor_state.i_abs_filter / conf_now->l_current_max);
		L -= L * comp_fact;
	} break;

	default:
		break;
	}

	// Temperature compensation
	if (conf_now->foc_temp_comp) {
		R = motor->m_res_temp_comp;
	}

	float ld_lq_diff = conf_now->foc_motor_ld_lq_diff;
	float id = motor->m_motor_state.id;
	float iq = motor->m_motor_state.iq;

	// Adjust inductance for saliency.
	if (fabsf(id) > 0.1 || fabsf(iq) > 0.1) {
		L = L - ld_lq_diff / 2.0 + ld_lq_diff * SQ(iq) / (SQ(id) + SQ(iq));
	}

	float L_ia = L * i_alpha;
	float L_ib = L * i_beta;
	const float R_ia = R * i_alpha;
	const float R_ib = R * i_beta;
	const float gamma_half = motor->m_gamma_now * 0.5;

	switch (conf_now->foc_observer_type) {
	case FOC_OBSERVER_ORTEGA_ORIGINAL: {
		float err = SQ(lambda) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));

		// Forcing this term to stay negative helps convergence according to
		//
		// http://cas.ensmp.fr/Publications/Publications/Papers/ObserverPermanentMagnet.pdf
		// and
		// https://arxiv.org/pdf/1905.00833.pdf
		if (err > 0.0) {
			err = 0.0;
		}

		float x1_dot = v_alpha - R_ia + gamma_half * (state->x1 - L_ia) * err;
		float x2_dot = v_beta - R_ib + gamma_half * (state->x2 - L_ib) * err;

		state->x1 += x1_dot * dt;
		state->x2 += x2_dot * dt;
	} break;

	case FOC_OBSERVER_MXLEMMING:
	case FOC_OBSERVER_MXLEMMING_LAMBDA_COMP:
		// LICENCE NOTE:
		// This function deviates slightly from the BSD 3 clause licence.
		// The work here is entirely original to the MESC FOC project, and not based
		// on any appnotes, or borrowed from another project. This work is free to
		// use, as granted in BSD 3 clause, with the exception that this note must
		// be included in where this code is implemented/modified to use your
		// variable names, structures containing variables or other minor
		// rearrangements in place of the original names I have chosen, and credit
		// to David Molony as the original author must be noted.

		state->x1 += (v_alpha - R_ia) * dt - L * (i_alpha - state->i_alpha_last);
		state->x2 += (v_beta - R_ib) * dt - L * (i_beta - state->i_beta_last);

		if (conf_now->foc_observer_type == FOC_OBSERVER_MXLEMMING_LAMBDA_COMP) {
			float err = SQ(state->lambda_est) - (SQ(state->x1) + SQ(state->x2));
			state->lambda_est += 0.1 * gamma_half * state->lambda_est * -err * dt;
			utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

			utils_truncate_number_abs(&(state->x1), state->lambda_est);
			utils_truncate_number_abs(&(state->x2), state->lambda_est);
		} else {
			utils_truncate_number_abs(&(state->x1), lambda);
			utils_truncate_number_abs(&(state->x2), lambda);
		}

		// Set these to 0 to allow using the same atan2-code as for Ortega
		L_ia = 0.0;
		L_ib = 0.0;
		break;

	case FOC_OBSERVER_ORTEGA_LAMBDA_COMP: {
		float err = SQ(state->lambda_est) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));

		// FLux linkage observer. See:
		// https://cas.mines-paristech.fr/~praly/Telechargement/Conferences/2017_IFAC_Bernard-Praly.pdf
		state->lambda_est += 0.2 * gamma_half * state->lambda_est * -err * dt;

		// Clamp the observed flux linkage (not sure if this is needed)
		utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

		if (err > 0.0) {
			err = 0.0;
		}

		float x1_dot = v_alpha - R_ia + gamma_half * (state->x1 - L_ia) * err;
		float x2_dot = v_beta - R_ib + gamma_half * (state->x2 - L_ib) * err;

		state->x1 += x1_dot * dt;
		state->x2 += x2_dot * dt;
	} break;

	case FOC_OBSERVER_MXV:
	case FOC_OBSERVER_MXV_LAMBDA_COMP:
	case FOC_OBSERVER_MXV_LAMBDA_COMP_LIN:
		state->x1 += (v_alpha - R_ia) * dt;
		state->x2 += (v_beta - R_ib) * dt;

		if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP ||
				conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
			if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP_LIN) {
				float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
				UTILS_LP_FAST(state->lambda_est, mag, 0.1 * gamma_half * dt * SQ(state->lambda_est));
				utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

				if (mag > state->lambda_est) {
					state->x1 = (state->x1 / mag) * state->lambda_est;
					state->x2 = (state->x2 / mag) * state->lambda_est;
				}
			} else if (conf_now->foc_observer_type == FOC_OBSERVER_MXV_LAMBDA_COMP) {
				float err = SQ(state->lambda_est) - (SQ(state->x1 - L_ia) + SQ(state->x2 - L_ib));
				state->lambda_est += 0.2 * gamma_half * state->lambda_est * -err * dt;
				utils_truncate_number(&(state->lambda_est), lambda * 0.3, lambda * 2.5);

				float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
				if (mag > state->lambda_est) {
					state->x1 = (state->x1 / mag) * state->lambda_est;
					state->x2 = (state->x2 / mag) * state->lambda_est;
				}
			}
		} else {
			float mag = NORM2_f(state->x1 - L_ia, state->x2 - L_ib);
			if (mag > lambda) {
				state->x1 = (state->x1 / mag) * lambda;
				state->x2 = (state->x2 / mag) * lambda;
			}
		}
		break;

	default:
		break;
	}

	state->i_alpha_last = i_alpha;
	state->i_beta_last = i_beta;

	UTILS_NAN_ZERO(state->x1);
	UTILS_NAN_ZERO(state->x2);

	// Prevent the magnitude from getting too low, as that makes the angle very unstable.
	float mag = NORM2_f(state->x1, state->x2);
	if (mag < (lambda * 0.5)) {
		state->x1 *= 1.1;
		state->x2 *= 1.1;
	}

	if (phase) {
		*phase = utils_fast_atan2(state->x2 - L_ib, state->x1 - L_ia);
	}

	// Can we clamp the flux in dq with q flux = 0 and d flux is lambda
	// Then the state->x1 and state->x2 (which are the alpha and beta fluxes) are set as lambda*sin and lambda*cos
	// The d flux each time would have a residual after transform from ab to dq. This can be used as an input to the flux estimator
}

void foc_pll_run(float phase, float dt, float *phase_var,
					float *speed_var, mc_configuration *conf) {
	UTILS_NAN_ZERO(*phase_var);
	float delta_theta = phase - *phase_var;
	utils_norm_angle_rad(&delta_theta);
	UTILS_NAN_ZERO(*speed_var);
	*phase_var += (*speed_var + conf->foc_pll_kp * delta_theta) * dt;
	utils_norm_angle_rad((float*)phase_var);
	*speed_var += conf->foc_pll_ki * delta_theta * dt;
}

/**
 * @brief svm Space vector modulation. Magnitude must not be larger than sqrt(3)/2, or 0.866 to avoid overmodulation.
 *        See https://github.com/vedderb/bldc/pull/372#issuecomment-962499623 for a full description.
 * @param alpha voltage
 * @param beta Park transformed and normalized voltage
 * @param PWMFullDutyCycle is the peak value of the PWM counter.
 * @param tAout PWM duty cycle phase A (0 = off all of the time, PWMFullDutyCycle = on all of the time)
 * @param tBout PWM duty cycle phase B
 * @param tCout PWM duty cycle phase C
 */
void foc_svm(float alpha, float beta, uint32_t PWMFullDutyCycle,
				uint32_t* tAout, uint32_t* tBout, uint32_t* tCout, uint32_t *svm_sector) {
	uint32_t sector;

	if (beta >= 0.0f) {
		if (alpha >= 0.0f) {
			//quadrant I
			if (ONE_BY_SQRT3 * beta > alpha) {
				sector = 2;
			} else {
				sector = 1;
			}
		} else {
			//quadrant II
			if (-ONE_BY_SQRT3 * beta > alpha) {
				sector = 3;
			} else {
				sector = 2;
			}
		}
	} else {
		if (alpha >= 0.0f) {
			//quadrant IV5
			if (-ONE_BY_SQRT3 * beta > alpha) {
				sector = 5;
			} else {
				sector = 6;
			}
		} else {
			//quadrant III
			if (ONE_BY_SQRT3 * beta > alpha) {
				sector = 4;
			} else {
				sector = 5;
			}
		}
	}

	// PWM timings
	uint32_t tA, tB, tC;

	switch (sector) {

	// sector 1-2
	case 1: {
		// Vector on-times
		uint32_t t1 = (alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t2 = (TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tA = (PWMFullDutyCycle + t1 + t2) / 2;
		tB = tA - t1;
		tC = tB - t2;

		break;
	}

	// sector 2-3
	case 2: {
		// Vector on-times
		uint32_t t2 = (alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t3 = (-alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tB = (PWMFullDutyCycle + t2 + t3) / 2;
		tA = tB - t3;
		tC = tA - t2;

		break;
	}

	// sector 3-4
	case 3: {
		// Vector on-times
		uint32_t t3 = (TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t4 = (-alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tB = (PWMFullDutyCycle + t3 + t4) / 2;
		tC = tB - t3;
		tA = tC - t4;

		break;
	}

	// sector 4-5
	case 4: {
		// Vector on-times
		uint32_t t4 = (-alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t5 = (-TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tC = (PWMFullDutyCycle + t4 + t5) / 2;
		tB = tC - t5;
		tA = tB - t4;

		break;
	}

	// sector 5-6
	case 5: {
		// Vector on-times
		uint32_t t5 = (-alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t6 = (alpha - ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tC = (PWMFullDutyCycle + t5 + t6) / 2;
		tA = tC - t5;
		tB = tA - t6;

		break;
	}

	// sector 6-1
	case 6: {
		// Vector on-times
		uint32_t t6 = (-TWO_BY_SQRT3 * beta) * PWMFullDutyCycle;
		uint32_t t1 = (alpha + ONE_BY_SQRT3 * beta) * PWMFullDutyCycle;

		// PWM timings
		tA = (PWMFullDutyCycle + t6 + t1) / 2;
		tC = tA - t1;
		tB = tC - t6;

		break;
	}
	}

	*tAout = tA;
	*tBout = tB;
	*tCout = tC;
	*svm_sector = sector;
}

void foc_run_pid_control_pos(bool index_found, float dt, motor_all_state_t *motor) {
	mc_configuration *conf_now = motor->m_conf;

	float angle_now = motor->m_pos_pid_now;
	float angle_set = motor->m_pos_pid_set;

	float p_term;
	float d_term;
	float d_term_proc;

	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_POS) {
		motor->m_pos_i_term = 0;
		motor->m_pos_prev_error = 0;
		motor->m_pos_prev_proc = angle_now;
		motor->m_pos_d_filter = 0.0;
		motor->m_pos_d_filter_proc = 0.0;
		return;
	}

	// Compute parameters
	float error = utils_angle_difference(angle_set, angle_now);
	float error_sign = 1.0;

	if (conf_now->m_sensor_port_mode != SENSOR_PORT_MODE_HALL) {
		if (conf_now->foc_encoder_inverted) {
			error_sign = -1.0;
		}
	}

	error *= error_sign;

	float kp = conf_now->p_pid_kp;
	float ki = conf_now->p_pid_ki;
	float kd = conf_now->p_pid_kd;
	float kd_proc = conf_now->p_pid_kd_proc;

	if (conf_now->p_pid_gain_dec_angle > 0.1) {
		float min_error = conf_now->p_pid_gain_dec_angle / conf_now->p_pid_ang_div;
		float error_abs = fabs(error);

		if (error_abs < min_error) {
			float scale = error_abs / min_error;
			kp *= scale;
			ki *= scale;
			kd *= scale;
			kd_proc *= scale;
		}
	}

	p_term = error * kp;
	motor->m_pos_i_term += error * (ki * dt);

	// Average DT for the D term when the error does not change. This likely
	// happens at low speed when the position resolution is low and several
	// control iterations run without position updates.
	// TODO: Are there problems with this approach?
	motor->m_pos_dt_int += dt;
	if (error == motor->m_pos_prev_error) {
		d_term = 0.0;
	} else {
		d_term = (error - motor->m_pos_prev_error) * (kd / motor->m_pos_dt_int);
		motor->m_pos_dt_int = 0.0;
	}

	// Filter D
	UTILS_LP_FAST(motor->m_pos_d_filter, d_term, conf_now->p_pid_kd_filter);
	d_term = motor->m_pos_d_filter;

	// Process D term
	motor->m_pos_dt_int_proc += dt;
	if (angle_now == motor->m_pos_prev_proc) {
		d_term_proc = 0.0;
	} else {
		d_term_proc = -utils_angle_difference(angle_now, motor->m_pos_prev_proc) * error_sign * (kd_proc / motor->m_pos_dt_int_proc);
		motor->m_pos_dt_int_proc = 0.0;
	}

	// Filter D process
	UTILS_LP_FAST(motor->m_pos_d_filter_proc, d_term_proc, conf_now->p_pid_kd_filter);
	d_term_proc = motor->m_pos_d_filter_proc;

	// I-term wind-up protection
	float p_tmp = p_term;
	utils_truncate_number_abs(&p_tmp, 1.0);
	utils_truncate_number_abs((float*)&motor->m_pos_i_term, 1.0 - fabsf(p_tmp));

	// Store previous error
	motor->m_pos_prev_error = error;
	motor->m_pos_prev_proc = angle_now;

	// Calculate output
	float output = p_term + motor->m_pos_i_term + d_term + d_term_proc;
	utils_truncate_number(&output, -1.0, 1.0);

	if (conf_now->m_sensor_port_mode != SENSOR_PORT_MODE_HALL) {
		if (index_found) {
			motor->m_iq_set = output * conf_now->l_current_max * conf_now->l_current_max_scale;;
		} else {
			// Rotate the motor with 40 % power until the encoder index is found.
			motor->m_iq_set = 0.4 * conf_now->l_current_max * conf_now->l_current_max_scale;;
		}
	} else {
		motor->m_iq_set = output * conf_now->l_current_max * conf_now->l_current_max_scale;;
	}
}

void foc_run_pid_control_speed(bool index_found, float dt, motor_all_state_t *motor) {
	mc_configuration *conf_now = motor->m_conf;
	float p_term;
	float d_term;

	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_SPEED) {
		motor->m_speed_i_term = 0.0;
		motor->m_speed_prev_error = 0.0;
		motor->m_speed_d_filter = 0.0;
		return;
	}

	if (conf_now->s_pid_ramp_erpms_s > 0.0) {
		utils_step_towards((float*)&motor->m_speed_pid_set_rpm, motor->m_speed_command_rpm, conf_now->s_pid_ramp_erpms_s * dt);
		if (!index_found) {
			utils_truncate_number_abs(&motor->m_speed_pid_set_rpm, conf_now->foc_openloop_rpm);
		}
		utils_truncate_number(&motor->m_speed_pid_set_rpm, conf_now->l_min_erpm, conf_now->l_max_erpm);
	}

	float rpm = 0.0;
	switch (conf_now->s_pid_speed_source) {
	case S_PID_SPEED_SRC_PLL:
		rpm = RADPS2RPM_f(motor->m_pll_speed);
		break;
	case S_PID_SPEED_SRC_FAST:
		rpm = RADPS2RPM_f(motor->m_speed_est_fast);
		break;
	case S_PID_SPEED_SRC_FASTER:
		rpm = RADPS2RPM_f(motor->m_speed_est_faster);
		break;
	}

	//start KF 2 state

	float angle_deg_now = encoder_read_deg();
	float angle_rad_now = angle_deg_now * (M_PI / 180.0f);  // Convert to radians

	// Compute delta using unwrap-safe function (handles wrap around 2π)
	float delta_rad = utils_angle_difference_rad(angle_rad_now, motor->kalman_last_angle_rad);

	// Store the new angle
	motor->kalman_last_angle_rad = angle_rad_now;

	// Accumulate unwrapped angle
	motor->unwrapped_theta += delta_rad;

	// Save delta for next iteration
	motor->kalman_last_delta_rad = delta_rad;


	// possibly changeable parameters
	float gear_ratio = motor->gear_ratio_bike;
	float incline    = 0.000f; // inclination angle, degrees
	float gearing    = (float)(motor->p_mech_gearing / gear_ratio);  // motor to wheel, /gearin. wheel to motor, *gearin, if speed
																	//if torque, motor to wheel * gearin, wheel to motor / gearin
	float slope      = incline * 3.14159265359f / 180.0f;            // radians

		// --- speed (fix: ERPM -> mech RPM -> rad/s -> m/s) ---
		//const float pole_pairs = motor->m_conf->si_motor_poles * 0.5f;   // poles -> pole-pairs
		//float omega_m   = motor->m_speed_est_fast / pole_pairs;            // mech motor RPM
		//float omega_w   = (rpm_m / gearing);    // wheel RPM (divide by Gm, multiply by Gbike)
		//float speed=omega_w*motor->p_wheel_radius
	float speed   = motor->d_speed_soll;                 // m/s  (kept name)

	// --- forces (same variable names, corrected formulas) ---
	float F_air     = 0.5f * motor->p_air_ro *0.9f*0.5f*speed*fabsf(speed);//* motor->p_c_air * motor->p_As * speed * speed*SIGN(speed);          // 0.5 * rho * Cd * A * v^2
	float F_roll    = smooth_force((motor->p_c_rr * motor->p_weight * 9.81f * cosf(slope)),speed,0.1f);                          // Crr * m g cos(theta)
	float F_incline = motor->p_weight * 9.81f * sinf(slope);                                          // m g sin(theta)  (set incline=0 if you want it off)
	// simple viscous bearing drag in force-domain (N·s/m). keep name, fix units:
	float F_bearings = smooth_force(((motor->p_c_bw * motor->p_k_v_bw) * speed),speed,0.1f);

	// F_res calculation
	float F_combine = F_air + F_roll + F_incline + F_bearings; // resistance force

	

	#define SCALE_INT 10000000.0f   // Float version for scaling
	//#define OBS_GAIN_FP 5000000000LL    // Observer gain (e.g. 0.5 scaled to 1e8)

	motor->omega_fp = (int_fast64_t)(motor->m_speed_est_fast_corrected * SCALE_INT);

	UTILS_LP_FAST_I64(&motor->omega_filtered_fp, motor->omega_fp, 30);

	// --- Calculate motor torque ---
	motor->te_calculated = motor->m_motor_state.iq * motor->p_kT + (motor->m_motor_state.iq * motor->m_motor_state.id) * (motor->p_ld - motor->p_lq);


						
	//Kalman Filter 3 Steps implemented


	//const float T_fric_ff = Tf_smooth(motor->kalman_x[1], -0.31f, 0.000149291f, 3.9f);
	motor->unwrapped_theta_filtered=UTILS_LP_FAST(motor->unwrapped_theta_filtered,motor->unwrapped_theta,0.05f);
	ekf3_step_simple(
    motor,
    dt,
    motor->te_calculated,   // Nm (you already compute this above)
    motor->unwrapped_theta // rad (your unwrapped mech angle)
	);

	// Use EKF outputs
	const float theta_hat = motor->kalman_x[0];
	const float omega_hat = motor->kalman_x[1];
	const float Tp_hat    = motor->kalman_x[2];

	motor->ekf_rpm= omega_hat * 9.54929f * (motor->m_conf->si_motor_poles / 2.0f);

	//float error = motor->m_speed_pid_set_rpm - rpm;
	//use the internal speed_soll directly as reference(no exchange to matlab needed, 40x faster)
	float error;


	error = motor->d_erpm_soll - motor->ekf_rpm;
	



	// Too low RPM set. Reset state, release motor and return.
	if (fabsf(motor->m_speed_pid_set_rpm) < conf_now->s_pid_min_erpm) {
		motor->m_speed_i_term = 0.0;
		motor->m_speed_prev_error = error;
		motor->m_iq_set = 0.0;
		return;
	}

	// Compute parameters
	p_term = error * conf_now->s_pid_kp * (1.0 / 20.0);
	d_term = (error - motor->m_speed_prev_error) * (conf_now->s_pid_kd / dt) * (1.0 / 20.0);

	// Filter D
	UTILS_LP_FAST(motor->m_speed_d_filter, d_term, conf_now->s_pid_kd_filter);
	d_term = motor->m_speed_d_filter;

	// Store previous error
	motor->m_speed_prev_error = error;
	// Keep your original semantics for tp_observed:


	motor->tp_observed = Tp_hat;

	motor->d_f_motor = (motor->te_calculated) / motor->p_wheel_radius * motor->p_mech_gearing;



	// --- Acceleration for integration (still uses full force model) ---
	motor->accel_ist = ((motor->tp_observed/(motor->p_wheel_radius)*gearing - (F_combine) )  ) * SCALE_INT/motor->p_weight; // m/s^2 scaled

	

	// --- Position integration (fixed-point) ---
	motor->integrated_value = (int_fast64_t)(((motor->last_accel + motor->accel_ist) * (dt * SCALE_INT)) / (SCALE_INT * 2)) + motor->integrated_value;
	//if (motor->integrated_value < 0) {
    //motor->integrated_value = 0;
	//}
	motor->last_accel = motor->accel_ist;
	motor->d_speed_soll = (float)(motor->integrated_value / SCALE_INT);


	//speed setpoint calculation

	float v = (float)motor->integrated_value / SCALE_INT; // linear speed in m/s
	float omega_mech = v *  (gearing) / motor->p_wheel_radius;  // rad/s
	motor->d_erpm_soll = omega_mech  * 9.54929 * (motor->m_conf->si_motor_poles/2.0f);

	// i_res calculation

	float T_res=F_combine*motor->p_wheel_radius/gearing;
	
	float i_res= -(T_res)/motor->p_kT;//motor->p_kT;
	
	//float lq= motor->m_conf->foc_motor_l+motor->m_conf->foc_motor_ld_lq_diff/2.0;
	float i_res_out= i_res/(conf_now->lo_current_max * conf_now->l_current_max_scale);


	motor->c_v_q_ff= i_res_out*motor->m_res_est; //+ motor->m_speed_est_fast*lq*i_res_out; 

	// Other motor variables remain unchanged
	motor->d_speed = motor->m_speed_est_fast*motor->p_wheel_radius/(motor->m_conf->si_motor_poles)/gearing;//speed in m/s;
	motor->d_f_air = F_air;
	motor->d_f_combine = F_combine;
	motor->d_f_bearings = F_bearings;
	motor->d_i_res= i_res;
	// Calculate output
	float output = p_term + motor->m_speed_i_term + d_term;
	utils_truncate_number_abs(&output, 1.0);

	// Integrator windup protection
	motor->m_speed_i_term += error * conf_now->s_pid_ki * dt * (1.0 / 20.0);
	utils_truncate_number_abs(&motor->m_speed_i_term, 1.0);

	if (conf_now->s_pid_ki < 1e-9) {
		motor->m_speed_i_term = 0.0;
	}

	// Optionally disable braking
	if (!conf_now->s_pid_allow_braking) {
		if (rpm > 20.0 && output < 0.0) {
			output = 0.0;
		}

		if (rpm < -20.0 && output > 0.0) {
			output = 0.0;
		}
	}

	motor->m_iq_set = output * conf_now->lo_current_max * conf_now->l_current_max_scale;
}



float foc_correct_encoder(float obs_angle, float enc_angle, float speed,
							 float sl_erpm, motor_all_state_t *motor) {
	float rpm_abs = fabsf(RADPS2RPM_f(speed));

	// Hysteresis 5 % of total speed
	float hyst = sl_erpm * 0.05;
	if (motor->m_using_encoder) {
		if (rpm_abs > (sl_erpm + hyst)) {
			motor->m_using_encoder = false;
		}
	} else {
		if (rpm_abs < (sl_erpm- hyst)) {
			motor->m_using_encoder = true;
		}
	}

	return motor->m_using_encoder ? enc_angle : obs_angle;
}

float foc_correct_hall(float angle, float dt, motor_all_state_t *motor, int hall_val) {
	mc_configuration *conf_now = motor->m_conf;
	motor->m_hall_dt_diff_now += dt;

	float rpm_abs = fabsf(RADPS2RPM_f(motor->m_pll_speed));
	float rad_per_sec_hall = (M_PI / 3.0) / motor->m_hall_dt_diff_last;
	float rpm_abs_hall = fabsf(RADPS2RPM_f(rad_per_sec_hall));

	motor->m_using_hall = rpm_abs < conf_now->foc_sl_erpm;
	float angle_old = angle;

	int ang_hall_int = conf_now->foc_hall_table[hall_val];

	// Only override the observer if the hall sensor value is valid.
	if (ang_hall_int < 201) {
		// Scale to the circle and convert to radians
		float ang_hall_now = ((float)ang_hall_int / 200.0) * 2.0 * M_PI;

		if (motor->m_ang_hall_int_prev < 0) {
			// Previous angle not valid
			motor->m_ang_hall_int_prev = ang_hall_int;
			motor->m_ang_hall = ang_hall_now;
		} else if (ang_hall_int != motor->m_ang_hall_int_prev) {
			int diff = ang_hall_int - motor->m_ang_hall_int_prev;
			if (diff > 100) {
				diff -= 200;
			} else if (diff < -100) {
				diff += 200;
			}

			// This is only valid if the direction did not just change. If it did, we use the
			// last speed together with the sign right now.
			if (SIGN(diff) == SIGN(motor->m_hall_dt_diff_last)) {
				if (diff > 0) {
					motor->m_hall_dt_diff_last = motor->m_hall_dt_diff_now;
				} else {
					motor->m_hall_dt_diff_last = -motor->m_hall_dt_diff_now;
				}
			} else {
				motor->m_hall_dt_diff_last = -motor->m_hall_dt_diff_last;
			}

			motor->m_hall_dt_diff_now = 0.0;

			// A transition was just made. The angle is in the middle of the new and old angle.
			int ang_avg = motor->m_ang_hall_int_prev + diff / 2;
			ang_avg %= 200;

			// Scale to the circle and convert to radians
			motor->m_ang_hall = ((float)ang_avg / 200.0) * 2.0 * M_PI;
		}

		motor->m_ang_hall_int_prev = ang_hall_int;

		if (RADPS2RPM_f((M_PI / 3.0) / fmaxf(fabsf(motor->m_hall_dt_diff_now),
				fabsf(motor->m_hall_dt_diff_last))) < conf_now->foc_hall_interp_erpm) {
			// Don't interpolate on very low speed, just use the closest hall sensor. The reason is that we might
			// get stuck at 60 degrees off if a direction change happens between two steps.
			motor->m_ang_hall = ang_hall_now;
		} else {
			// Interpolate
			float diff = utils_angle_difference_rad(motor->m_ang_hall, ang_hall_now);
			if (fabsf(diff) < ((2.0 * M_PI) / 12.0) || SIGN(diff) != SIGN(rad_per_sec_hall)) {
				// Do interpolation
				motor->m_ang_hall += rad_per_sec_hall * dt;
			} else {
				// We are too far away with the interpolation
				motor->m_ang_hall -= diff * 0.01;
			}
		}

		utils_norm_angle_rad((float*)&motor->m_ang_hall);

		// Limit hall sensor rate of change. This will reduce current spikes in the current controllers when the angle estimation
		// changes fast.
		float angle_step = (fmaxf(rpm_abs_hall, conf_now->foc_hall_interp_erpm) / 60.0) * 2.0 * M_PI * dt * 1.5;
		float angle_diff = utils_angle_difference_rad(motor->m_ang_hall, motor->m_ang_hall_rate_limited);
		if (fabsf(angle_diff) < angle_step) {
			motor->m_ang_hall_rate_limited = motor->m_ang_hall;
		} else {
			motor->m_ang_hall_rate_limited += angle_step * SIGN(angle_diff);
		}

		utils_norm_angle_rad((float*)&motor->m_ang_hall_rate_limited);

		if (motor->m_using_hall) {
			angle = motor->m_ang_hall_rate_limited;
		}
	} else {
		// Invalid hall reading. Don't update angle.
		motor->m_ang_hall_int_prev = -1;

		// Also allow open loop in order to behave like normal sensorless
		// operation. Then the motor works even if the hall sensor cable
		// gets disconnected (when the sensor spacing is 120 degrees).
		if (motor->m_phase_observer_override && motor->m_state == MC_STATE_RUNNING) {
			angle = motor->m_phase_now_observer_override;
		}
	}

	// Map output angle between hall angle and observer angle in transition region to make
	// a smooth transition.
	if (angle_old != angle) {
		float weight_hall = utils_map(rpm_abs, conf_now->foc_sl_erpm_start, conf_now->foc_sl_erpm, 1.0, 0.0);
		utils_truncate_number(&weight_hall, 0.0, 1.0);
		angle = utils_interpolate_angles_rad(angle, angle_old, weight_hall);
	}

	return angle;
}

void foc_run_fw(motor_all_state_t *motor, float dt) {
	if (motor->m_conf->foc_fw_current_max < fmaxf(motor->m_conf->cc_min_current, 0.001)) {
		return;
	}

	// Field Weakening
	// FW is used in the current and speed control modes. If a different mode is used
	// this code also runs if field weakening was active before. This allows
	// changing control mode even while in field weakening.
	if (motor->m_state == MC_STATE_RUNNING &&
			(motor->m_control_mode == CONTROL_MODE_CURRENT ||
					motor->m_control_mode == CONTROL_MODE_CURRENT_BRAKE ||
					motor->m_control_mode == CONTROL_MODE_SPEED ||
					motor->m_i_fw_set > motor->m_conf->cc_min_current)) {
		float fw_current_now = 0.0;
		float duty_abs = motor->m_duty_abs_filtered;

		if (motor->m_conf->foc_fw_duty_start < 0.99 &&
				duty_abs > motor->m_conf->foc_fw_duty_start * motor->m_conf->l_max_duty) {
			fw_current_now = utils_map(duty_abs,
					motor->m_conf->foc_fw_duty_start * motor->m_conf->l_max_duty,
					motor->m_conf->l_max_duty,
					0.0, motor->m_conf->foc_fw_current_max);

			// m_current_off_delay is used to not stop the modulation too soon after leaving FW. If axis decoupling
			// is not working properly an oscillation can occur on the modulation when changing the current
			// fast, which can make the estimated duty cycle drop below the FW threshold long enough to stop
			// modulation. When that happens the body diodes in the MOSFETs can see a lot of current and unexpected
			// braking happens. Therefore the modulation is left on for some time after leaving FW to give the
			// oscillation a chance to decay while the MOSFETs are still driven.
			motor->m_current_off_delay = 1.0;
		}

		if (motor->m_conf->foc_fw_ramp_time < dt) {
			motor->m_i_fw_set = fw_current_now;
		} else {
			utils_step_towards((float*)&motor->m_i_fw_set, fw_current_now,
					(dt / motor->m_conf->foc_fw_ramp_time) * motor->m_conf->foc_fw_current_max);
		}
	}
}

void foc_hfi_adjust_angle(float ang_err, motor_all_state_t *motor, float dt) {
	mc_configuration *conf = motor->m_conf;
	utils_truncate_number_abs(&ang_err, conf->foc_hfi_max_err);

	// TODO: Check if ratio between these is sane or introduce separate gains
	const float gain_int = 4000.0 * conf->foc_hfi_gain;
	const float gain_int2 = 10.0 * conf->foc_hfi_gain;
	motor->m_hfi.double_integrator += ang_err * gain_int2;
	utils_truncate_number_abs(&motor->m_hfi.double_integrator, fabsf(motor->m_speed_est_fast));
	motor->m_hfi.angle -= dt * (gain_int * ang_err + motor->m_hfi.double_integrator);
	utils_norm_angle_rad((float*)&motor->m_hfi.angle);
	motor->m_hfi.ready = true;
}

void foc_precalc_values(motor_all_state_t *motor) {
	const mc_configuration *conf_now = motor->m_conf;
	motor->p_lq = conf_now->foc_motor_l + conf_now->foc_motor_ld_lq_diff * 0.5;
	motor->p_ld = conf_now->foc_motor_l - conf_now->foc_motor_ld_lq_diff * 0.5;
	motor->p_inv_ld_lq = (1.0 / motor->p_lq - 1.0 / motor->p_ld);
	motor->p_v2_v3_inv_avg_half = (0.5 / motor->p_lq + 0.5 / motor->p_ld) * 0.9; // With the 0.9 we undo the adjustment from the detection
	motor->m_observer_state.lambda_est = conf_now->foc_motor_flux_linkage;
}


inline float sgn_db(float x, float dead) {
    if (x >  dead) return  1.0f;
    if (x < -dead) return -1.0f;
    return 0.0f;
}


static inline float dTf_domega(float omega, float b) {
    // derivative of b*|omega|; Coulomb derivative ~0 (a.e.)
    return (fabsf(omega) < 1e-4f) ? 0.0f : b * (omega > 0.0f ? 1.0f : -1.0f);
}
// --- EKF step: predict + angle update -------------------------------------
inline void ekf3_step_simple(
    motor_all_state_t *m,
    float dt,
    float Te,             // motor torque [Nm]
    float theta_meas      // unwrapped mechanical angle [rad]
) {
    // ---- Unpack ----
    float *x      = m->kalman_x;   // x = [theta, omega, T_p]
    float (*P)[3] = m->kalman_P;
    float (*Q)[3] = m->kalman_Q;
    const float J = m->p_J;

    if (dt <= 0.0f) return;

    // ---- Friction: Stribeck (smooth, consistent Jacobian) ----
    float Tf = 0.0f, dTf = 0.0f;
    stribeck_tf_and_dtf(m, x[1], &Tf, &dTf);

    // ---- Predict (nonlinear) ----
    const float invJ = 1.0f / J;
    const float b_dt = dt * invJ;

    const float th_n = x[0] + dt * x[1];
    const float om_n = x[1] + b_dt * (Te + x[2] - Tf);
    const float Tp_n = x[2];

    x[0] = th_n;
    x[1] = om_n;
    x[2] = Tp_n;

    // ---- Covariance predict with consistent F ----
    // F = [1, dt, 0;
    //      0, 1 - dt/J * dTf/dω, dt/J;
    //      0, 0, 1]
    const float a   = dt * invJ * dTf;
    const float F00 = 1.0f, F01 = dt,      F02 = 0.0f;
    const float F10 = 0.0f, F11 = 1.0f - a, F12 = b_dt;
    const float F20 = 0.0f, F21 = 0.0f,    F22 = 1.0f;

    // FP = F*P
    float FP00 = F00*P[0][0] + F01*P[1][0] + F02*P[2][0];
    float FP01 = F00*P[0][1] + F01*P[1][1] + F02*P[2][1];
    float FP02 = F00*P[0][2] + F01*P[1][2] + F02*P[2][2];

    float FP10 = F10*P[0][0] + F11*P[1][0] + F12*P[2][0];
    float FP11 = F10*P[0][1] + F11*P[1][1] + F12*P[2][1];
    float FP12 = F10*P[0][2] + F11*P[1][2] + F12*P[2][2];

    float FP20 = F20*P[0][0] + F21*P[1][0] + F22*P[2][0];
    float FP21 = F20*P[0][1] + F21*P[1][1] + F22*P[2][1];
    float FP22 = F20*P[0][2] + F21*P[1][2] + F22*P[2][2];

    // P = FP*F' + Q
    float P00 = FP00*F00 + FP01*F01 + FP02*F02 + Q[0][0];
    float P01 = FP00*F10 + FP01*F11 + FP02*F12 + Q[0][1];
    float P02 = FP00*F20 + FP01*F21 + FP02*F22 + Q[0][2];

    float P10 = FP10*F00 + FP11*F01 + FP12*F02 + Q[1][0];
    float P11 = FP10*F10 + FP11*F11 + FP12*F12 + Q[1][1];
    float P12 = FP10*F20 + FP11*F21 + FP12*F22 + Q[1][2];

    float P20 = FP20*F00 + FP21*F01 + FP22*F02 + Q[2][0];
    float P21 = FP20*F10 + FP21*F11 + FP22*F12 + Q[2][1];
    float P22 = FP20*F20 + FP21*F21 + FP22*F22 + Q[2][2];

    P[0][0]=P00; P[0][1]=P01; P[0][2]=P02;
    P[1][0]=P10; P[1][1]=P11; P[1][2]=P12;
    P[2][0]=P20; P[2][1]=P21; P[2][2]=P22;

    // ---- Angle update (H = [1 0 0]) ----
    {
        const float y    = theta_meas - x[0];
        const float S    = P[0][0] + m->kalman_R;
        const float invS = 1.0f / S;

        const float r0 = P[0][0], r1 = P[0][1], r2 = P[0][2];
        const float K0 = P[0][0] * invS;
        const float K1 = P[1][0] * invS;
        const float K2 = P[2][0] * invS;

        x[0] += K0 * y;
        x[1] += K1 * y;
        x[2] += K2 * y;

        P[0][0] -= K0 * r0;  P[0][1] -= K0 * r1;  P[0][2] -= K0 * r2;
        P[1][0] -= K1 * r0;  P[1][1] -= K1 * r1;  P[1][2] -= K1 * r2;
        P[2][0] -= K2 * r0;  P[2][1] -= K2 * r1;  P[2][2] -= K2 * r2;

        // keep symmetry
        P[1][0] = P[0][1];
        P[2][0] = P[0][2];
        P[2][1] = P[1][2];
    }

    // ---- Guards ----
    const float Pmin = 1e-12f;
    if (P[0][0] < Pmin) P[0][0] = Pmin;
    if (P[1][1] < Pmin) P[1][1] = Pmin;
    if (P[2][2] < Pmin) P[2][2] = Pmin;

    // export (optional)
    m->Tf_hat = Tf;
}

static inline void stribeck_tf_and_dtf(
    const motor_all_state_t *m,
    float omega,
    float *Tf_out,
    float *dTf_out
) {
    const float B      = m->fric_B;        // Nm·s/rad
    const float Tc     = m->fric_Tc;       // Nm
    const float Ts     = m->fric_Ts;       // Nm
    const float vs     = m->fric_vs;       // rad/s
    const float alpha  = m->fric_alpha;    // -
    const float eps    = m->fric_eps;      // rad/s   (for tanh)
    const float delta  = m->fric_delta;    // rad/s   (for |w| smoothing)

    // Smooth sign and smooth absolute value
    const float s      = tanhf(omega / eps);           // in (-1,1)
    const float sech2  = 1.0f - s * s;                 // d/dx tanh = sech^2 = 1 - tanh^2
    const float dsdw   = (1.0f / eps) * sech2;         // ds/dω

    const float wabs   = sqrtf(omega*omega + delta*delta); // ≥ delta
    // exp term e = exp( - ( (|w|/vs)^alpha ) )
    const float r      = wabs / vs;
    const float z      = powf(r, alpha);
    const float e      = expf(-z);

    // Amplitude A(w) = Tc + (Ts - Tc)*e
    const float A      = Tc + (Ts - Tc) * e;

    // dA/dω = (Ts - Tc) * de/dω
    // de/dω = -e * (alpha / vs^alpha) * wabs^(alpha-2) * omega
    //       = -e * alpha * omega * (wabs)^(alpha-2) / (vs^alpha)
    float dAdw = 0.0f;
    {
        const float vs_a   = powf(vs, alpha);
        const float w_pow  = powf(wabs, alpha - 2.0f);   // safe via delta
        dAdw = (Ts - Tc) * (-e) * alpha * omega * (w_pow / vs_a);
    }

    // Tf = B*omega + A(w)*s
    const float Tf  = B * omega + A * s;

    // dTf/dω = B + dA/dω * s + A * ds/dω
    const float dTf = B + dAdw * s + A * dsdw;

    if (Tf_out)  { *Tf_out  = Tf; }
    if (dTf_out) { *dTf_out = dTf; }
}

float smooth_force(float mag, float v, float v_eps) {
    return mag * tanhf(v / v_eps); // smoothly goes negative if v<0
}