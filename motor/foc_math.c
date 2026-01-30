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
	float kp, ki, kd, kd_proc;
		
	if (motor->p_kp_pos ==0){
		kp = conf_now->p_pid_kp;
		ki = conf_now->p_pid_ki;
		kd = conf_now->p_pid_kd;
		kd_proc = conf_now->p_pid_kd_proc;
	}
	else{
		kp = conf_now->p_pid_kp;
		ki = conf_now->p_pid_ki;
		kd = conf_now->p_pid_kd;
		kd_proc = conf_now->p_pid_kd_proc;

	}
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
	float p_term_pos;
	float d_term_pos;
	float pos_error;
	float d_term_proc_pos;

	//First treat position eror
	float pos_kp = motor->p_kp_pos;
	float pos_ki = motor->p_ki_pos;
	float pos_kd = 0.0f;// motor->p_kd_pos;
	float pos_kd_proc = conf_now->p_pid_kd_proc;



	//ADD GEAR CHANGE LOGIC HERE -> stop control and adapt setpoints to gear ratio for pos control
	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_SPEED) {
		motor->m_speed_i_term = 0.0f;
		motor->m_speed_prev_error = 0.0f;
		motor->m_speed_d_filter = 0.0f;
		motor->model_pos_i_term = 0.0f;
	//	motor->integrated_value= 0.0f;
	//	motor->unwrapped_theta=motor->model_pos_set_model;
	//	motor->unwrapped_theta_filtered=motor->model_pos_set_model;
	//	motor->model_pos_set_model=0.0f;
	//	motor->model_pos_d_filter=0.0f;

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

	float erpm = fabsf(motor->d_erpm_soll);

 	// piecewise -> then smoothstep each region
	float thr;

	if (erpm < 200.0f) {
		thr = 0.5f;
	} else if (erpm < 800.0f) {
		float x = (erpm - 200.0f) / 600.0f;          // 0..1
		x = x * x * (3.0f - 2.0f * x);
		thr = 0.5f + (3.0f - 0.5f) * x;
	} else if (erpm < 2000.0f) {
		float x = (erpm - 800.0f) / 1200.0f;         // 0..1
		x = x * x * (3.0f - 2.0f * x);
		thr = 3.0f + (10.0f - 3.0f) * x;
	} else {
		thr = 10.0f;
	}

	motor->rpm_inc_filter_th = thr;
	switch (conf_now->s_pid_speed_source) {
	case S_PID_SPEED_SRC_PLL:
		rpm = RADPS2RPM_f(motor->m_pll_speed);
		if (rpm  > motor->last_rpm +motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		}
		else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;
	case S_PID_SPEED_SRC_FAST:
		rpm = RADPS2RPM_f(motor->m_speed_est_fast);
		if (rpm  > motor->last_rpm +motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		}
		else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;
	case S_PID_SPEED_SRC_FASTER:
		rpm =RADPS2RPM_f(motor->m_speed_est_faster);
		if (rpm  > motor->last_rpm +motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		}
		else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;
	}

	motor->last_rpm = rpm;

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

	float speed_ratio = motor->d_erpm_soll / motor->p_speed_limit_pos_control_activation;
	float smooth_factor = tanhf(speed_ratio);

	// possibly changeable parameters
	float gear_ratio = motor->gear_ratio_bike;
	float incline    = 0.000f; // inclination angle, degrees
	float gearing    = (float)(motor->p_mech_gearing / gear_ratio);  // motor to wheel, /gearin. wheel to motor, *gearin, if speed
																	//if torque, motor to wheel * gearin, wheel to motor / gearin
	float slope      = incline * 3.14159265359f / 180.0f;            // radians

	float speed=rpm/(9.54929f*(motor->m_conf->si_motor_poles/2))*motor->p_wheel_radius/gearing; 
	
	
	// --- forces (same variable names, corrected formulas) ---
	float Area_s= motor->p_k_area*motor->p_height*motor->p_height;	//calculate section area
	float F_air     = 0.5f * motor->p_air_ro *motor->p_c_air*Area_s*speed*fabsf(speed);
	float F_roll    = smooth_factor*(motor->p_c_rr * motor->p_weight * 9.81f * cosf(slope));                    // Crr * m g cos(theta)
	float F_incline = motor->p_weight * 9.81f * sinf(slope);                                          // m g sin(theta)  (set incline=0 if you want it off)
	// simple viscous bearing drag in force-domain (N·s/m). keep name, fix units:
	float F_bearings = smooth_factor*(motor->p_c_bw * motor->p_k_v_bw) * speed;

	// F_res calculation	
	float F_combine = (F_air + F_roll + F_incline + F_bearings); // resistance force

	#define SCALE_INT 10000000.0f   // Float version for scaling

	motor->omega_fp = (int_fast64_t)(motor->m_speed_est_fast_corrected * SCALE_INT);

	UTILS_LP_FAST_I64(&motor->omega_filtered_fp, motor->omega_fp, 30);

	// --- Calculate motor torque ---
	motor->te_calculated = motor->m_motor_state.iq * motor->p_kT + (motor->m_motor_state.iq * motor->m_motor_state.id) * (motor->p_ld - motor->p_lq);
				
	


	motor->unwrapped_theta_filtered_prev = motor->unwrapped_theta_filtered;
	UTILS_LP_FAST(motor->unwrapped_theta_filtered,motor->unwrapped_theta,0.05f);
/*
	float inc_filter_th=motor->d_erpm_soll * (2.0f * M_PI / 60.0f) / motor->m_conf->si_motor_poles/2.0f +3.0f; // 5 times per dt

	if (motor->unwrapped_theta_filtered >= motor->unwrapped_theta_filtered_prev + inc_filter_th*dt){
		motor->unwrapped_theta_filtered = motor->unwrapped_theta_filtered_prev + inc_filter_th*dt;
	}
	else if (motor->unwrapped_theta_filtered <= motor->unwrapped_theta_filtered_prev - inc_filter_th*dt){
		motor->unwrapped_theta_filtered = motor->unwrapped_theta_filtered_prev - inc_filter_th*dt;
	}
*/	


	leso3_step(
    motor,
    dt,
    motor->te_calculated,          // motor torque [Nm] (applied/estimated)
    motor->unwrapped_theta_filtered);// unwrapped mechanical angle [rad]


	float error;


	motor->d_f_motor = (motor->te_calculated) / motor->p_wheel_radius * motor->p_mech_gearing;
	
	//motor->simulated_erpm = 360.0f*motor->erpm_time + sin(motor->erpm_time/0.5f)*600.0f;

	//if (motor->simulated_erpm >9000.0f){
	//	motor->simulated_erpm=9000.0f + sin(motor->erpm_time/0.5f)*800.0f;
	//}
	//motor->erpm_time+=dt;
	float wheel_erpm = motor->d_erpm_soll ;

	

	// Slip on motor side: how much faster the wheel is than the motor
	// (already both in RPM now, so direct compare)

	
	float slip_rpm = (wheel_erpm - rpm)/ (motor->m_conf->si_motor_poles / 2.0f);
	
	

	// Freewheeling logic
	// If the motor is being back-driven by the rider, and the rider torque exceeds the motor torque by a certain threshold,
	// then we disengage the motor (set iq to 0) and let it freewheel.
	// The motor will re-engage when the rider torque drops below the motor torque by a certain hysteresis threshold,
	// or if the rider applies a certain amount of positive torque (pedal push) to re-engage.
	// === FREEWHEEL: engage/disengage using RPMs ===
	// Wheel RPM (from virtual plant setpoint)

	// thresholds (tune or move to config)
	float FW_SLIP_ON_RPM,FW_SLIP_REENG,FW_T_DISENG,FW_T_REENG,FW_T_DISENG_FORCED; 
	FW_SLIP_ON_RPM =30.0f; // disengage if wheel outruns by >20 rpm
	FW_SLIP_REENG  = 20.0f;   // disengage if wheel outruns by >20 rpm
	FW_T_REENG      = 0.30f;   // Nm rider push to re-engage
	FW_T_DISENG = -0.20f;

	if (motor->freewheel_enabled || motor->forced_freewheel) {

		if (rpm < motor->p_speed_limit_pos_control_activation){
			FW_T_REENG=0.1f;
		}
		else if (motor->forced_freewheel)
		{
			FW_T_REENG=0.6f;
		}	
		else{
			FW_T_REENG=0.3f;
		}

		if (!motor->freewheel_active) {
			if ((motor->tp_observed <= FW_T_DISENG)||( motor->forced_freewheel)) {
				motor->freewheel_active = true;
			}
		} else {
			if (motor->tp_observed > FW_T_REENG && fabsf(slip_rpm)<FW_SLIP_REENG) {
				motor->freewheel_active = false;
				motor->forced_freewheel= false;
			}
		}
	}else
	{
		motor->freewheel_active = false;
		motor->forced_freewheel= false;
	}


	// --- Acceleration for integration (still uses full force model) ---
	motor->accel_ist = ((motor->Text_ext_hat_f/(motor->p_wheel_radius)*gearing - (F_combine) )  ) * SCALE_INT/motor->p_weight; // m/s^2 scaled

	
	// --- Position integration (fixed-point) ---
	motor->integrated_value = (int_fast64_t)(((motor->last_accel + motor->accel_ist) * (dt * SCALE_INT)) / (SCALE_INT * 2)) + motor->integrated_value;
	//if (motor->integrated_value < 0) {
    //motor->integrated_value = 0;
	//}
	motor->last_accel = motor->accel_ist;
	motor->d_speed_soll = speed;//(float)(motor->integrated_value / SCALE_INT);
	
	//speed setpoint calculation

	float v = (float)motor->integrated_value / SCALE_INT; // linear speed in m/s
	float omega_mech = v *  (gearing) / motor->p_wheel_radius;  // rad/s
	
	//should I also multiply by gear ratio here??? it seems no.

	motor->model_pos_set_model = motor->model_pos_set_model + omega_mech * dt; // rad, motor shaft
	motor->d_erpm_soll = omega_mech  * 9.54929 * (motor->m_conf->si_motor_poles/2.0f);

	pos_error=  motor->model_pos_set_model - motor->unwrapped_theta_filtered ; //desired position in rad, motor shaft
	p_term_pos = pos_error * pos_kp;
	motor->model_pos_i_term += dt * (pos_ki * pos_error - pos_ki * motor->model_pos_i_term);






	motor->model_pos_dt_int += dt;
	if (error == motor->model_pos_prev_error) {
		d_term_pos = 0.0;
	} else {
		d_term_pos = (pos_error - motor->model_pos_prev_error) * (pos_kd *pos_kd_proc/ motor->model_pos_dt_int);
		motor->m_pos_dt_int = 0.0;
	}
	// Filter D
	UTILS_LP_FAST(motor->model_pos_d_filter, d_term_pos, conf_now->p_pid_kd_filter);
	d_term_pos = motor->model_pos_d_filter;


	motor->model_pos_prev_error = pos_error;

	utils_truncate_number_abs((float*)&motor->model_pos_i_term, 1.0f - fabsf(p_term_pos));//windup protection
	//end pos control	

	float pos_output = p_term_pos + motor->model_pos_i_term ;
	utils_truncate_number(&pos_output, -1.0f, 1.0f);
	float speed_set_rpm = pos_output * conf_now->l_max_erpm;
	//end pos control

	//error = speed_set_rpm - rpm;

	error = motor->m_speed_pid_set_rpm - rpm;//motor->ekf_rpm; //original speed control without pos control and using kalman filter.

	p_term = error * conf_now->s_pid_kp * (1.0 / 20.0);
	d_term = (error - motor->m_speed_prev_error) * (conf_now->s_pid_kd / dt) * (1.0 / 20.0);

	// Filter D
	UTILS_LP_FAST(motor->m_speed_d_filter, d_term, conf_now->s_pid_kd_filter);
	d_term = motor->m_speed_d_filter;

	// Store previous error
	motor->m_speed_prev_error = error;

	// T_res calculation

	float T_res=F_combine*motor->p_wheel_radius/gearing;

	// ---------------- TORQUE FEEDFORWARD (ALL MOTOR DOMAIN, Method B, no prediction) ----------------

	float Te_ff = (-motor->Text_ext_hat_f);//- T_res);

	
	if (motor->d_erpm_soll<motor->p_speed_limit_pos_control_activation)
	{
		Te_ff = Te_ff*smooth_factor; //reduce feedforward torque at low speed to avoid sudden starts
	}
	// Store outputs
	motor->Te_set    = Te_ff;         // Te* (FF-only for now)
	motor->iq_set_ff = Te_ff / motor->p_kT;    // current FF [A]




	float iq_ff      = motor->iq_set_ff;  // from Te_ff above
	float iq_ff_norm = iq_ff / (conf_now->lo_current_max * conf_now->l_current_max_scale)*conf_now->p_pid_kd;


	motor->c_v_q_ff= iq_ff_norm*motor->m_res_est; //+ motor->m_speed_est_fast*lq*i_res_out; 

	// Other motor variables remain unchanged
	motor->d_speed = motor->m_speed_est_fast*motor->p_wheel_radius/(motor->m_conf->si_motor_poles)/gearing;//speed in m/s;
	motor->d_f_air = F_air;
	motor->d_f_combine = F_combine;
	motor->d_f_bearings = F_bearings;
	// For visibility (optional telemetry)
	motor->d_i_res = -iq_ff;          // store FF current (A)


	float pos_i_term;
	// Calculate output using smooth tanh linearization

	pos_i_term = motor->model_pos_i_term ;
	p_term_pos = p_term_pos;

	float output = p_term+ motor->m_speed_i_term  + d_term + iq_ff_norm + p_term_pos + pos_i_term;// + motor->model_pos_d_filter; 
	utils_truncate_number_abs(&output, 1.0);

	// Integrator windup protection
   // === FREEWHEEL: integrator handling
	float i_inc = error * conf_now->s_pid_ki * dt * (1.0f / 20.0f);
	bool wants_accel = (i_inc > 0.0f); // proxy for "controller wants to speed up"

	if (conf_now->s_pid_ki < 1e-9f) {
		motor->m_speed_i_term = 0.0f;
	} else {
		if (motor->freewheel_active && wants_accel) {
			// freeze or gently bleed the I-term while propulsion is blocked
			motor->m_speed_i_term *= 0.98f; // small decay; tune as needed
		} else {
			motor->m_speed_i_term += i_inc;
		}
		utils_truncate_number_abs(&motor->m_speed_i_term, 1.0f);
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
	
	if (motor->freewheel_active) {
    	output = 0.0f;  // no motor braking while freewheeling
		motor->c_v_q_ff=0.0f;
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


float smooth_force(float mag, float v, float v_eps) {
    return mag * tanhf(v / v_eps); // smoothly goes negative if v<0
}

// ===== Linear ESO (LESO) for external torque =============================
// States in m: leso_th [rad], leso_om [rad/s], leso_z [rad/s^2]
// Outputs: m->Text_ext_hat [Nm], m->Text_ext_hat_f [Nm]
// Inputs: theta_meas = UNWRAPPED angle [rad], Te_meas = kT * iq_meas_f [Nm], dt [s]

/*
inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas
){
	
    if (!m || dt <= 0.0f) return;
    const float J = m->p_J; if (!(J > 0.0f)) return;

    // ---------- Tuning (linear ESO) ----------
    const float fo_hz   = m->p_fo_hz;          // observer bandwidth (try 10–18 Hz)//100Hz //8Hz for small motor
    const float gz_hz   = m->p_gz_hz;           // tiny leak on z to suppress random-walk hiss (0–0.7 Hz)
    const float fc_TLPF = m->p_fc_TLPF;          // LPF for output (control/display)//500Hz //40 Hz for small motor
    // ----------------------------------------

    // Gains from bandwidth: standard cubic (s+wo)^3
    const float wo = 2.0f * (float)M_PI * fo_hz;
    const float b1 = 3.0f * wo;
    const float b2 = 3.0f * wo * wo;
    const float b3 =        wo * wo * wo;
    const float b0 = 1.0f / J;

    // Innovation (linear, no deadband, no e-LPF)
    const float e = theta_meas - m->leso_th;

    // LESO dynamics
    const float th_dot = m->leso_om + b1 * e;
    const float om_dot = b0 * Te_meas + m->leso_z + b2 * e;
    float       z_dot  = b3 * e - (2.0f*(float)M_PI*gz_hz) * m->leso_z; // small leak

	float Tdot_max = 250.0f; // Nm/s, max allowed rate of change of estimated torque
	// ------- RATE CLAMP on z (maps to motor torque via T = J*z) -------
	{
    // Max allowed change in z per step that corresponds to |dT/dt| <= Tdot_max
    const float dz_max = (Tdot_max / J) * dt;      // [1/s] * s = [1]
    float dz = dt * z_dot;

    // clamp dz to ±dz_max
    if (dz >  dz_max) dz =  dz_max;
    if (dz < -dz_max) dz = -dz_max;

    m->leso_z += dz;
	}
    // Integrate (Euler is fine at ≥1 kHz; use Tustin if you want extra smoothness)
    m->leso_th += dt * th_dot;
    m->leso_om += dt * om_dot;
    m->leso_z  += dt * z_dot;

    // External torque (motor side)
    const float Text_hat = J * m->leso_z;
    m->Text_ext_hat = Text_hat;

    // Output LPF (single-pole)
    const float aT = expf(-2.0f * (float)M_PI * fc_TLPF * dt);
    m->Text_ext_hat_f = aT * m->Text_ext_hat_f + (1.0f - aT) * Text_hat;

}*/
// Add this to your motor_all_state_t (or equivalent):
// float leso_th_ref;

// Add these to motor_all_state_t:
// float leso_th_ref;
// float leso_theta_prev;
// float leso_Te_prev;
/*
inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas
){
    if (!m || dt <= 0.0f) return;

    const float J = m->p_J;
    if (!(J > 0.0f)) return;

    // ---------- Tuning (linear ESO) ----------
    const float fo_hz   = m->p_fo_hz;
    const float gz_hz   = m->p_gz_hz;
    const float fc_TLPF = m->p_fc_TLPF;
    // ----------------------------------------

    const float B  = m->p_B;   // viscous [Nm/(rad/s)] (as you wrote)
    const float b0 = 1.0f / J;

    // Gains from bandwidth: (s + wo)^3
    const float wo = 2.0f * (float)M_PI * fo_hz;
    const float b1 = 3.0f * wo;
    const float b2 = 3.0f * wo * wo;
    const float b3 =        wo * wo * wo;

    const float gz = 2.0f * (float)M_PI * gz_hz;

    // --------- Reference-centering + init ----------
    if (!isfinite(m->leso_th_ref)) {
        m->leso_th_ref     = theta_meas;
        m->leso_th         = 0.0f;
        m->leso_om         = 0.0f;
        m->leso_z          = 0.0f;
        m->leso_theta_prev = theta_meas;   // init trapezoid memory
        m->leso_Te_prev    = Te_meas;
    }

    // Centered angle
    float theta_k   = theta_meas - m->leso_th_ref;
    float theta_km1 = m->leso_theta_prev - m->leso_th_ref;

    // Trapezoidal (midpoint) inputs for Tustin consistency
    const float theta_bar = 0.5f * (theta_k + theta_km1);
    const float Te_bar    = 0.5f * (Te_meas + m->leso_Te_prev);

    // Update memories for next call
    m->leso_theta_prev = theta_meas;
    m->leso_Te_prev    = Te_meas;

    const float h = 0.5f * dt;

    // M = I - hA
    const float m11 = 1.0f + h * b1;
    const float m12 = -h;
    const float m21 =  h * b2;
    const float m22 = 1.0f + h * (b0 * B);
    const float m23 = -h;
    const float m31 =  h * b3;
    const float m33 = 1.0f + h * gz;

    // N = I + hA
    const float n11 = 1.0f - h * b1;
    const float n12 =  h;
    const float n21 = -h * b2;
    const float n22 = 1.0f - h * (b0 * B);
    const float n23 =  h;
    const float n31 = -h * b3;
    const float n33 = 1.0f - h * gz;

    // affine input c (use trapezoidal averaged inputs)
    const float c1 = b1 * theta_bar;
    const float c2 = b2 * theta_bar + b0 * Te_bar;
    const float c3 = b3 * theta_bar;

    // RHS r = N*xk + dt*c
    const float thk = m->leso_th;
    const float omk = m->leso_om;
    const float zk  = m->leso_z;

    const float r1 = n11 * thk + n12 * omk + dt * c1;
    const float r2 = n21 * thk + n22 * omk + n23 * zk  + dt * c2;
    const float r3 = n31 * thk + n33 * zk  + dt * c3;

    // Solve M x_{k+1} = r exploiting sparsity
    const float inv_m11 = 1.0f / m11;
    const float inv_m33 = 1.0f / m33;

    const float th_const = r1 * inv_m11;
    const float th_om    = (-m12) * inv_m11;

    const float z_const  = (r3 - m31 * th_const) * inv_m33;
    const float z_om     = (-(m31 * th_om)) * inv_m33;

    const float const2 = m21 * th_const + m23 * z_const;
    const float coeff2 = m21 * th_om    + m22        + m23 * z_om;

    float om_new = omk;
    if (fabsf(coeff2) > 1e-12f) {
        om_new = (r2 - const2) / coeff2;
    }

    const float th_new = th_const + th_om * om_new;
    const float z_new  = z_const  + z_om  * om_new;

    // Commit
    m->leso_th = th_new;
    m->leso_om = om_new;
    m->leso_z  = z_new;

    // External torque (motor side)
    const float Text_hat = J * m->leso_z;
    m->Text_ext_hat = Text_hat;

    // Output LPF (single-pole)
    const float aT = expf(-2.0f * (float)M_PI * fc_TLPF * dt);
    m->Text_ext_hat_f = aT * m->Text_ext_hat_f + (1.0f - aT) * Text_hat;
}
*/
// Requires in state:
// m->leso_th, m->leso_om, m->leso_z

// Requires in state:
// m->leso_th, m->leso_om, m->leso_z
// Add to motor_all_state_t:
// float leso_th_prev;   // previous theta_meas [rad]
// float leso_om_fd;     // filtered derivative omega [rad/s]

// Add to motor_all_state_t:
// float leso_th_prev;   // previous theta_meas [rad]
// float leso_om_fd;     // filtered derivative omega [rad/s]

inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas
){
    if (!m || dt <= 0.0f) return;
    const float J = m->p_J; if (!(J > 0.0f)) return;

    // Tuning
    const float fo_hz   = m->p_fo_hz;
    const float gz_hz   = m->p_gz_hz;
    const float fc_TLPF = m->p_fc_TLPF;

    const float B  = m->p_B;            // Nm/(rad/s), can be 0
    const float b0 = 1.0f / J;

    const float wo = 2.0f * (float)M_PI * fo_hz;
    const float b1 = 3.0f * wo;
    const float b2 = 3.0f * wo * wo;
    const float b3 =        wo * wo * wo;
    const float gz = 2.0f * (float)M_PI * gz_hz;

    const float h = 0.5f * dt;

    // Shorthand old states
    const float thk = m->leso_th;
    const float omk = m->leso_om;
    const float zk  = m->leso_z;

    // e_k uses th_k
    const float ek = theta_meas - thk;

    // We solve a 3x3 linear system for [th_{k+1}, om_{k+1}, z_{k+1}]
    // derived from trapezoidal integration with e_{k+1} = theta_meas - th_{k+1}.
    //
    // The resulting linear system can be written as:
    // A * [th1, om1, z1]^T = rhs
    //
    // where A is sparse and depends only on (dt, b1,b2,b3,b0,B,gz).

    // Build A coefficients
    // From:
    // th1 = thk + h*( (omk + b1*ek) + (om1 + b1*(theta - th1)) )
    // om1 = omk + h*( (b0*(Te - B*omk) + zk + b2*ek) + (b0*(Te - B*om1) + z1 + b2*(theta - th1)) )
    // z1  = zk  + h*( (b3*ek - gz*zk) + (b3*(theta - th1) - gz*z1) )

    const float a11 = 1.0f + h * b1;          // th1 coeff in eq(th)
    const float a12 = -h;                      // om1 coeff in eq(th)
    const float a13 = 0.0f;

    const float a21 =  h * b2;                 // th1 coeff in eq(om)  (moved to LHS)
    const float a22 = 1.0f + h * (b0 * B);     // om1 coeff
    const float a23 = -h;                      // z1 coeff

    const float a31 =  h * b3;                 // th1 coeff in eq(z)
    const float a32 = 0.0f;
    const float a33 = 1.0f + h * gz;           // z1 coeff

    // RHS
    // eq(th): thk + h*(omk + b1*ek) + h*(b1*theta)  BUT theta term is included via (theta - th1) already on LHS
    // When expanding carefully, RHS becomes:
    const float rhs1 = thk + h * (omk + b1 * ek) + h * (b1 * theta_meas);

    // eq(om):
    const float fk_om = b0 * (Te_meas - B * omk) + zk + b2 * ek;
    const float rhs2  = omk + h * fk_om + h * (b0 * Te_meas + b2 * theta_meas);

    // eq(z):
    const float fk_z = b3 * ek - gz * zk;
    const float rhs3 = zk + h * fk_z + h * (b3 * theta_meas);

    // Solve sparse system like before:
    // a11*th1 + a12*om1          = rhs1
    // a21*th1 + a22*om1 + a23*z1 = rhs2
    // a31*th1          + a33*z1  = rhs3

    const float inv_a11 = 1.0f / a11;
    const float inv_a33 = 1.0f / a33;

    // th1(om1) = (rhs1 - a12*om1)/a11
    const float th_const = rhs1 * inv_a11;
    const float th_om    = (-a12) * inv_a11;

    // z1(om1) = (rhs3 - a31*th1)/a33
    const float z_const  = (rhs3 - a31 * th_const) * inv_a33;
    const float z_om     = (-(a31 * th_om)) * inv_a33;

    // Plug into eq2:
    // a21*th1 + a22*om1 + a23*z1 = rhs2
    const float const2 = a21 * th_const + a23 * z_const;
    const float coeff2 = a21 * th_om + a22 + a23 * z_om;

    float om1 = omk;
    if (fabsf(coeff2) > 1e-12f) {
        om1 = (rhs2 - const2) / coeff2;
    }

    const float th1 = th_const + th_om * om1;
    const float z1  = z_const  + z_om  * om1;

    // Commit
    m->leso_th = th1;
    m->leso_om = om1;
    m->leso_z  = z1;

    // External torque (motor side)
    const float Text_hat = J * m->leso_z;
    m->Text_ext_hat = Text_hat;

    // LPF for output torque
    const float aT = expf(-2.0f * (float)M_PI * fc_TLPF * dt);
    m->Text_ext_hat_f = aT * m->Text_ext_hat_f + (1.0f - aT) * Text_hat;
}
