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

	// First treat position error
	float pos_kp = motor->p_kp_pos;
	float pos_ki = motor->p_ki_pos;
	float pos_kd = 0.0f; // motor->p_kd_pos;
	float pos_kd_proc = conf_now->p_pid_kd_proc;

	// PID is off. Return.
	if (motor->m_control_mode != CONTROL_MODE_SPEED) {
		motor->m_speed_i_term = 0.0f;
		motor->m_speed_prev_error = 0.0f;
		motor->m_speed_d_filter = 0.0f;
		motor->model_pos_prev_error = 0.0f;
		motor->model_pos_d_filter = 0.0f;
		motor->model_pos_i_term = 0.0f;
		motor->Text_ext_hat_f = 0.0f;
		return;
	}

	if (conf_now->s_pid_ramp_erpms_s > 0.0) {
		utils_step_towards((float*)&motor->m_speed_pid_set_rpm, motor->m_speed_command_rpm,
				conf_now->s_pid_ramp_erpms_s * dt);
		if (!index_found) {
			utils_truncate_number_abs(&motor->m_speed_pid_set_rpm, conf_now->foc_openloop_rpm);
		}
		utils_truncate_number(&motor->m_speed_pid_set_rpm, conf_now->l_min_erpm, conf_now->l_max_erpm);
	}

	float rpm = 0.0f;

	// ----- Speed measurement slew limiter schedule (based on setpoint magnitude) -----
	float erpm = fabsf(motor->d_erpm_soll);

	float thr;
	if (erpm < 500.0f) {
		thr = 0.5f;
	} else if (erpm < 800.0f) {
		float x = (erpm - 500.0f) / 300.0f; // 0..1
		x = x * x * (3.0f - 2.0f * x);
		thr = 0.5f + (3.0f - 0.5f) * x;
	} else if (erpm < 2000.0f) {
		float x = (erpm - 800.0f) / 1200.0f; // 0..1
		x = x * x * (3.0f - 2.0f * x);
		thr = 3.0f + (10.0f - 3.0f) * x;
	} else {
		thr = 10.0f;
	}

	motor->rpm_inc_filter_th = thr;

	switch (conf_now->s_pid_speed_source) {
	case S_PID_SPEED_SRC_PLL:
		rpm = RADPS2RPM_f(motor->m_pll_speed);
		if (rpm > motor->last_rpm + motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		} else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;

	case S_PID_SPEED_SRC_FAST:
		rpm = RADPS2RPM_f(motor->m_speed_est_fast);
		if (rpm > motor->last_rpm + motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		} else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;

	case S_PID_SPEED_SRC_FASTER:
		rpm = RADPS2RPM_f(motor->m_speed_est_faster);
		if (rpm > motor->last_rpm + motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm + motor->rpm_inc_filter_th;
		} else if (rpm < motor->last_rpm - motor->rpm_inc_filter_th) {
			rpm = motor->last_rpm - motor->rpm_inc_filter_th;
		}
		break;
	}

	motor->last_rpm = rpm;

	// ================= Gain scheduling (bidirectional-safe) =================
	// Floors are ZERO where requested; saturation point is 2x activation threshold.
	const float erpm_act = motor->p_speed_limit_pos_control_activation*4; // ~500x4.0
	const float erpm_sat = 2.0f * erpm_act;                              // ~1000x4.0

	// Use magnitude for scheduling only; keep signed signals for control.
	float erpm_abs_meas = 0.5f * (fabsf(rpm) + fabsf(motor->d_erpm_soll));


	// g = 0 at 0, g = 1 at erpm_sat
	float g = (erpm_sat > 1e-6f) ? (erpm_abs_meas / erpm_sat) : 1.0f;
	if (g < 0.0f) g = 0.0f;
	if (g > 1.0f) g = 1.0f;
	g = g * g * (3.0f - 2.0f * g); // smoothstep
	g = tanhf(2.0f * g);

	// gpos = 0 at <=5 ERPM, gpos = 1 at >=erpm_sat
	float gpos;
	if (erpm_abs_meas <= 5.0f) {
		gpos = 0.0f;
	} else {
		const float denom = (erpm_sat - 5.0f);
		gpos = (denom > 1e-6f) ? ((erpm_abs_meas - 5.0f) / denom) : 1.0f;
		if (gpos < 0.0f) gpos = 0.0f;
		if (gpos > 1.0f) gpos = 1.0f;
		gpos = gpos * gpos * (3.0f - 2.0f * gpos);
		gpos = tanhf(2.0f * gpos);
	}

	// ki_pos_scale: increases at half speed (saturates at erpm_act instead of erpm_sat)
	float gi_pos;
	if (erpm_abs_meas <= 5.0f) {
		gi_pos = 0.0f;
	} else {
		const float denom_i = (erpm_act - 5.0f);
		gi_pos = (denom_i > 1e-6f) ? ((erpm_abs_meas - 5.0f) / denom_i) : 1.0f;
		if (gi_pos < 0.0f) gi_pos = 0.0f;
		if (gi_pos > 1.0f) gi_pos = 1.0f;
		gi_pos = tanhf(2.0f * gi_pos);
	}

	// Speed: Kp 0.3 -> 1.0 over 0..erpm_sat; Ki 0 -> 1 over 0..erpm_sat
	const float kp_speed_scale = 0.3f + 0.7f * g;
	const float ki_speed_scale = gi_pos;

	// Position: Kp 0.1 -> 1.0 over 5..erpm_sat; Ki 0 -> 1 over 5..erpm_act (saturates faster)
	const float kp_pos_scale   = 0.1f + 0.9f * gpos;
	const float ki_pos_scale   = gi_pos;

	const float pos_kp_eff = pos_kp * kp_pos_scale;
	const float pos_ki_eff = pos_ki * ki_pos_scale;

	const float sp_kp_eff  = conf_now->s_pid_kp * kp_speed_scale;
	const float sp_ki_eff  = conf_now->s_pid_ki * ki_speed_scale;

	// Keep I-term from "remembering" nonsense at standstill (symmetric)
	if (erpm_abs_meas < 2.0f) {
		motor->m_speed_i_term = 0.0f;
	}

	// ---------------- Encoder unwrap ----------------
	float angle_deg_now = encoder_read_deg();
	float angle_rad_now = angle_deg_now * (M_PI / 180.0f);

	float delta_rad = utils_angle_difference_rad(angle_rad_now, motor->kalman_last_angle_rad);
	motor->kalman_last_angle_rad = angle_rad_now;

	motor->unwrapped_theta += delta_rad;
	motor->kalman_last_delta_rad = delta_rad;

	float speed_ratio = motor->d_erpm_soll / motor->p_speed_limit_pos_control_activation;
	float smooth_factor = tanhf(speed_ratio);

	// Plant parameters
	float gear_ratio = motor->gear_ratio_bike;
	float incline = 0.000f;
	float gearing = (float)(motor->p_mech_gearing / gear_ratio);
	float slope = incline * 3.14159265359f / 180.0f;

	float speed = rpm / (9.54929f * (motor->m_conf->si_motor_poles / 2)) * motor->p_wheel_radius / gearing;

	// Forces
	float Area_s = motor->p_k_area * motor->p_height * motor->p_height;
	float F_air = 0.5f * motor->p_air_ro * motor->p_c_air * Area_s * speed * fabsf(speed);
	float F_roll = smooth_factor * (motor->p_c_rr * motor->p_weight * 9.81f * cosf(slope));
	float F_incline = motor->p_weight * 9.81f * sinf(slope);
	float F_bearings = smooth_factor * (motor->p_c_bw * motor->p_k_v_bw) * speed;

	float F_combine = (F_air + F_roll + F_incline + F_bearings);

	#define SCALE_INT 10000000.0f

	// Motor torque estimate
	motor->te_calculated =
		motor->m_motor_state.iq * motor->p_kT +
		(motor->m_motor_state.iq * motor->m_motor_state.id) * (motor->p_ld - motor->p_lq);

	motor->unwrapped_theta_filtered_prev = motor->unwrapped_theta_filtered;
	UTILS_LP_FAST(motor->unwrapped_theta_filtered, motor->unwrapped_theta, 0.05f);
	float omega_for_leso= rpm / (9.54929f * (motor->m_conf->si_motor_poles / 2));

	leso3_step(
		motor,
		dt,
		motor->te_calculated,
		motor->unwrapped_theta_filtered,
		omega_for_leso
	);

	float error;

	motor->d_f_motor = (motor->te_calculated) / motor->p_wheel_radius * motor->p_mech_gearing;

	float wheel_erpm = motor->d_erpm_soll;
	float slip_rpm = (wheel_erpm - rpm) / (motor->m_conf->si_motor_poles / 2.0f);

	{
		float FW_SLIP_ON_RPM, FW_SLIP_REENG, FW_T_DISENG, FW_T_REENG, FW_T_DISENG_FORCED;
		FW_SLIP_ON_RPM = 30.0f;
		FW_SLIP_REENG  = 20.0f;
		FW_T_REENG     = 0.30f;
		FW_T_DISENG    = -0.20f;

		if (motor->freewheel_enabled || motor->forced_freewheel) {

			if (rpm < motor->p_speed_limit_pos_control_activation) {
				FW_T_REENG = 0.1f;
			} else if (motor->forced_freewheel) {
				FW_T_REENG = 0.6f;
			} else {
				FW_T_REENG = 0.3f;
			}

			if (!motor->freewheel_active) {
				if ((motor->tp_observed <= FW_T_DISENG) || (motor->forced_freewheel)) {
					motor->freewheel_active = true;
				}
			} else {
				if (motor->tp_observed > FW_T_REENG && fabsf(slip_rpm) < FW_SLIP_REENG) {
					motor->freewheel_active = false;
					motor->forced_freewheel = false;
				}
			}
		} else {
			motor->freewheel_active = false;
			motor->forced_freewheel = false;
		}
	}

	// Acceleration for integration
	motor->accel_ist = ((motor->Text_ext_hat_f / (motor->p_wheel_radius) * gearing - (F_combine))) *
	                   SCALE_INT / motor->p_weight;

	// Position integration (fixed-point)
	motor->integrated_value =
		(int_fast64_t)(((motor->last_accel + motor->accel_ist) * (dt * SCALE_INT)) / (SCALE_INT * 2)) +
		motor->integrated_value;

	motor->last_accel = motor->accel_ist;
	motor->d_speed_soll = speed;

	// Speed setpoint calculation
	float v = (float)motor->integrated_value / SCALE_INT;
	float omega_mech = v * (gearing) / motor->p_wheel_radius;

	motor->model_pos_set_model = motor->model_pos_set_model + omega_mech * dt;
	motor->d_erpm_soll = omega_mech * 9.54929f * (motor->m_conf->si_motor_poles / 2.0f);

	pos_error = motor->model_pos_set_model - motor->unwrapped_theta_filtered;

	// Position P/I with scheduled gains
	p_term_pos = pos_error * pos_kp_eff;
	motor->model_pos_i_term += dt * (pos_ki_eff * pos_error - pos_ki_eff * motor->model_pos_i_term);

	// Position D (pos_kd == 0 currently)
	motor->model_pos_dt_int += dt;
	if (pos_error == motor->model_pos_prev_error) {
		d_term_pos = 0.0f;
	} else {
		d_term_pos = (pos_error - motor->model_pos_prev_error) * (pos_kd * pos_kd_proc / motor->model_pos_dt_int);
		motor->m_pos_dt_int = 0.0f;
	}

	UTILS_LP_FAST(motor->model_pos_d_filter, d_term_pos, conf_now->p_pid_kd_filter);
	d_term_pos = motor->model_pos_d_filter;

	motor->model_pos_prev_error = pos_error;

	utils_truncate_number_abs((float*)&motor->model_pos_i_term, 1.0f - fabsf(p_term_pos));

	float pos_output = p_term_pos + motor->model_pos_i_term;
	utils_truncate_number(&pos_output, -1.0f, 1.0f);

	// Speed loop error (signed)
	error = motor->m_speed_pid_set_rpm - rpm;

	// Speed P/D with scheduled Kp (Ki is applied to i_inc below)
	p_term = error * sp_kp_eff * (1.0f / 20.0f);
	d_term = (error - motor->m_speed_prev_error) * (conf_now->s_pid_kd / dt) * (1.0f / 20.0f);

	UTILS_LP_FAST(motor->m_speed_d_filter, d_term, conf_now->s_pid_kd_filter);
	d_term = motor->m_speed_d_filter;

	motor->m_speed_prev_error = error;

	// Resistance torque
	float T_res = F_combine * motor->p_wheel_radius / gearing;

	// Torque feedforward
	float Te_ff = (-motor->Text_ext_hat_f - T_res);

	if (motor->d_erpm_soll < motor->p_speed_limit_pos_control_activation) {
		Te_ff = Te_ff * smooth_factor;
	}

	motor->Te_set = Te_ff;
	motor->iq_set_ff = Te_ff / motor->p_kT;

	float iq_ff = motor->iq_set_ff;

	// Keep your original scaling (note: you might want to remove *p_pid_kd later)
	float iq_ff_norm = iq_ff / (conf_now->lo_current_max * conf_now->l_current_max_scale) * conf_now->p_pid_kd;

	motor->c_v_q_ff = iq_ff_norm * motor->m_res_est;

	motor->d_speed = motor->m_speed_est_fast * motor->p_wheel_radius / (motor->m_conf->si_motor_poles) / gearing;
	motor->d_f_air = F_air;
	motor->d_f_combine = F_combine;
	motor->d_f_bearings = F_bearings;
	motor->d_i_res = -iq_ff;

	float pos_i_term = motor->model_pos_i_term;

	// Controller output (your original summation)
	float output = p_term + motor->m_speed_i_term + d_term + iq_ff_norm + p_term_pos + pos_i_term;
	utils_truncate_number_abs(&output, 1.0f);

	// Integrator windup protection / update (scheduled Ki)
	float i_inc = error * sp_ki_eff * dt * (1.0f / 20.0f);
	bool wants_accel = (i_inc > 0.0f);

	if (sp_ki_eff < 1e-9f) {
		motor->m_speed_i_term = 0.0f;
	} else {
		if (motor->freewheel_active && wants_accel) {
			motor->m_speed_i_term *= 0.98f;
		} else {
			motor->m_speed_i_term += i_inc;
		}
		utils_truncate_number_abs(&motor->m_speed_i_term, 1.0f);
	}

	// Optionally disable braking
	if (!conf_now->s_pid_allow_braking) {
		if (rpm > 20.0f && output < 0.0f) {
			output = 0.0f;
		}

		if (rpm < -20.0f && output > 0.0f) {
			output = 0.0f;
		}
	}

	if (motor->freewheel_active) {
		output = 0.0f;
		motor->c_v_q_ff = 0.0f;
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

inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
    float omega_meas   // mechanical rad/s (measured/derived)
){
    if (!m || dt <= 0.0f) return;
    const float J = m->p_J; if (!(J > 0.0f)) return;

    // ---------- Tuning ----------
    const float fo_hz   = m->p_fo_hz;
    const float gz_hz   = m->p_gz_hz;
    const float fc_TLPF = m->p_fc_TLPF;
    // ----------------------------

    const float B  = m->p_B;            // Nm/(rad/s), can be 0
    const float b0 = 1.0f / J;

    // ---------- Coulomb friction model (NEW) ----------
    // Tc: Coulomb friction amplitude [Nm]
    // ws: smoothing speed [rad/s] to avoid sign chatter around zero
    const float Tc = m->p_Tc;           // <-- add this parameter to motor_all_state_t
    float ws = m->p_Tc_ws;              // <-- add this parameter too (rad/s)
    if (!(ws > 1e-6f)) ws = 1.0f;       // fallback (tune, e.g. 1..5 rad/s)

    // Use omega_meas for a stable sign near zero (less self-excited flip-flop)
    const float Tc_term = Tc * tanhf(omega_meas / ws);

    // Effective electromagnetic torque used in the model
    const float Te_eff = Te_meas - Tc_term;
    // -------------------------------------------------

    const float wo = 2.0f * (float)M_PI * fo_hz;
    const float b1 = 3.0f * wo;
    const float b2 = 3.0f * wo * wo;
    const float b3 =        wo * wo * wo;
    const float gz = 2.0f * (float)M_PI * gz_hz;

    const float h = 0.5f * dt;

    // Old states
    const float thk = m->leso_th;
    const float omk = m->leso_om;
    const float zk  = m->leso_z;

    const float ek = theta_meas - thk;

    // --- Linear system coefficients ---
    const float a11 = 1.0f + h * b1;
    const float a12 = -h;

    const float a21 =  h * b2;
    const float a22 = 1.0f + h * (b0 * B);
    const float a23 = -h;

    const float a31 =  h * b3;
    const float a33 = 1.0f + h * gz;

    const float rhs1 = thk + h * (omk + b1 * ek) + h * (b1 * theta_meas);

    // Use Te_eff (NEW)
    const float fk_om = b0 * (Te_eff - B * omk) + zk + b2 * ek;
    const float rhs2  = omk + h * fk_om + h * (b0 * Te_eff + b2 * theta_meas);

    const float fk_z  = b3 * ek - gz * zk;
    const float rhs3  = zk + h * fk_z + h * (b3 * theta_meas);

    // Solve sparse system
    const float inv_a11 = 1.0f / a11;
    const float inv_a33 = 1.0f / a33;

    const float th_const = rhs1 * inv_a11;
    const float th_om    = (-a12) * inv_a11;

    const float z_const  = (rhs3 - a31 * th_const) * inv_a33;
    const float z_om     = (-(a31 * th_om)) * inv_a33;

    const float const2 = a21 * th_const + a23 * z_const;
    const float coeff2 = a21 * th_om + a22 + a23 * z_om;

    float om1 = omk;
    if (fabsf(coeff2) > 1e-12f) {
        om1 = (rhs2 - const2) / coeff2;
    }

    float th1 = th_const + th_om * om1;
    float z1  = z_const  + z_om  * om1;

    // ====================== CLAMPS (post-solve) ======================

    // 1) Physical torque clamp -> clamp z magnitude
    float Te_max = 50.0f; // Nm fallback
    if (m->m_conf) {
        const float Imax = m->m_conf->lo_current_max * m->m_conf->l_current_max_scale; // A
        Te_max = fabsf(Imax) * fabsf(m->p_kT) * 1.2f + 0.5f; // Nm
    }
    const float z_abs_max = Te_max / J; // since Text_hat = J*z
    if (z1 >  z_abs_max) z1 =  z_abs_max;
    if (z1 < -z_abs_max) z1 = -z_abs_max;

    // 2) z increment clamp (limits Text_dot)
    const float Tdot_max = 4000.0f;                 // Nm/s (tune)
    const float dz_max   = (Tdot_max / J) * dt;
    {
        float dz = z1 - zk;
        if (dz >  dz_max) dz =  dz_max;
        if (dz < -dz_max) dz = -dz_max;
        z1 = zk + dz;
    }

    // 3) ω plausibility clamp around measured mechanical omega
    {
        const float om_meas = omega_meas;    // mechanical rad/s
        const float om_floor = 10.0f;        // rad/s floor (~95 rpm). Tune 5..20.
        float dom_allow = 0.20f * fabsf(om_meas) + om_floor;

        float om_abs_max = 500.0f; // rad/s fallback (~4775 rpm mech)
        if (m->m_conf) {
            const float pp = 0.5f * m->m_conf->si_motor_poles;
            const float erpm_max = (float)m->m_conf->l_max_erpm;
            const float om_mech_max = (erpm_max / pp) * (2.0f * (float)M_PI / 60.0f);
            om_abs_max = 1.2f * om_mech_max;
        }

        if (dom_allow > om_abs_max) dom_allow = om_abs_max;

        const float om_lo = om_meas - dom_allow;
        const float om_hi = om_meas + dom_allow;

        float om_target = om1;
        if (om_target < om_lo) om_target = om_lo;
        if (om_target > om_hi) om_target = om_hi;

        const float dom_rate = 550.0f; // rad/s^2 (tune)
        float dom = om_target - omk;
        const float dom_max = dom_rate * dt;
        if (dom >  dom_max) dom =  dom_max;
        if (dom < -dom_max) dom = -dom_max;
        om1 = omk + dom;

        if (om1 >  om_abs_max) om1 =  om_abs_max;
        if (om1 < -om_abs_max) om1 = -om_abs_max;
    }

    // ==================== Commit ====================
    m->leso_th = th1;
    m->leso_om = om1;
    m->leso_z  = z1;

    // External torque (motor side)
    float Text_hat = J * z1;
    if (Text_hat >  Te_max) Text_hat =  Te_max;
    if (Text_hat < -Te_max) Text_hat = -Te_max;

    m->Text_ext_hat = Text_hat;

    // LPF for output torque
    const float aT = expf(-2.0f * (float)M_PI * fc_TLPF * dt);
    m->Text_ext_hat_f = aT * m->Text_ext_hat_f + (1.0f - aT) * Text_hat;
}


static inline float falf(float e, float alpha, float delta) {
    float ae = fabsf(e);
    if (ae <= delta) {
        // Linear in the boundary layer, scaled for continuity
        return e / powf(delta, 1.0f - alpha);
    } else {
        return copysignf(powf(ae, alpha), e);
    }
}
