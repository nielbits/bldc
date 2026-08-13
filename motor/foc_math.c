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

/*
 * Bicycle torque-observer implementation:
 *
 * The original third-order linear extended-state observer (LESO) is active.
 * The simplified speed-differentiation estimator is retained under #if 0
 * for later diagnostic comparison. No state-structure changes are required.
 */

/*
 * Coordinate convention for the bicycle model.
 *
 * On this machine:
 *
 *     native VESC positive rotation = physical bicycle reverse
 *
 * The bicycle model uses:
 *
 *     bicycle positive rotation = physical bicycle forward
 *
 * Therefore the conversion sign is -1.
 *
 * Since the conversion is only a sign reversal, the inverse conversion is
 * identical:
 *
 *     bicycle_value = BIKE_DIRECTION_SIGN * vesc_value
 *     vesc_value    = BIKE_DIRECTION_SIGN * bicycle_value
 *
 * Keep foc_encoder_inverted at the value required for correct commutation.
 * This high-level sign is independent of motor-commutation configuration.
 */
#define BIKE_DIRECTION_SIGN (-1.0f)

static inline float bike_from_vesc(float value) {
    return BIKE_DIRECTION_SIGN * value;
}

static inline float vesc_from_bike(float value) {
    return BIKE_DIRECTION_SIGN * value;
}

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
    float pos_error;

    /*
     * Do not require a physical incremental-encoder index pulse for the
     * custom bicycle-controller state machine.
     *
     * The configured encoder is absolute and already supplies a usable
     * mechanical angle. encoder_index_found() can therefore remain false
     * permanently even while encoder_read_deg() and FOC commutation work.
     *
     * Preserve every available readiness indication:
     *   1. the value supplied by the caller,
     *   2. the native VESC index-found flag,
     *   3. a finite absolute-encoder angle.
     */
    const bool index_found_from_caller = index_found;
    const bool native_index_found = encoder_index_found();
    bool encoder_angle_valid = false;

    /*
     * Keep mechanical position tracking alive independently of whether
     * the speed controller is enabled.
     */
    {
        const float angle_deg_now = encoder_read_deg();

        if (isfinite(angle_deg_now)) {
            encoder_angle_valid = true;

            const float angle_rad_now =
                angle_deg_now * ((float)M_PI / 180.0f);

            float delta_rad =
                utils_angle_difference_rad(
                    angle_rad_now,
                    motor->kalman_last_angle_rad
                );

            motor->kalman_last_angle_rad = angle_rad_now;

            if (conf_now->foc_encoder_inverted) {
                delta_rad = -delta_rad;
            }

            /*
             * encoder_read_deg() is first corrected into the native VESC
             * mechanical direction above. Convert that increment into the
             * bicycle convention, where physical forward is positive.
             */
            delta_rad = bike_from_vesc(delta_rad);

            motor->kalman_last_delta_rad = delta_rad;
            motor->unwrapped_theta += delta_rad;
        }
    }

    const bool encoder_ready =
        index_found_from_caller ||
        native_index_found ||
        encoder_angle_valid;

    float pos_kp = conf_now->p_pid_kp;
    float pos_ki = conf_now->p_pid_ki;

    bool ctrl_enabled = false;

    if (motor->m_control_mode != CONTROL_MODE_SPEED) {
        motor->m_speed_i_term = 0.0f;
        motor->m_speed_prev_error = 0.0f;
        motor->m_speed_d_filter = 0.0f;
        motor->model_pos_prev_error = 0.0f;
        motor->model_pos_d_filter = 0.0f;
        motor->model_pos_i_term = 0.0f;

        /*
         * Reset the LESO safely while keeping its position and speed states
         * synchronized with the current bicycle-coordinate measurements.
         *
         * leso_z4 is not part of the third-order LESO. It is reused only as
         * an initialization flag:
         *
         *     0.0f = initialize on the next valid LESO call
         *     1.0f = initialized
         */
        motor->leso_th = motor->unwrapped_theta;
        motor->leso_om =
            bike_from_vesc(
                motor->m_pll_speed *
                motor->c_inv_pole_pairs
            );
        motor->leso_z = 0.0f;
        motor->leso_z4 = 0.0f;

        motor->Tf_hat = 0.0f;
        motor->Text_ext_hat = 0.0f;
        motor->Text_ext_hat_f = 0.0f;
        motor->tp_observed = 0.0f;

        motor->freewheel_active = false;
        motor->forced_freewheel = false;

        /*
         * Do NOT force ctrl_sm_state back to START merely because the VESC
         * control mode is temporarily not CONTROL_MODE_SPEED.
         *
         * A short mode transition would otherwise destroy the controller
         * enable state and force the full START -> INDEX_FOUND -> ENABLE
         * sequence again.
         *
         * The state machine is allowed to return to START only when all encoder
         * readiness indications are lost in the state-machine code below.
         */

        motor->status_bits = 0;
        motor->status_bits |=
            (motor->forced_freewheel) << STATUS_BIT_FORCED_FREEWHEEL;
        motor->status_bits |=
            (motor->ctrl_sm_state == CTRL_SM_START)
            << STATUS_BIT_CTRL_SM_START;
        motor->status_bits |=
            (motor->ctrl_sm_state == CTRL_SM_INDEX_FOUND)
            << STATUS_BIT_CTRL_SM_INDEX_FOUND;
        motor->status_bits |=
            (motor->ctrl_sm_state == CTRL_SM_ENABLE)
            << STATUS_BIT_CTRL_SM_ENABLE;

        return;
    }

    const float radps_to_rpm       = motor->c_radps_to_rpm;
    const float mech_radps_to_erpm = motor->c_mech_radps_to_erpm;
    const float erpm_to_mech_radps = motor->c_erpm_to_mech_radps;
    const float inv_pole_pairs     = motor->c_inv_pole_pairs;

    float rpm = 0.0f;
    float t_rpm = 0.0f;

    const float clamp = 30000.0f;

    switch (conf_now->s_pid_speed_source) {
    case S_PID_SPEED_SRC_PLL:
        t_rpm = motor->m_pll_speed * radps_to_rpm;
        break;
    case S_PID_SPEED_SRC_FAST:
        t_rpm = motor->m_speed_est_fast * radps_to_rpm;
        break;
    case S_PID_SPEED_SRC_FASTER:
        t_rpm = motor->m_speed_est_faster * radps_to_rpm;
        break;
    }

    if (fabsf(t_rpm) > clamp) {
        /*
         * last_rpm is stored in bicycle coordinates.
         */
        rpm = motor->last_rpm;
    } else {
        /*
         * t_rpm is native VESC ERPM. Convert it once here. From this point
         * onward, every speed, position, model and controller quantity in
         * this function uses bicycle coordinates.
         */
        rpm = bike_from_vesc(t_rpm);
    }

    motor->last_rpm = rpm;
    ctrl_sm_state_t prev = motor->ctrl_sm_state;

    // ================= Controller enable state machine =================
    {
        const float rpm_deadband = 5.0f;
        const uint32_t still_req_cycles = 5000;

        const float abs_rpm = fabsf(motor->m_pll_speed * radps_to_rpm);

        switch (motor->ctrl_sm_state) {
        default:
        case CTRL_SM_START:
            motor->ctrl_sm_still_cycles = 0;

            if (encoder_ready) {
                motor->ctrl_sm_state = CTRL_SM_INDEX_FOUND;
                //motor->unwrapped_theta = 0.0f;
            }
            break;

        case CTRL_SM_INDEX_FOUND:
            if (!encoder_ready) {
                motor->ctrl_sm_state = CTRL_SM_START;
                motor->ctrl_sm_still_cycles = 0;
                motor->leso_z4 = 0.0f;
                break;
            }

            if (abs_rpm <= rpm_deadband) {
                if (motor->ctrl_sm_still_cycles < still_req_cycles) {
                    motor->ctrl_sm_still_cycles++;
                }
                if (motor->ctrl_sm_still_cycles >= still_req_cycles) {
                    motor->ctrl_sm_state = CTRL_SM_ENABLE;
                }
            } else {
                motor->ctrl_sm_still_cycles = 0;
            }
            break;

        case CTRL_SM_ENABLE:
            if (!encoder_ready) {
                motor->ctrl_sm_state = CTRL_SM_START;
                motor->ctrl_sm_still_cycles = 0;
                motor->leso_z4 = 0.0f;
            }
            break;
        }

        ctrl_enabled = (motor->ctrl_sm_state == CTRL_SM_ENABLE);

        bool entered_enable = (prev != CTRL_SM_ENABLE) && (motor->ctrl_sm_state == CTRL_SM_ENABLE);
     
        if (entered_enable) {
            motor->model_v = 0.0f;
            motor->model_accel_prev = 0.0f;

            motor->model_pos_set_model = motor->unwrapped_theta_filtered;

            motor->model_pos_i_term = 0.0f;
            motor->m_speed_i_term = 0.0f;
            motor->m_speed_prev_error = 0.0f;
            motor->m_speed_d_filter = 0.0f;
        }
    }

   // ================= Gain scheduling =================
    const float erpm_abs_meas = 0.5f * (fabsf(rpm) + fabsf(motor->d_erpm_soll));

    // User parameters
    const float spd_floor    = motor->p_sched_spd_floor;
    const float pos_floor    = motor->p_sched_pos_floor;
    const float pos_dead     = motor->p_sched_pos_dead_erpm;
    const float spd_sat_erpm = motor->p_sched_spd_sat_erpm;
    const float pos_sat_erpm = motor->p_sched_pos_sat_erpm;


    // Speed schedule
    const float m_spd = ramp_rational_p2(erpm_abs_meas, spd_sat_erpm);
    const float g_spd = map_floor_local(m_spd, spd_floor);

    // Position schedule
    float x_pos = erpm_abs_meas - pos_dead;
    if (x_pos < 0.0f) {
        x_pos = 0.0f;
    }

    const float m_pos = ramp_rational_p2(x_pos, pos_sat_erpm);
    const float g_pos = map_floor_local(m_pos, pos_floor);

    // Effective gains
    const float pos_kp_eff = pos_kp * g_pos * 100.0f;
    const float pos_ki_eff = pos_ki * g_pos * 100.0f;

    const float sp_kp_eff  = conf_now->s_pid_kp * g_spd;
    const float sp_ki_eff  = conf_now->s_pid_ki * g_spd;

    if (erpm_abs_meas < 2.0f) {
        motor->m_speed_i_term = 0.0f;
    }

    /*
     * Encoder unwrapping is performed at the beginning of this function,
     * before the control-mode early return.
     */

    // Plant parameters
    UTILS_LP_FAST(motor->p_gear_ratio_filtered, motor->gear_ratio_bike, 0.001f);
    float gear_ratio = motor->p_gear_ratio_filtered;
    float gearing = (float)(motor->p_mech_gearing / gear_ratio);
    UTILS_LP_FAST(motor->p_incline_filtered, motor->p_incline_deg, 0.01f);

    float incline_base_deg = motor->p_incline_filtered;
    float incline_deg_cmd = incline_base_deg;

    if (motor->pumptrack_enabled) {
        float T_sec = motor->pumptrack_period_min * 60.0f;
        if (T_sec < 1e-4f) T_sec = 1e-4f;

        motor->pumptrack_time += dt;
        if (motor->pumptrack_time >= T_sec) {
            motor->pumptrack_time -= T_sec * floorf(motor->pumptrack_time / T_sec);
        }

        float phase = 2.0f * (float)M_PI * motor->pumptrack_time / T_sec - 0.5f * (float)M_PI;
        incline_deg_cmd = 0.5f * incline_base_deg * (1.0f + sinf(phase));
    } else {
        incline_deg_cmd = incline_base_deg;
        motor->pumptrack_time = 0.0f;
    }

    motor->incline_result = incline_deg_cmd;

    /*
     * Sign-audited virtual-bicycle force model.
     *
     * IMPORTANT:
     * All forces that update model_v must be computed from model_v, not from
     * the measured motor speed. Mixing measured speed with virtual-model
     * direction can make aerodynamic, rolling and incline terms disagree in
     * sign whenever the physical plant temporarily does not track the model.
     *
     * F_combine is a SIGNED load-force quantity used as:
     *
     *     m * dv/dt = F_drive - F_combine
     *
     * Therefore:
     *
     *     forward model motion  -> F_combine > 0
     *     reverse model motion  -> F_combine < 0
     *
     * A negative F_combine in reverse operation is not an inverted loss:
     * subtracting it produces a positive force that opposes reverse motion.
     */
    const float motor_mech_radps =
        rpm * erpm_to_mech_radps;

    const float model_speed =
        motor->model_v;

    const float model_erpm_for_losses =
        model_speed *
        gearing *
        motor->c_wheel_radius_inv *
        mech_radps_to_erpm;

    float model_direction = 0.0f;

    /*
     * Small deadband prevents resistance-direction chatter around zero.
     * d_erpm_soll is used only as a fallback and is itself a model signal.
     */
    const float direction_deadband_mps = 1.0e-4f;

    if (model_speed > direction_deadband_mps) {
        model_direction = 1.0f;
    } else if (model_speed < -direction_deadband_mps) {
        model_direction = -1.0f;
    } else {
        model_direction = SIGN(motor->d_erpm_soll);
    }

    /*
     * The existing incline semantics are preserved:
     * a positive incline command represents an uphill load in the current
     * model-travel direction. Thus it always resists travel rather than
     * becoming an assisting downhill force when the model reverses.
     */
    const float slope =
        incline_deg_cmd *
        ((float)M_PI / 180.0f) *
        model_direction;

    // Smooth rolling/grade forces near zero virtual speed.
    float speed_ratio = 0.0f;

    if (motor->p_speed_limit_pos_control_activation > 1.0e-6f) {
        speed_ratio =
            fabsf(model_erpm_for_losses) /
            motor->p_speed_limit_pos_control_activation;
    }

    float smooth_factor = speed_ratio;
    if (smooth_factor > 1.0f) smooth_factor = 1.0f;
    if (smooth_factor < 0.0f) smooth_factor = 0.0f;

    /*
     * Signed load terms. Every velocity-dependent term has the same sign as
     * virtual motion, so subtracting F_combine always dissipates energy.
     */
    const float F_air =
        0.5f *
        motor->p_air_ro *
        motor->p_c_air *
        motor->c_area_s *
        model_speed *
        fabsf(model_speed);

    const float F_roll =
        smooth_factor *
        motor->p_c_rr *
        motor->p_weight *
        9.81f *
        cosf(slope) *
        model_direction;

    const float F_incline =
        smooth_factor *
        motor->p_weight *
        9.81f *
        sinf(slope);

    const float F_bearings =
        smooth_factor *
        (motor->p_c_bw * motor->p_k_v_bw) *
        model_speed;

    const float F_combine =
        F_air +
        F_roll +
        F_incline +
        F_bearings;

    motor->d_f_combine = F_combine;

    /*
     * Electromagnetic torque estimate.
     *
     * Match the March implementation by using the instantaneous dq currents,
     * not iq_filter. Keep the current bicycle-coordinate sign convention at
     * the boundary so the rest of the controller/observer stays consistent.
     */
    const float Te_native =
        motor->m_motor_state.iq * motor->p_kT +
        (motor->m_motor_state.iq * motor->m_motor_state.id) *
        (motor->p_ld - motor->p_lq);

    motor->te_calculated = bike_from_vesc(Te_native);

    // theta filter
    motor->unwrapped_theta_filtered_prev = motor->unwrapped_theta_filtered;
    UTILS_LP_FAST(motor->unwrapped_theta_filtered, motor->unwrapped_theta, 0.15f);

    const float omega_for_leso = motor_mech_radps;

    /*
     * Third-order LESO inputs are all expressed in bicycle coordinates:
     *
     *     theta = unwrapped mechanical position [rad]
     *     omega = mechanical speed [rad/s]
     *     Te    = electromagnetic motor-shaft torque [Nm]
     */
    leso3_step(
        motor,
        dt,
        motor->te_calculated,
        motor->unwrapped_theta,
        omega_for_leso
    );
/*
    // ================= FREEWHEEL =================
    const float FW_SLIP_REENG      = 20.0f;
    const float FW_T_REENG         = 2.0f;
    const float FW_T_DISENG        = 0.0f;
    const float FW_T_DISENG_FORCED = -20.0f;

    const float erpm_ratio_freewheel_disengage = 0.90f;
    const float erpm_ratio_forced_disengage    = 0.70f;

    float wheel_erpm = motor->d_erpm_soll;
    float speed_dir = SIGN(wheel_erpm);

    float motor_erpm = motor->m_pll_speed * radps_to_rpm;
    float slip_rpm = (motor_erpm - wheel_erpm) * inv_pole_pairs;

    float Tproj = speed_dir * motor->tp_observed;

    if (!(motor->forced_freewheel || motor->freewheel_active) &&
        (fabsf(motor_erpm) <= fabsf(wheel_erpm) * erpm_ratio_forced_disengage) &&
        (motor->ctrl_sm_state == CTRL_SM_ENABLE) &&
        (Tproj < FW_T_DISENG_FORCED)) {
        motor->forced_freewheel = true;
    }

    if (motor->freewheel_enabled || motor->forced_freewheel) {
        if (!motor->freewheel_active) {
            if ((Tproj <= FW_T_DISENG) ||
                (fabsf(motor_erpm) <= fabsf(wheel_erpm) * erpm_ratio_freewheel_disengage)) {
                motor->freewheel_active = true;
            }
        } else {
            if ((Tproj > FW_T_REENG) &&
                (fabsf(slip_rpm) < FW_SLIP_REENG)) {
                motor->freewheel_active = false;
                motor->forced_freewheel = false;
                motor->model_pos_i_term = 0.0f;
                motor->model_pos_set_model = motor->unwrapped_theta_filtered;
            }
        }
    } else {
        motor->freewheel_active = false;
        motor->forced_freewheel = false;
    }
*/
    /*
     * Freewheel state-machine code is disabled above. Force the associated
     * states off so stale values cannot disable the model or PI loop.
     */
    motor->freewheel_active = false;
    motor->forced_freewheel = false;
    const bool fw = false;

    // ================= Model integration =================
    /*
     * LESO external torque is in bicycle coordinates:
     *
     *     positive tp_observed = rider drives physical bicycle forward
     *
     * Hence positive F_drive increases positive virtual-bicycle speed.
     */
    const float F_drive =
        (motor->tp_observed * motor->c_wheel_radius_inv) *
        gearing;

    float accel_now;
    if (fw) {
        accel_now = (-F_combine) / motor->p_weight;
    } else {
        accel_now = (F_drive - F_combine) / motor->p_weight;
    }

    motor->accel_ist = accel_now;

    motor->model_v += 0.5f * (accel_now + motor->model_accel_prev) * dt;
    motor->model_accel_prev = accel_now;

    const float omega_mech = motor->model_v * gearing * motor->c_wheel_radius_inv;

    motor->model_pos_set_model += omega_mech * dt;
    motor->d_erpm_soll = omega_mech * mech_radps_to_erpm;

    // ====================== CASCADE ======================
    pos_error = motor->model_pos_set_model - motor->unwrapped_theta_filtered;

    p_term_pos = pos_error * pos_kp_eff;

    const float pos_corr_erpm_max = 20000.0f;

    float pos_i = motor->model_pos_i_term;
    float pos_corr_raw = p_term_pos + pos_i;
    float pos_corr_sat = pos_corr_raw;

    utils_truncate_number_abs(&pos_corr_sat, pos_corr_erpm_max);

    const float Tt = 0.10f;
    const float kaw = (Tt > 1e-6f) ? (1.0f / Tt) : 0.0f;

    pos_i += dt * (pos_ki_eff * pos_error + kaw * (pos_corr_sat - pos_corr_raw));

    motor->model_pos_i_term = pos_i;

    pos_corr_raw = p_term_pos + motor->model_pos_i_term;
    pos_corr_sat = pos_corr_raw;
    utils_truncate_number_abs(&pos_corr_sat, pos_corr_erpm_max);
    motor->speed_out_pos_controller = pos_corr_sat;

    float base_speed_ref_erpm = omega_mech * mech_radps_to_erpm;
    float speed_ref_erpm = base_speed_ref_erpm + pos_corr_sat;
    utils_truncate_number(&speed_ref_erpm, conf_now->l_min_erpm, conf_now->l_max_erpm);

    // ==================== Speed loop ====================
    float error = speed_ref_erpm - rpm;
    motor->speed_error = error;

    p_term = error * sp_kp_eff * (1.0f / 20.0f);

    d_term = (error - motor->m_speed_prev_error) * (conf_now->s_pid_kd / dt) * (1.0f / 20.0f);
    UTILS_LP_FAST(motor->m_speed_d_filter, d_term, conf_now->s_pid_kd_filter);
    d_term = motor->m_speed_d_filter;

    motor->m_speed_prev_error = error;

    motor->T_f_combine = F_combine * motor->p_wheel_radius / gearing;

    float speed_ratio2 = erpm_abs_meas / 100.0f;
    float smooth_factor2 = speed_ratio2;
    if (smooth_factor2 > 1.0f) smooth_factor2 = 1.0f;

    float Te_ff = (-motor->Text_ext_hat_f) * smooth_factor2 * motor->p_adrc_scale;

    motor->Te_set = Te_ff;
    motor->iq_set_ff = Te_ff / motor->p_kT;
  

    float iq_ff = motor->iq_set_ff;
    utils_truncate_number_abs(&iq_ff, conf_now->l_current_max);

    float iq_ff_norm = iq_ff * motor->c_iq_norm_inv;

    /*
     * iq_ff is in bicycle coordinates. c_v_q_ff is consumed by the native
     * VESC q-axis current controller, so convert it back.
     */
    motor->c_v_q_ff =
        vesc_from_bike(iq_ff) *
        motor->m_res_est;

    float i_inc = error * sp_ki_eff * dt * (1.0f / 20.0f);
    bool wants_accel = (i_inc > 0.0f);

    if (sp_ki_eff < 1e-9f) {
        motor->m_speed_i_term = 0.0f;
    } else {
        if ((motor->freewheel_active || motor->forced_freewheel) && wants_accel) {
            motor->m_speed_i_term *= 0.98f;
        } else {
            motor->m_speed_i_term += i_inc;
        }
        utils_truncate_number_abs(&motor->m_speed_i_term, 1.0f);
    }

    motor->d_speed_soll = motor->model_v;

    float output = p_term + motor->m_speed_i_term + d_term + iq_ff_norm;
    utils_truncate_number_abs(&output, 1.0f);

    if (!conf_now->s_pid_allow_braking) {
        if (rpm > 20.0f && output < 0.0f) output = 0.0f;
        if (rpm < -20.0f && output > 0.0f) output = 0.0f;
    }

    if (fw) {
        motor->model_pos_i_term = 0.0f;
        motor->m_speed_i_term   = 0.0f;
        motor->m_speed_prev_error = 0.0f;
        motor->m_speed_d_filter = 0.0f;
        output = 0.0f;
        motor->c_v_q_ff = 0.0f;
    }

    if (!ctrl_enabled) {
        motor->m_speed_i_term = 0.0f;
        motor->m_speed_prev_error = 0.0f;
        motor->m_speed_d_filter = 0.0f;
        motor->model_pos_i_term = 0.0f;

        output = 0.0f;
        motor->c_v_q_ff = 0.0f;
        motor->iq_set_ff = 0.0f;
        motor->Te_set = 0.0f;
    }

    /*
     * Diagnostic value kept in bicycle coordinates.
     */
    motor->d_i_res = -iq_ff;

    /*
     * output is a bicycle-coordinate normalized torque/current command.
     * m_iq_set is native VESC q-current, so this is the final inverse
     * coordinate conversion.
     */
    motor->m_iq_set =
        vesc_from_bike(
            output *
            conf_now->lo_current_max *
            conf_now->l_current_max_scale
        );

    motor->status_bits = 0;
    motor->status_bits |= true << STATUS_BIT_SPEED_CONTROL_ACTIVE;
    motor->status_bits |= (motor->forced_freewheel) << STATUS_BIT_FORCED_FREEWHEEL;
    motor->status_bits |= (motor->ctrl_sm_state == CTRL_SM_START) << STATUS_BIT_CTRL_SM_START;
    motor->status_bits |= (motor->ctrl_sm_state == CTRL_SM_INDEX_FOUND) << STATUS_BIT_CTRL_SM_INDEX_FOUND;
    motor->status_bits |= (motor->ctrl_sm_state == CTRL_SM_ENABLE) << STATUS_BIT_CTRL_SM_ENABLE;
    motor->status_bits |= ((uint32_t)1U) << 5;
    motor->status_bits |= ((uint32_t)1U) << 6;
    motor->status_bits |= ((uint32_t)1U) << 7;
    motor->status_bits |= ((uint32_t)1U) << 16;
    motor->status_bits |= ((uint32_t)1U) << 31;

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



float smooth_force(float mag, float v, float v_eps) {
    return mag * (v / v_eps); // smoothly goes negative if v<0
}

// ===== Linear ESO (LESO) for external torque =============================
// States in m: leso_th [rad], leso_om [rad/s], leso_z [rad/s^2]
// Outputs: m->Text_ext_hat_f [Nm] (filtered total torque estimate)
// Inputs: theta_meas = UNWRAPPED angle [rad], Te_meas = kT * iq_meas_f [Nm], dt [s]
/*
 * Active third-order LESO.
 *
 * States:
 *     leso_th = estimated unwrapped mechanical angle [rad]
 *     leso_om = estimated mechanical speed [rad/s]
 *     leso_z  = estimated external acceleration disturbance [rad/s^2]
 *
 * leso_z4 is not an observer state here; it is only an initialization flag.
 */
inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
    float omega_meas
) {
    if (!m || !(dt > 0.0f)) {
        return;
    }

    if (!(m->p_J > 0.0f)) {
        return;
    }

    if (!isfinite(Te_meas) ||
        !isfinite(theta_meas) ||
        !isfinite(omega_meas)) {
        return;
    }

    /*
     * Initialize directly from the measured bicycle-coordinate position and
     * speed. This prevents a large innovation when entering speed mode with
     * an already nonzero unwrapped position.
     */
    if (m->leso_z4 != 1.0f ||
        !isfinite(m->leso_th) ||
        !isfinite(m->leso_om) ||
        !isfinite(m->leso_z) ||
        !isfinite(m->Text_ext_hat_f)) {

        m->leso_th = theta_meas;
        m->leso_om = omega_meas;
        m->leso_z = 0.0f;
        m->leso_z4 = 1.0f;

        m->Tf_hat = m->p_B * omega_meas;
        m->Text_ext_hat = 0.0f;
        m->Text_ext_hat_f = 0.0f;
        m->tp_observed = 0.0f;
        return;
    }

    const float B =
        (isfinite(m->p_B) && m->p_B >= 0.0f) ?
        m->p_B :
        0.0f;
    const float b0 = m->c_b0;
    const float b1 = m->c_b1;
    const float b2 = m->c_b2;
    const float b3 = m->c_b3;
    const float gz = m->c_gz;

    const float h = 0.5f * dt;

    // Previous observer states
    const float thk = m->leso_th;
    const float omk = m->leso_om;
    const float zk  = m->leso_z;

    // Position innovation
    const float ek = theta_meas - thk;

    /*
     * LESO model:
     *
     * theta_dot = omega + b1 * e
     *
     * omega_dot =
     *     b0 * (Te - B * omega)
     *     + z
     *     + b2 * e
     *
     * z_dot =
     *     b3 * e
     *     - gz * z
     */

    // Implicit trapezoidal coefficients
    const float a11 = 1.0f + h * b1;
    const float a12 = -h;

    const float a21 = h * b2;
    const float a22 = 1.0f + h * b0 * B;
    const float a23 = -h;

    const float a31 = h * b3;
    const float a33 = 1.0f + h * gz;

    // Explicit part evaluated at state k
    const float fk_th =
        omk +
        b1 * ek;

    const float fk_om =
        b0 * (Te_meas - B * omk) +
        zk +
        b2 * ek;

    const float fk_z =
        b3 * ek -
        gz * zk;

    // Right-hand sides
    const float rhs1 =
        thk +
        h * fk_th +
        h * b1 * theta_meas;

    const float rhs2 =
        omk +
        h * fk_om +
        h * (b0 * Te_meas + b2 * theta_meas);

    const float rhs3 =
        zk +
        h * fk_z +
        h * b3 * theta_meas;

    // Sparse implicit solve
    const float inv_a11 = 1.0f / a11;
    const float inv_a33 = 1.0f / a33;

    const float th_const =
        rhs1 * inv_a11;

    const float th_om =
        -a12 * inv_a11;

    const float z_const =
        (rhs3 - a31 * th_const) * inv_a33;

    const float z_om =
        -(a31 * th_om) * inv_a33;

    const float const2 =
        a21 * th_const +
        a23 * z_const;

    const float coeff2 =
        a21 * th_om +
        a22 +
        a23 * z_om;

    if (!isfinite(coeff2) || fabsf(coeff2) <= 1e-9f) {
        return;
    }

    const float om1 =
        (rhs2 - const2) / coeff2;

    const float th1 =
        th_const + th_om * om1;

    const float z1 =
        z_const + z_om * om1;

    if (!isfinite(th1) ||
        !isfinite(om1) ||
        !isfinite(z1)) {
        return;
    }

    // Commit pure LESO states: no clamps or slew limits
    m->leso_th = th1;
    m->leso_om = om1;
    m->leso_z  = z1;

    // Diagnostic friction estimate
    m->Tf_hat = B * omega_meas;

    /*
     * Friction is already explicitly included in the observer plant:
     *
     *     Te - B * omega
     *
     * Therefore J*z is the external torque referred to the motor shaft.
     */
    const float Text_motor_hat =
        m->p_J * z1;

    m->Text_ext_hat = Text_motor_hat;

    // Retain the torque-output LPF.
    if (isfinite(m->c_fc_2pi) && m->c_fc_2pi > 0.0f) {
        const float aT =
            expf(-m->c_fc_2pi * dt);

        m->Text_ext_hat_f =
            aT * m->Text_ext_hat_f +
            (1.0f - aT) * Text_motor_hat;
    } else {
        m->Text_ext_hat_f =
            Text_motor_hat;
    }

    // Filtered observer torque used by the bicycle model/controller
    m->tp_observed =
        m->Text_ext_hat_f;
}

#if 0
/*
 * Simplified speed-based external-torque estimator retained for diagnostic
 * comparison. It is not compiled while the third-order LESO is active.
 *
 * Mechanical model:
 *
 *     J * alpha = Te - B * omega + Text
 *
 * therefore:
 *
 *     Text = J * alpha - Te + B * omega
 *
 * The estimator uses only measured mechanical speed and electromagnetic
 * torque. theta_meas is intentionally ignored.
 *
 * Existing fields are reused, so motor_all_state_t does not need changes:
 *
 *     leso_th = initialization flag, 0.0f or 1.0f
 *     leso_om = filtered measured mechanical speed [rad/s]
 *     leso_z  = corrected external disturbance acceleration Text / J
 *     leso_z4 = learned constant torque bias [Nm]
 */
inline void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
    float omega_meas
) {
    (void)theta_meas;

    if (!m) {
        return;
    }

    /*
     * Reject invalid or unexpectedly long intervals. Reinitialize on the
     * next valid call instead of differentiating across a timing gap.
     */
    if (!(dt > 0.0f) ||
        dt > 0.020f ||
        !isfinite(dt) ||
        !isfinite(Te_meas) ||
        !isfinite(omega_meas) ||
        !(m->p_J > 0.0f) ||
        !isfinite(m->p_J)) {

        m->leso_th = 0.0f;
        m->leso_z = 0.0f;

        m->Tf_hat = 0.0f;
        m->Text_ext_hat = 0.0f;
        m->Text_ext_hat_f = 0.0f;
        m->tp_observed = 0.0f;
        return;
    }

    /*
     * Use zero viscous friction if B is invalid. A negative B would model
     * anti-friction and must not be used in the reconstruction.
     */
    const float B =
        (isfinite(m->p_B) && m->p_B >= 0.0f) ?
        m->p_B :
        0.0f;

    /*
     * Speed-filter bandwidth.
     *
     * Reuse the existing LESO parameter p_fo_hz so it can be changed
     * through the existing Python/parameter interface.
     *
     * In this simplified estimator:
     *
     *     p_fo_hz = measured-speed LPF cutoff [Hz]
     *
     * This is no longer an LESO observer bandwidth.
     */
    float f_speed_hz = m->p_fo_hz;

    /*
     * Safe fallback for invalid parameter values.
     */
    if (!isfinite(f_speed_hz) || !(f_speed_hz > 0.0f)) {
        f_speed_hz = 8.0f;
    }

    const float wc_speed =
        2.0f * (float)M_PI * f_speed_hz;

    const float a_speed =
        expf(-wc_speed * dt);

    /*
     * Initialize directly from measured speed and generate no acceleration
     * or feedforward impulse on the first valid iteration.
     */
    if (m->leso_th != 1.0f ||
        !isfinite(m->leso_om) ||
        !isfinite(m->Text_ext_hat_f)) {

        m->leso_th = 1.0f;
        m->leso_om = omega_meas;
        m->leso_z = 0.0f;

        m->Tf_hat = B * omega_meas;
        m->Text_ext_hat = 0.0f;
        m->Text_ext_hat_f = 0.0f;
        m->tp_observed = 0.0f;
        return;
    }

    const float omega_f_old =
        m->leso_om;

    const float omega_f_new =
        a_speed * omega_f_old +
        (1.0f - a_speed) * omega_meas;

    /*
     * Derivative of filtered speed: a first-order low-pass
     * differentiator.
     */
    const float alpha_f =
        (omega_f_new - omega_f_old) / dt;

    const float Tf_hat =
        B * omega_f_new;

    const float Text_motor_raw =
        m->p_J * alpha_f
        - Te_meas
        + Tf_hat;

    if (!isfinite(omega_f_new) ||
        !isfinite(alpha_f) ||
        !isfinite(Text_motor_raw)) {

        m->leso_th = 0.0f;
        m->leso_z = 0.0f;

        m->Tf_hat = 0.0f;
        m->Text_ext_hat = 0.0f;
        m->Text_ext_hat_f = 0.0f;
        m->tp_observed = 0.0f;
        return;
    }

    /*
     * Automatic constant-bias estimation without adding a new state.
     *
     * Existing field reuse:
     *
     *     leso_z4 = estimated zero-torque bias [Nm]
     *
     * Existing parameter reuse:
     *
     *     p_gz_hz = bias-estimator LPF cutoff [Hz]
     *
     * Update the bias only while the controller is waiting in
     * CTRL_SM_INDEX_FOUND and the mechanism is nearly stationary.
     * Once CTRL_SM_ENABLE is entered, the bias is frozen so sustained
     * rider torque cannot be learned away.
     */
    float f_bias_hz = m->p_gz_hz;

    if (!isfinite(f_bias_hz) || !(f_bias_hz > 0.0f)) {
        f_bias_hz = 1.0f;
    }

    const bool bias_calibration_active =
        (m->ctrl_sm_state == CTRL_SM_INDEX_FOUND) &&
        (fabsf(omega_f_new) < 0.5f) &&
        (fabsf(Te_meas) < 0.5f);

    if (bias_calibration_active) {
        const float a_bias =
            expf(
                -2.0f *
                (float)M_PI *
                f_bias_hz *
                dt
            );

        m->leso_z4 =
            a_bias * m->leso_z4 +
            (1.0f - a_bias) * Text_motor_raw;
    }

    if (!isfinite(m->leso_z4)) {
        m->leso_z4 = 0.0f;
    }

    const float Text_motor_corrected =
        Text_motor_raw - m->leso_z4;

    m->leso_om = omega_f_new;
    m->leso_z = Text_motor_corrected / m->p_J;

    m->Tf_hat = Tf_hat;
    m->Text_ext_hat = Text_motor_corrected;

    /*
     * Retain the existing torque-output LPF.
     * c_fc_2pi = 2*pi*p_fc_TLPF.
     */
    if (isfinite(m->c_fc_2pi) && m->c_fc_2pi > 0.0f) {
        const float a_torque =
            expf(-m->c_fc_2pi * dt);

        m->Text_ext_hat_f =
            a_torque * m->Text_ext_hat_f +
            (1.0f - a_torque) * Text_motor_corrected;
    } else {
        m->Text_ext_hat_f =
            Text_motor_corrected;
    }

    /*
     * Motor-shaft external torque used by the bicycle model and
     * feedforward path. There is no estimator-state clamp.
     */
    m->tp_observed =
        m->Text_ext_hat_f;
}
#endif


float fal_gain(float e, float alpha, float delta, float g0) {
    float ae = fabsf(e);
    if (ae <= delta) {
        return g0 * e;  // slope = g0 at origin
    } else {
        float scale = g0 * powf(delta, 1.0f - alpha);
        return copysignf(scale * powf(ae, alpha), e);
    }
}


inline float falf(float e, float alpha, float delta) {
    float ae = fabsf(e);
    if (ae <= delta) {
        // Linear in the boundary layer, scaled for continuity
        return e / powf(delta, 1.0f - alpha);
    } else {
        return copysignf(powf(ae, alpha), e);
    }
}


inline float ramp_rational_ref(float x, float x_ref, float p) {
    // Choose how close to 1 you want at x_ref. 0.999 = "practically saturated".
    const float target = 0.999f;

    if (!(x > 0.0f)) return 0.0f;
    if (!(x_ref > 1e-6f)) return 1.0f;
    if (!(p > 1e-3f)) p = 1.0f;

    // x0 = x_ref / ( (target/(1-target))^(1/p) )
    const float ratio = target / (1.0f - target);
    const float x0    = x_ref / powf(ratio, 1.0f / p);

    float r = powf(x / x0, p);
    float m = r / (1.0f + r);

    // numeric safety
    if (m < 0.0f) m = 0.0f;
    if (m > 1.0f) m = 1.0f;
    return m;
}

inline float map_floor(float m, float floor) {
    // maps m in [0,1] -> scale in [floor, 1]
    if (floor < 0.0f) floor = 0.0f;
    if (floor > 1.0f) floor = 1.0f;
    return floor + (1.0f - floor) * m;
}

inline float ramp_rational_x0(float x, float x0, float p) {
	if (!(x > 0.0f)) return 0.0f;
	if (!(x0 > 1e-6f)) return 1.0f;
	if (!(p > 1e-3f)) p = 1.0f;

	float r = powf(x / x0, p);
	float m = r / (1.0f + r);

	if (m < 0.0f) m = 0.0f;
	if (m > 1.0f) m = 1.0f;
	return m;
}

inline float fal_nleso_erpm(float e_th_rad,
                                  float alpha,
                                  float delta_erpm,
                                  float dt,
                                  float pole_pairs)
{
    // Convert theta residual -> "equivalent ERPM residual"
    // e_erpm ≈ (e_th/dt) [rad/s] * k_erpm
    // k_erpm = (60/2π)*pp
    const float k_erpm = (60.0f / (2.0f * (float)M_PI)) * pole_pairs;

    // protect against dt=0 (already guarded by caller, but keep it safe)
    if (!(dt > 0.0f) || !(k_erpm > 0.0f)) {
        return e_th_rad; // fallback: linear
    }

    const float e_erpm = (e_th_rad / dt) * k_erpm;

    const float ae = fabsf(e_erpm);
    if (ae <= delta_erpm) {
        // unit small-signal gain: ek = e_th
        return e_th_rad;
    } else {
        // ek/e_th gain = (delta_erpm/|e_erpm|)^(1-alpha)
        // implement in ERPM domain, then map back to theta (rad)
        const float scale = powf(delta_erpm, 1.0f - alpha);
        const float e_erpm_nl = copysignf(scale * powf(ae, alpha), e_erpm);

        // back to theta residual: e_th_nl = (e_erpm_nl / k_erpm) * dt
        return (e_erpm_nl / k_erpm) * dt;
    }
}

// linear map |omega| in [0..w1] to rate in [rate0..rate1]
inline float rate_from_abs_omega(float om_abs, float w1,
                                        float rate0, float rate1) {
    if (!(w1 > 1e-6f)) return rate1;
    float t = clampf(om_abs / w1, 0.0f, 1.0f);
    return rate0 + (rate1 - rate0) * t;
}

inline float slew_limit(float x, float x_prev, float rate, float dt) {
    const float dx_max = rate * dt;
    return x_prev + clampf(x - x_prev, -dx_max, dx_max);
}

inline float clampf(float x, float lo, float hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}
// ===== Nonlinear 4th-order ESO (NLESO4) for external torque ==================
// States in m:
//   leso_th [rad], leso_om [rad/s], leso_z [rad/s^2], leso_z4 [rad/s^3]   <-- ADD leso_z4
// Output:
//   m->Text_ext_hat_f [Nm] (filtered total torque estimate)
// Inputs:
//   theta_meas = UNWRAPPED angle [rad], Te_meas [Nm], omega_meas [rad/s], dt [s]
//
// Notes:
// - z4 is INTERNAL ONLY: it shapes z (z3). We still use ONLY z (z3) for torque estimation.
// - A z4 clamp is included but COMMENTED OUT for now (as requested).
void nleso4_step_ext_torque(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
    float omega_meas
){
    if (!m || dt <= 0.0f) return;
    const float J = m->p_J; if (!(J > 0.0f)) return;

    const float fo_hz   = m->p_fo_hz;
    const float gz_hz   = m->p_gz_hz;      // leak on z3 (can cause bias if large and no outer I)
    const float fc_TLPF = m->p_fc_TLPF;

    const float B  = m->p_B;
    const float b0 = 1.0f / J;

    // ---------- Coulomb friction ----------
    const float Tc = m->p_Tc;
    float ws = m->p_Tc_ws;
    if (!(ws > 1e-6f)) ws = 1.0f;

    const float Tc_term = Tc * (omega_meas / ws);
    m->Tf_hat = Tc_term + B * omega_meas;

    const float Te_eff = Te_meas - Tc_term;
    // -------------------------------------

    const float wo = 2.0f * (float)M_PI * fo_hz;

    // 4th-order canonical observer gains ( (s+wo)^4 )
    const float b1 = 4.0f * wo;
    const float b2 = 6.0f * wo * wo;
    const float b3 = 4.0f * wo * wo * wo;
    const float b4 =        wo * wo * wo * wo;

    const float gz = 2.0f * (float)M_PI * gz_hz;

    // z4 leak/damping: keep it fairly strong to prevent noise accumulation in the extra state.
    // Start with 4*wo; if z becomes noisy, increase this multiplier.
    const float g4 = 8.0f * wo;

    const float h = 0.5f * dt;

    // Old states
    const float thk = m->leso_th;
    const float omk = m->leso_om;
    const float zk  = m->leso_z;
    const float z4k = m->leso_z4;   // <-- NEW STATE

    // ===================== NLESO innovation (ERPM delta) =====================
    const float e_th = theta_meas - thk;

    float pole_pairs = 23.0f; // fallback
    if (m->m_conf) {
        pole_pairs = 0.5f * (float)m->m_conf->si_motor_poles;
        if (!(pole_pairs > 0.0f)) pole_pairs = 23.0f;
    }

    const float alpha_nl   = 0.8f;     // tune 0.4..0.8 (lower = more spike suppression)
    const float delta_erpm = 100.0f;   // start around 50..150 ERPM for spike shaping

    const float ek = fal_nleso_erpm(e_th, alpha_nl, delta_erpm, dt, pole_pairs);
    // =========================================================================

    // --- Linear system coefficients (implicit for th/om/z, as before) ---
    const float a11 = 1.0f + h * b1;
    const float a12 = -h;

    const float a21 =  h * b2;
    const float a22 = 1.0f + h * (b0 * B);
    const float a23 = -h;

    const float a31 =  h * b3;
    const float a33 = 1.0f + h * gz;

    const float rhs1 = thk + h * (omk + b1 * ek) + h * (b1 * theta_meas);

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
    if (fabsf(coeff2) > 1e-6f) {
        om1 = (rhs2 - const2) / coeff2;
    }

    float th1 = th_const + th_om * om1;
    float z1  = z_const  + z_om  * om1;

    // ====================== CLAMPS (post-solve) ======================

    // Torque-based z bound
    float Te_max = 80.0f; // Nm fallback
    if (m->m_conf) {
        const float Imax = m->m_conf->lo_current_max * m->m_conf->l_current_max_scale; // A
        Te_max = fabsf(Imax) * fabsf(m->p_kT) * 1.2f + 0.5f; // Nm
    }
    const float z_abs_max_torque = Te_max / J;

    // --- ω plausibility clamp + slew on LESO omega state (om1) ---
    {
        const float om_meas = omega_meas;
        const float omk_loc = omk;

        float om_abs_max = 500.0f; // rad/s fallback
        if (m->m_conf) {
            const float pp = 0.5f * m->m_conf->si_motor_poles;
            const float erpm_max = (float)m->m_conf->l_max_erpm;
            const float om_mech_max = (erpm_max / pp) * (2.0f * (float)M_PI / 60.0f);
            om_abs_max = 1.2f * om_mech_max;
        }

        const float om_floor = 5.0f;   // rad/s
        const float rel_band = 0.20f;  // fraction of |omega|

        float dom_allow = om_floor + rel_band * fabsf(om_meas);
        if (dom_allow > om_abs_max) dom_allow = om_abs_max;

        float om_target = om1;
        const float om_lo = om_meas - dom_allow;
        const float om_hi = om_meas + dom_allow;
        if (om_target < om_lo) om_target = om_lo;
        if (om_target > om_hi) om_target = om_hi;

        // slew: slow near zero, fast at speed
        const float w_s    = 20.0f;   // rad/s
        const float rate0w = 300.0f;  // rad/s^2
        const float rate1w = 6000.0f; // rad/s^2

        float s = fabsf(om_meas) / w_s;
        if (s > 1.0f) s = 1.0f;
        s = s * s * (3.0f - 2.0f * s); // smoothstep

        const float dom_rate = rate0w + (rate1w - rate0w) * s;

        float dom = om_target - omk_loc;
        const float dom_max = dom_rate * dt;
        if (dom >  dom_max) dom =  dom_max;
        if (dom < -dom_max) dom = -dom_max;

        om1 = omk_loc + dom;

        if (om1 >  om_abs_max) om1 =  om_abs_max;
        if (om1 < -om_abs_max) om1 = -om_abs_max;

        th1 = th_const + th_om * om1;
    }

    // ===================== NEW: z4 dynamics (internal only) =====================
    // Update z4 as a fast, damped state driven by the same innovation.
    // This gives z a "2nd-order" behavior without feeding z4 to torque.
    float z4_1 = z4k + dt * (b4 * ek - g4 * z4k);

    // --- z4 clamp (DISABLED FOR NOW) ---
    // Uncomment if z4 becomes noisy. This is a physics-ish jerk bound: allow z to traverse full range in ~t_ramp.
    /*
    const float t_ramp = 0.020f; // 20 ms
    const float z4_abs_max = z_abs_max_torque / fmaxf(t_ramp, 1e-3f); // (rad/s^2)/s = rad/s^3
    utils_truncate_number_abs(&z4_1, z4_abs_max);
    */

    // Let z4 shape z: add a slope component to z (post-solve).
    z1 = z1 + dt * z4_1;

    // z clamp (keep existing torque-based bound)
    utils_truncate_number_abs(&z1, z_abs_max_torque);
    // ===========================================================================

    // --- Slew-limit omega used for reconstruction only ---
    const float w1_recon    = 1000.0f;
    const float rate0_recon = 50.0f;
    const float rate1_recon = 2000.0f;

    const float rate_recon = rate_from_abs_omega(fabsf(omega_meas), w1_recon, rate0_recon, rate1_recon);
    m->leso_omega_in = slew_limit(omega_meas, m->leso_omega_in, rate_recon, dt);
    const float omega_leso = m->leso_omega_in;

    // ==================== Commit ====================
    m->leso_th  = th1;
    m->leso_om  = om1;
    m->leso_z   = z1;
    m->leso_z4  = z4_1;   // <-- NEW STATE COMMIT

    const float Tpedal_ext_hat = J * z1;
    m->tp_observed = Tpedal_ext_hat;

    const float Text_ext_hat = Tpedal_ext_hat + Tc_term + B * omega_leso;

    const float aT = expf(-2.0f * (float)M_PI * fc_TLPF * dt);
    m->Text_ext_hat_f = aT * m->Text_ext_hat_f + (1.0f - aT) * Text_ext_hat;
}

float map_floor_local(float m, float floor) {
    utils_truncate_number(&floor, 0.0f, 1.0f);

    if (m < 0.0f) m = 0.0f;
    if (m > 1.0f) m = 1.0f;

    return floor + (1.0f - floor) * m;
}

float ramp_rational_p2(float x, float x_sat) {
    if (x <= 0.0f) {
        return 0.0f;
    }

    if (x_sat < 1e-6f) {
        return 1.0f;
    }

    // Define x_sat as the point where m = 0.8
    // For m = x^2 / (x^2 + x0^2), that gives x0 = 0.5 * x_sat
    const float x0 = 0.5f * x_sat;

    const float x2  = x * x;
    const float x02 = x0 * x0;

    return x2 / (x2 + x02);
}