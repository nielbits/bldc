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

#ifndef FOC_MATH_H_
#define FOC_MATH_H_

#include "datatypes.h"

// Types
typedef struct {
	float va;
	float vb;
	float vc;
	float v_mag_filter;
	float mod_alpha_filter;
	float mod_beta_filter;
	float mod_alpha_measured;
	float mod_beta_measured;
	float mod_alpha_raw;
	float mod_beta_raw;
	float id_target;
	float iq_target;
	float max_duty;
	float duty_now;
	float phase;
	float phase_cos;
	float phase_sin;
	float i_alpha;
	float i_beta;
	float i_abs;
	float i_abs_filter;
	float i_bus;
	float v_bus;
	float v_alpha;
	float v_beta;
	float mod_d;
	float mod_q;
	float mod_q_filter;
	float id;
	float iq;
	float id_filter;
	float iq_filter;
	float vd;
	float vq;
	float vd_int;
	float vq_int;
	uint32_t svm_sector;
	bool is_using_phase_filters;
} motor_state_t;

typedef struct {
	int sample_num;
	float avg_current_tot;
	float avg_voltage_tot;
} mc_sample_t;

typedef struct {
	void(*fft_bin0_func)(float*, float*, float*);
	void(*fft_bin1_func)(float*, float*, float*);
	void(*fft_bin2_func)(float*, float*, float*);

	int samples;
	int table_fact;
	float buffer[32];
	float buffer_current[32];
	bool ready;
	int ind;
	bool is_samp_n;
	float sign_last_sample;
	float cos_last, sin_last;
	float prev_sample;
	float angle;
	float double_integrator;
	int est_done_cnt;
	float observer_zero_time;
	int flip_cnt;
} hfi_state_t;

typedef struct {
	float x1;
	float x2;
	float lambda_est;
	float i_alpha_last;
	float i_beta_last;
} observer_state;

#define MC_AUDIO_CHANNELS	4

typedef enum {
	MC_AUDIO_OFF = 0,
	MC_AUDIO_TABLE,
	MC_AUDIO_SAMPLED,
} mc_audio_mode;

typedef struct {
	mc_audio_mode mode;

	const float *table[MC_AUDIO_CHANNELS];
	int table_len[MC_AUDIO_CHANNELS];
	float table_voltage[MC_AUDIO_CHANNELS];
	float table_freq[MC_AUDIO_CHANNELS];
	float table_pos[MC_AUDIO_CHANNELS];

	// Double-buffered sampled audio
	const int8_t *sample_table[2];
	int sample_table_len[2];
	bool sample_table_filled[2];
	int sample_table_now;
	float sample_freq;
	float sample_pos;
	float sample_voltage;
} mc_audio_state;

typedef enum {
	FOC_PWM_DISABLED = 0,
	FOC_PWM_ENABLED,
	FOC_PWM_FULL_BRAKE
} foc_pwm_mode;
// Put this near your other local defines / enums, e.g. in the relevant motor control source file
typedef enum {
    STATUS_BIT_SPEED_CONTROL_ACTIVE   = 0,
    STATUS_BIT_FORCED_FREEWHEEL       = 1,
    STATUS_BIT_CTRL_SM_START          = 2,
    STATUS_BIT_CTRL_SM_INDEX_FOUND    = 3,
    STATUS_BIT_CTRL_SM_ENABLE         = 4
} status_bits_t;

typedef enum {
	CTRL_SM_START = 0,
	CTRL_SM_INDEX_FOUND = 1,
	CTRL_SM_ENABLE = 2
} ctrl_sm_state_t;

typedef struct {
	mc_configuration *m_conf;
	mc_state m_state;
	mc_control_mode m_control_mode;
	motor_state_t m_motor_state;
	float m_curr_unbalance;
	float m_currents_adc[3];
	bool m_phase_override;
	float m_phase_now_override;
	float m_duty_cycle_set;
	float m_id_set;
	float m_iq_set;
	float m_i_fw_set;
	float m_current_off_delay;
	float m_openloop_speed;
	float m_openloop_phase;
	foc_pwm_mode m_pwm_mode;
	float m_pos_pid_set;
	float m_speed_pid_set_rpm;
	float m_speed_command_rpm;
	float m_phase_now_observer;
	float m_phase_now_observer_override;
	float m_observer_x1_override;
	float m_observer_x2_override;
	bool m_phase_observer_override;
	float m_phase_now_encoder;
	float m_phase_now_encoder_no_index;
	observer_state m_observer_state;
	float m_pll_phase;
	float m_pll_speed;
	float m_speed_est_fast;
	float m_speed_est_fast_corrected; // Same as m_speed_est_fast, but always based on the corrected position
	float m_speed_est_faster;
	mc_sample_t m_samples;
	int m_tachometer;
	int m_tachometer_abs;
	float m_pos_pid_now;
	float m_gamma_now;
	bool m_using_encoder;
	int m_duty1_next, m_duty2_next, m_duty3_next;
	bool m_duty_next_set;
	float m_i_alpha_sample_next;
	float m_i_beta_sample_next;
	float m_i_alpha_sample_with_offset;
	float m_i_beta_sample_with_offset;
	float m_i_alpha_beta_has_offset;
	hfi_state_t m_hfi;
	int m_hfi_plot_en;
	float m_hfi_plot_sample;

	// Audio Modulation
	mc_audio_state m_audio;

	// For braking
	float m_br_speed_before;
	float m_br_vq_before;
	int m_br_no_duty_samples;

	float m_duty_abs_filtered;
	float m_duty_filtered;
	bool m_was_control_duty;
	float m_duty_i_term;
	bool duty_was_pi;
	float duty_pi_duty_last;
	float m_openloop_angle;
	float m_x1_prev;
	float m_x2_prev;
	float m_phase_before_speed_est;
	float m_phase_before_speed_est_corrected;
	int m_tacho_step_last;
	float m_pid_div_angle_last;
	float m_pid_div_angle_accumulator;
	float m_min_rpm_hyst_timer;
	float m_min_rpm_timer;
	bool m_cc_was_hfi;
	float m_pos_i_term;
	float m_pos_prev_error;
	float m_pos_dt_int;
	float m_pos_prev_proc;
	float m_pos_dt_int_proc;
	float m_pos_d_filter;
	float m_pos_d_filter_proc;
	float m_speed_i_term;
	float m_speed_prev_error;
	float m_speed_d_filter;
	int m_ang_hall_int_prev;
	bool m_using_hall;
	float m_ang_hall;
	float m_ang_hall_rate_limited;
	float m_hall_dt_diff_last;
	float m_hall_dt_diff_now;
	bool m_motor_released;

	// Resistance observer
	float m_res_est;
	float m_r_est_state;

	// Temperature-compensated parameters
	float m_res_temp_comp;
	float m_current_ki_temp_comp;

	// Pre-calculated values
	float p_lq;
	float p_ld;
	float p_inv_ld_lq; // (1.0/lq - 1.0/ld)
	float p_v2_v3_inv_avg_half; // (0.5/ld + 0.5/lq)


	//values_for_debugging
	float d_speed;
	float d_f_air;
	float d_f_combine;
	float d_f_bearings;
	float d_f_roll;	
	float d_erpm_soll; //desired rpm
	float d_i_res;


	//filtered speed error and resp. parameters

	float gear_ratio_bike;

	//feedforward control

	float c_v_q_ff;

	//trapezoidal integration and tp observer int


	float accel_ist;

	int_fast64_t tp_observed_fp;  

	//soll speed (model speed)
	float d_speed_soll;
	float d_f_motor;
	//parameters fixed
	float p_air_ro; //air density
	float p_c_rr; //rolling friction
	float p_weight; 
	float p_As; //section area
	float p_c_air; //drag coefficient
	float p_c_bw;
	float p_c_wl;//air resistance coefficient
	float p_wheel_radius; //bike wheel radius;
	float p_mech_gearing;//mechanical gearing from motor to crank = 240/90
	float p_r_bearings;
	float p_k_v_bw;
	float p_kT;
	float p_J;
	float p_B;
	float p_k_area;
	float p_height;
	float p_speed_limit_pos_control_activation;
	float p_incline_deg;
	float p_incline_filtered;
	float p_gear_ratio_filtered;


	//adrc tunables
	float p_fo_hz;      // = 40.0f;          // observer bandwidth (try 10–18 Hz)//100Hz //8Hz for small motor
    float p_gz_hz;      // = 0.00f;           // tiny leak on z to suppress random-walk hiss (0–0.7 Hz)
    float p_fc_TLPF; 	  // = 200.0f;   
	float p_adrc_scale; // =1.0f;
	float p_Tc;           // <-- add this parameter to motor_all_state_t
    float p_Tc_ws; 	// <-- add this parameter too (rad/s)
	//control tunables
	float p_kp_pos;
	float p_ki_pos;
	float p_kd_pos;
	float leso_z4; // [rad/s^3]
	float T_f_combine;//torque caused by_ combined forces (air, rolling, incline, bearings) at speed limit activation point. Used for feedforward compensation in position control.

	bool pumptrack_enabled;
	float pumptrack_time;
	float pumptrack_period_min;
	float incline_result;
	
	// --- freewheel state ---
	bool forced_freewheel;
	bool freewheel_enabled;
	float leso_omega_in;

	//Tp observer calculations
	float tp_observed;
	float te_calculated;
	
	float  model_v;
	float model_accel_prev;


	int32_t status_bits;

	float Tf_hat;


	ctrl_sm_state_t ctrl_sm_state;

	uint32_t ctrl_sm_still_cycles;

	bool freewheel_active;
	float fw_timer_s;


	//position setpoint for cascade control.
	float model_pos_set_model;
	float model_pos_i_term;
	float model_pos_prev_error;
	float model_pos_prev_proc;
	float model_pos_d_filter;
	float model_pos_dt_int;
	float last_tp;
	float last_rpm;
	float rpm_inc_filter_th;
	
	float unwrapped_theta;
	float unwrapped_theta_filtered;
	float unwrapped_theta_filtered_prev;
	float last_angle_rad;
	float last_delta_rad;
	//erpm simulation for frequency response test.

	float simulated_erpm;
	float erpm_time;


	float leso_th;    // theta_hat [rad]
	float leso_om;    // omega_hat [rad/s]
	float leso_z;     // disturbance acceleration z_hat [rad/s^2]

	float Tdist_total_hat;     // Nm  (total matched disturbance torque estimate)
	float Tdist_total_hat_f;   // Nm  (LPF version for cancellation)

	// Te_feed_forward
	float Te_set;         // (optional) motor torque setpoint [Nm] for logging
	float iq_set_ff;      // (optional) current feedforward [A]


	// LESO outputs
	float Text_ext_hat;    // [Nm] instantaneous external torque estimate
	float Text_ext_hat_f;  // [Nm] LPF'd external torque estimate

	//kalman removed previously



	// --- Angle Unwrapping using encoder_read_deg()
	float kalman_last_angle_rad;          // Last raw mechanical angle in degrees
	float kalman_last_delta_rad;          // Last raw mechanical angle in degrees

	float speed_out_pos_controller;
	float speed_error;


	int param_index;
	float param_value;


	// ================= Cached constants / gains =================
	float c_pole_pairs;
	float c_inv_pole_pairs;

	float c_radps_to_rpm;
	float c_rpm_to_radps;
	float c_mech_radps_to_erpm;
	float c_erpm_to_mech_radps;

	float c_area_s;              // p_k_area * p_height^2
	float c_wheel_radius_inv;    // 1 / p_wheel_radius

	float c_iq_norm_inv;         // 1 / (lo_current_max * l_current_max_scale)
	float c_z_abs_max;           // cached LESO z clamp
	float c_om_abs_max;          // cached omega plausibility clamp

	// Speed scheduling cache
	float c_erpm_act;
	float c_erpm_sat;
	float c_inv_erpm_sat;

	float c_pos_dead;
	float c_pos_floor;
	float c_spd_floor;
	float c_ref_pos;
	float c_inv_ref_pos;

	// LESO cache
	float c_b0;
	float c_b1;
	float c_b2;
	float c_b3;
	float c_gz;
	float c_fc_2pi;             // 2*pi*fc_TLPF

} motor_all_state_t;


// Functions
void foc_observer_update(float v_alpha, float v_beta, float i_alpha, float i_beta,
		float dt, observer_state *state, float *phase, motor_all_state_t *motor);
void foc_pll_run(float phase, float dt, float *phase_var,
		float *speed_var, mc_configuration *conf);
void foc_svm(float alpha, float beta, uint32_t PWMFullDutyCycle,
		uint32_t* tAout, uint32_t* tBout, uint32_t* tCout, uint32_t *svm_sector);
void foc_run_pid_control_pos(bool index_found, float dt, motor_all_state_t *motor);
void foc_run_pid_control_speed(bool index_found, float dt, motor_all_state_t *motor);
float foc_correct_encoder(float obs_angle, float enc_angle, float speed, float sl_erpm, motor_all_state_t *motor);
float foc_correct_hall(float angle, float dt, motor_all_state_t *motor, int hall_val);
void foc_run_fw(motor_all_state_t *motor, float dt);
void foc_hfi_adjust_angle(float ang_err, motor_all_state_t *motor, float dt);
void foc_precalc_values(motor_all_state_t *motor);
float band_pass_filter(float input,float centerFreq,   float bandwidth, float sampleRate, motor_all_state_t *motor);
//static inline float friction_T(float omega, float Tc, float b);
//static inline float sgn_db(float x, float dead);
//static inline float dTf_domega(float omega, float b);
//static inline void stribeck_tf_and_dtf(    const motor_all_state_t *m,    float omega,    float *Tf_out,    float *dTf_out);
//inline float Tf_smooth(float omega, float Tc, float B, float omega_s);
void leso3_step(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
	float omega
);
float clampf(float x, float lo, float hi);
float slew_limit(float x, float x_prev, float rate, float dt);
float rate_from_abs_omega(float om_abs, float w1,float rate0, float rate1);
float map_floor(float m, float floor);
float ramp_rational_x0(float x, float x0, float p);

/*
inline float falf(float e, float alpha, float delta);
inline float ramp_rational_ref(float x, float x_ref, float p);



float fal_gain(float e, float alpha, float delta, float g0);
float smooth_force(float mag, float v, float v_eps);
inline float rate_from_abs_omega(float om_abs, float w1,float rate0, float rate1);
inline float fal_nleso_erpm(float e_th_rad,
                                  float alpha,
                                  float delta_erpm,
                                  float dt,
                                  float pole_pairs);
*/

								  /*
void nleso4_step_ext_torque(
    motor_all_state_t *m,
    float dt,
    float Te_meas,
    float theta_meas,
    float omega_meas
);

*/
#endif /* FOC_MATH_H_ */