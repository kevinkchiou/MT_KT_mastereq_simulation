#include "params.h"
#include "protos1.h"
#include <iostream>
#include <math.h>
#include <vector>


using namespace std;

extern double free_mono;
extern double prev_free_mono;
extern double total_conc;

extern double zero_tension_binding_frac;
extern double minimum_tension_for_det_mod;

extern unsigned long long sim_step;

extern FILE *ktposfile;
extern FILE *mtlengthfile;
extern FILE *rightmtlengthfile;
extern FILE *tensionfile;
extern FILE *distribfile;
extern char ktposfilename[128];
extern char mtlengthfilename[128];
extern char rightmtlengthfilename[128];
extern char tensionname[128];
extern char distribname[128];


extern int shift_counter;//this will count the number of shifts (positive or negative) to the kinetochore box
		//note positive shifts wil make shift_counter decrease by 1 because shift_counter is used to determine the true position during output

        extern int comm_line_trial;
	extern double comm_line_time;
        extern long comm_line_numsteps;
        extern double comm_line_free_mono;
	extern int comm_line_init_fil_length;
        extern double comm_line_depolym;
        extern double comm_line_polym;
        extern double comm_line_cat;
        extern double comm_line_rescue;
	extern double comm_line_att_rate_rescue;
	extern double comm_line_det_rate_rescue;
	extern double comm_line_det_rate_cat_factor;
        extern double comm_line_p_exp, comm_line_comp_p_exp;
        extern double comm_line_d_exp, comm_line_comp_d_exp;
        extern double comm_line_r_exp, comm_line_comp_r_exp;
        extern double comm_line_c_exp, comm_line_comp_c_exp;
	extern double comm_line_a_exp;
	extern double comm_line_det_exp;
	extern int comm_line_lz;
	extern int comm_line_numskip;
	extern double comm_line_load;
	extern int comm_line_restart;
	extern double comm_line_kt_mt_spring;
	extern double comm_line_kt_kt_spring;
	extern double comm_line_epsilon;
	extern double comm_line_left_r0;
	extern double comm_line_det_exp_cat;
	extern double comm_line_attachment_range;
	extern double comm_line_tdelay;
	extern bool comm_line_dash_pot;
	extern double comm_line_relative_velocity_drag;
	extern bool comm_line_sigmoidal_c, comm_line_sigmoidal_p;
	extern double comm_line_sig_cat_coeff, comm_line_sig_polym_coeff, comm_line_sig_depolym_coeff, comm_line_sig_res_coeff;
	extern bool comm_line_max_force_test;
	extern double comm_line_max_force_test_spring;

//max force test
extern bool spring_started;
extern double spring_location;
