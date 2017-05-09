#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

void initialize()
{

/*
//130620 -- this is removed because I have comm_line_det_rate_cat_factor to take me from growing det rate to shrinking det rate
if(!EXPONENTIAL_DET_RATE)
  comm_line_det_rate_cat = comm_line_det_rate_rescue; 
*/

//max force test
spring_started = false;
spring_location = comm_line_init_fil_length;


int ii;
double init_ra,init_ca, init_rd, init_cd;

total_conc = (0.5*LEFT_R_0_TOT + 0.5*LEFT_C_0_TOT)*(comm_line_init_fil_length) + FREE_MONO_0;
#if TWO_KINETO
total_conc += (0.5*RIGHT_R_0_TOT + 0.5*RIGHT_C_0_TOT)*comm_line_init_fil_length;
#endif
prev_free_mono = FREE_MONO_0;
free_mono = prev_free_mono;

zero_tension_binding_frac = 1./(1. + exp(-1.*EPSILON)/KbT);
minimum_tension_for_det_mod = 0.;

shift_counter = 0;

for(ii = 0; ii < LZ+1; ii++)//box is length LZ, need LZ+1 grid pts to include length=0
{
	concentration temp_conc1(0., true);
	concentration temp_conc2(0., false);

	left_res.push_back(temp_conc1);
	left_cat.push_back(temp_conc2);

	#if FOUR_POPULATIONS
	left_res_det.push_back(temp_conc1);
	left_cat_det.push_back(temp_conc2);
	
	#if TWO_KINETO
	right_res.push_back(temp_conc1);
	right_cat.push_back(temp_conc2);
	right_res_det.push_back(temp_conc1);
	right_cat_det.push_back(temp_conc2);
	#endif
	#endif

}

if(comm_line_init_fil_length > LZ)
{
	fprintf(stderr, "Initial filament length is longer than box size.  Exiting\n");
	exit(1);
}


#if START_FILAMENTS_AT_KT
#if !FOUR_POPULATIONS
left_res[comm_line_init_fil_length].set_conc(0.5*LEFT_R_0_TOT);//0.5's for these two lines added on 121001
left_cat[comm_line_init_fil_length].set_conc(0.5*LEFT_C_0_TOT);
left_res[comm_line_init_fil_length].set_prev_conc();
left_cat[comm_line_init_fil_length].set_prev_conc();
#elif !TWO_KINETO//i.e., 4 pops model with just a single kinetochore
left_res[comm_line_init_fil_length].set_conc(0.5*LEFT_R_0_TOT);
left_cat[comm_line_init_fil_length].set_conc(0.5*LEFT_C_0_TOT);
left_res[comm_line_init_fil_length].set_prev_conc();
left_cat[comm_line_init_fil_length].set_prev_conc();
/*left_res_det[comm_line_init_fil_length].set_conc(0.);//0.5*LEFT_R_0_TOT);
left_cat_det[comm_line_init_fil_length].set_conc(0.);//0.5*LEFT_C_0_TOT);
left_res_det[comm_line_init_fil_length].set_prev_conc();
left_cat_det[comm_line_init_fil_length].set_prev_conc();
*/ // don't need to set concentrations to 0. they are initialized to zero!!!!!
#else //two kineto

//all detached filaments start at 0
#if !ASYMMETRIC

symmetric_initial_condition1();

#else

//initial_condition1();
//initial_condition2();
//initial_condition3();
initial_condition4();

#endif

#endif//options
#endif//START_FILAMENTS_AT_KT

#if !START_FILAMENTS_AT_KT
//fprintf(stderr, "comm_line_init_fil_length is %i\n", comm_line_init_fil_length);
for(ii = 0; ii <= comm_line_init_fil_length; ii++)
{
#if !FOUR_POPULATIONS
left_res[ii].set_conc(0.5*LEFT_R_0_TOT/(comm_line_init_fil_length+1.));
left_cat[ii].set_conc(0.5*LEFT_C_0_TOT/(comm_line_init_fil_length+1.));
#else
/*NOTE: following initialization is only ok if kt starts at comm_line_init_fil_length*/
/*if(ii > comm_line_init_fil_length - ATTACHMENT_RANGE)
{
 init_ra = LEFT_R_0_TOT/(comm_line_init_fil_length+1.);
 init_ca = LEFT_C_0_TOT/(comm_line_init_fil_length+1.);
 init_rd = 0.;
 init_cd = 0.;
}
else
{*/

//this starts with 0 filaments attached and all filaments detached
 init_ra = 0.;
 init_ca = 0.;
 init_rd = 0.5*LEFT_R_0_TOT/(comm_line_init_fil_length+1.);
 init_cd = 0.5*LEFT_C_0_TOT/(comm_line_init_fil_length+1.);

//fprintf(stderr, "%g %g\n", init_rd, init_cd);
//}

left_res[ii].set_conc(init_ra);
left_cat[ii].set_conc(init_ca);
left_res_det[ii].set_conc(init_rd);
left_cat_det[ii].set_conc(init_cd);

left_res_det[ii].set_prev_conc();
left_cat_det[ii].set_prev_conc();

#if TWO_KINETO
right_res[ii].set_conc(init_ra);
right_cat[ii].set_conc(init_ca);
right_res_det[ii].set_conc(init_rd);
right_cat_det[ii].set_conc(init_cd);

right_res_det[ii].set_prev_conc();
right_cat_det[ii].set_prev_conc();
right_res[ii].set_prev_conc();
right_cat[ii].set_prev_conc();
#endif

#endif

left_res[ii].set_prev_conc();
left_cat[ii].set_prev_conc();
}
#endif//if !START_FILAMENTS_AT_KT

#if (ASYMMETRIC && FOUR_POPULATIONS && TWO_KINETO && !START_FILAMENTS_AT_KT)
//fprintf(stderr, "asymm\n");
int asymmetry_width = 50;
double unadjusted_conc;
double extra_conc = 0.;
for(ii = 0; ii < asymmetry_width; ii++)
{
	extra_conc += right_res[comm_line_init_fil_length - ii].get_conc() + right_cat[comm_line_init_fil_length - ii].get_conc() + right_res_det[comm_line_init_fil_length - ii].get_conc() + right_cat_det[comm_line_init_fil_length - ii].get_conc();
	right_res[comm_line_init_fil_length - ii].set_conc(0.);
	right_cat[comm_line_init_fil_length - ii].set_conc(0.);
	right_res_det[comm_line_init_fil_length - ii].set_conc(0.);
	right_cat_det[comm_line_init_fil_length - ii].set_conc(0.);
        right_res[comm_line_init_fil_length - ii].set_prev_conc();
        right_cat[comm_line_init_fil_length - ii].set_prev_conc();
        right_res_det[comm_line_init_fil_length - ii].set_prev_conc();
        right_cat_det[comm_line_init_fil_length - ii].set_prev_conc();
}
extra_conc *= 0.5 / ( (double) comm_line_init_fil_length + 1. - (double) asymmetry_width);
for(ii = 0; ii <= comm_line_init_fil_length - asymmetry_width; ii++)
{
        unadjusted_conc = right_res_det[ii].get_conc();
	right_res_det[ii].set_conc(unadjusted_conc + extra_conc);
        unadjusted_conc = right_cat_det[ii].get_conc();
	right_cat_det[ii].set_conc(unadjusted_conc + extra_conc);

right_res_det[ii].set_prev_conc();
right_cat_det[ii].set_prev_conc();
right_res[ii].set_prev_conc();
right_cat[ii].set_prev_conc();
}
for(ii = 0; ii < comm_line_init_fil_length + asymmetry_width; ii++)
{ 
 init_ra = 0.;
 init_ca = 0.;
 init_rd = 0.5*LEFT_R_0_TOT/(comm_line_init_fil_length+asymmetry_width);
 init_cd = 0.5*LEFT_C_0_TOT/(comm_line_init_fil_length+asymmetry_width);

//fprintf(stderr, "%g %g\n", init_rd, init_cd);
//}

left_res[ii].set_conc(init_ra);
left_cat[ii].set_conc(init_ca);
left_res_det[ii].set_conc(init_rd);
left_cat_det[ii].set_conc(init_cd);

left_res_det[ii].set_prev_conc();
left_cat_det[ii].set_prev_conc();
left_res[ii].set_prev_conc();
left_cat[ii].set_prev_conc();
}

#endif


//for(ii = 0; ii < LZ; ii++)
//fprintf(stderr, "%i %g %g %g %g\n", ii, left_res[ii].get_prev_conc(), right_res[ii].get_prev_conc(), left_cat_det[ii].get_prev_conc(), right_cat_det[ii].get_prev_conc());


////////////////////kinetochore/////////////////////////////
#if (!(FOUR_POPULATIONS && TWO_KINETO && ASYMMETRIC))
standard_kt_initial_condition();
#else
//standard_kt_initial_condition();
stretched_kt_initial_condition();
#endif

#if (ASYMMETRIC && FOUR_POPULATIONS && TWO_KINETO && !START_FILAMENTS_AT_KT)
//fprintf(stderr, "asym2\n");
	double right_kt_pos = kt_list[1].get_prev_pos();
	kt_list[1].set_pos(right_kt_pos - asymmetry_width);
	kt_list[1].set_prev_pos();
	kt_list[0].set_pos(kt_list[1].get_pos()+KT_KT_SPRING_LENGTH);
	kt_list[0].set_prev_pos();
#endif

initialize_files();


fprintf(stderr, "size = %i\n", left_res.size());
}


void initialize_files()
{
char command0[96];
char command1[96];
char command2[96];
sprintf(command0, "mkdir /var/tmp/new_and_improved_scratch/");
system(command0);
sprintf(command1, "mkdir /var/tmp/new_and_improved_scratch/master_eq");
system(command1);
sprintf(command2, "mkdir /var/tmp/new_and_improved_scratch/master_eq/output");
system(command2);

sprintf(mtlengthfilename, "/var/tmp/new_and_improved_scratch/master_eq/output/mtlengths%6.6i", TRIALNUMBER);
mtlengthfile = fopen(mtlengthfilename, "w");

#if TWO_KINETO
sprintf(rightmtlengthfilename, "/var/tmp/new_and_improved_scratch/master_eq/output/mtlengths_r%6.6i", TRIALNUMBER);
rightmtlengthfile = fopen(rightmtlengthfilename, "w"); 
#endif 

sprintf(ktposfilename, "/var/tmp/new_and_improved_scratch/master_eq/output/ktpos%6.6i", TRIALNUMBER);
ktposfile = fopen(ktposfilename, "w");

sprintf(tensionname, "/var/tmp/new_and_improved_scratch/master_eq/output/tension%6.6i", TRIALNUMBER);
tensionfile = fopen(tensionname, "w");

sprintf(distribname, "/var/tmp/new_and_improved_scratch/master_eq/output/distrib%6.6i", TRIALNUMBER);
distribfile = fopen(distribname, "w");


//write_restart();
//exit(1);

}


void initialize_command_line_variables()
{
	comm_line_trial = 999999;
	comm_line_restart = false;
	comm_line_time = 5.e4;
	comm_line_numsteps = (long) comm_line_time / dt +1;
	comm_line_free_mono = 1.;//1000.;

//commlinelz must be defined before LZ is used
        comm_line_lz = 2000;

	#if TWO_KINETO
	comm_line_init_fil_length = (LZ - (int)KT_KT_SPRING_LENGTH)/2;
	#else
	comm_line_init_fil_length = LZ/2;
	#endif

	comm_line_depolym = 2.0;
	comm_line_polym = 0.07;//1.e-3;
	comm_line_cat = 0.0002;
	comm_line_rescue = 0.003;
	comm_line_att_rate_rescue = 0.02;
	comm_line_det_rate_rescue = 0.02;
	comm_line_det_rate_cat_factor = 1.0;
	comm_line_p_exp = 0.05;
	comm_line_comp_p_exp = comm_line_p_exp;
	comm_line_d_exp = -0.1;
	comm_line_comp_d_exp = comm_line_d_exp;
	comm_line_r_exp = 0.2;
	comm_line_comp_r_exp = comm_line_r_exp;
	comm_line_c_exp = -0.3;
	comm_line_comp_c_exp = comm_line_c_exp;
	comm_line_a_exp = 0.;//-1.0;//assumption
	comm_line_det_exp = 0.;//0.1;//0.1 for attached filaments, -0.1 for detached filaments
	comm_line_det_exp_cat = comm_line_det_exp;//-0.1;
	comm_line_numskip = 20000;
	comm_line_kt_mt_spring = 10.*KbT/MONO_DIAM/MONO_DIAM;
        comm_line_kt_kt_spring = 50.*KbT/MONO_DIAM/MONO_DIAM;
        comm_line_epsilon = 20.*KbT;
	comm_line_left_r0 = 25.;
	comm_line_attachment_range = 2.0*MONO_DIAM;

	comm_line_load = 0.;

#if INSTANT_FORCE
	comm_line_tdelay = 0.;
#else
	comm_line_tdelay = 1000.;
#endif


	comm_line_dash_pot = false;
	comm_line_relative_velocity_drag = 1.0;
	comm_line_sigmoidal_c = false;
	comm_line_sigmoidal_p = false;
	comm_line_sig_cat_coeff = 1.0;
	comm_line_sig_polym_coeff = 1.;
	comm_line_sig_depolym_coeff = 1.;
	comm_line_sig_res_coeff = 1.;

	comm_line_max_force_test = false;
	comm_line_max_force_test_spring = 5.0*KbT/MONO_DIAM/MONO_DIAM;

}





