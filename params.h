#define TRIALNUMBER comm_line_trial
#define RESTART ((bool) comm_line_restart)
#define PERTURB_RESTART false
#define INSTANT_FORCE false
#define FOUR_POPULATIONS true 


//num objects, initial concs
#define TWO_KINETO false
#if (TWO_KINETO && FOUR_POPULATIONS)
#define NUMBER_OF_KTS 2
#define ASYMMETRIC true
#else
#define NUMBER_OF_KTS 1
#endif


#define DASH_POT comm_line_dash_pot


//concentrations in micromole
#define LEFT_R_0_TOT comm_line_left_r0//10.//~5nM? --- 30 filaments, but tip conc is very low!
#define LEFT_C_0_TOT LEFT_R_0_TOT
#define RIGHT_R_0_TOT LEFT_R_0_TOT
#define RIGHT_C_0_TOT LEFT_R_0_TOT
#define START_FILAMENTS_AT_KT true
#define FREE_MONO_0 comm_line_free_mono//~5micromole? in num monos per box, box assumed to be 1micron x 1 micron x 10micron
//1000
#define CONST_FREE_MONO true 

//kinetic rates, note that rates are multiplied by dt in monomer object
//rates in nanomole / tau = nanomole per .5millisec
#define BASE_POLYM_RATE comm_line_polym//1.e-3//4.e-4// num monos per 0.5 ms//micromole/sec = 3.e-5//need free*polym = 1.16e-4//growth polym rate // 7 nm/s
#define BASE_DEPOLYM_RATE comm_line_depolym//1.e-2// 200 nm/s
//1.e-3
#define BASE_CAT_RATE comm_line_cat//1.e-6//1.e-6//catastrophe rate // 8/hr
//1.e-5
#define BASE_RES_RATE comm_line_rescue//1.e-5//1.e-5 // 90/hr
#define BASE_ATT_RATE_RESCUE comm_line_att_rate_rescue//0.02s^-1
#define BASE_DET_RATE_RESCUE comm_line_det_rate_rescue//0.015s^-1
#define BASE_ATT_RATE_CAT BASE_ATT_RATE_RESCUE
#define BASE_DET_RATE_CAT (comm_line_det_rate_cat_factor*BASE_DET_RATE_RESCUE)
#define ATTACHMENT_RANGE comm_line_attachment_range


// these constants are measured in inverse force units, so they only depend on KbT / MONO_DIAM
#define POLYM_EXP_FACTOR comm_line_p_exp//1.e-2///3.e-2///nM/.5ms/pN //rate up to 60nm/s with force of 16pN... .13/pN const = ln(y2/y1) / (x2-x1)
#define COMPRESSION_POLYM_EXP_FACTOR comm_line_comp_p_exp
#define DEPOLYM_EXP_FACTOR comm_line_d_exp//(-5.e-2)//-6.e-2////rate down to 20 nm/s with force of 7.5pN... -.31/pN
#define COMPRESSION_DEPOLYM_EXP_FACTOR comm_line_comp_d_exp
#define RES_EXP_FACTOR comm_line_r_exp//1.e-1//8.e-2 //rate up to 600/hr w/ force of 10pN... .43/pN
#define COMPRESSION_RES_EXP_FACTOR comm_line_comp_r_exp
#define CAT_EXP_FACTOR comm_line_c_exp//(-1.e-1)// //rate down to 0.1/hr w/ force of 9.5pN... -.71/pN
#define COMPRESSION_CAT_EXP_FACTOR comm_line_comp_c_exp
#define ATT_RATE_EXP_FACTOR_RESCUE comm_line_a_exp
#define DET_RATE_EXP_FACTOR_RESCUE comm_line_det_exp// 0.7/hr @ 0pN up to 20/hr @ 13pN... .26/pN
#define ATT_RATE_EXP_FACTOR_CAT ATT_RATE_EXP_FACTOR_RESCUE
#define DET_RATE_EXP_FACTOR_CAT comm_line_det_exp_cat//100/hr @ 0pN down to 20/hr @ 7pN... -0.23/pN

#define SIGMOIDAL_POLYM_RATE comm_line_sigmoidal_p
#define SIGMOIDAL_DEPOLYM_RATE SIGMOIDAL_POLYM_RATE
#define SIGMOIDAL_CAT_RATE comm_line_sigmoidal_c
#define SIGMOIDAL_RES_RATE SIGMOIDAL_CAT_RATE
#define SIGMOID_POLYM_COEFF comm_line_sig_polym_coeff
#define SIGMOID_DEPOLYM_COEFF comm_line_sig_depolym_coeff
#define SIGMOID_CAT_COEFF comm_line_sig_cat_coeff
#define SIGMOID_RES_COEFF comm_line_sig_res_coeff



//These determine the functional form of the on/off rates as a func of tension/ mt dist from kt
#define EXPONENTIAL_ATT_RATE false
#define LINEAR_ATT_RATE false

#define KRAMERS_DET false
#define EXPONENTIAL_DET_RATE false

#define CONST_ATTDET_RATES true//applies to both 2 and four pop model
#define INF_DETACH false//if attachment only allowed within attachment range

//for two pop. model
#define EXPONENTIAL_ATTACHMENT_PROFILE false 
#define FRAC_ATT_EXP_FACTOR (-1./(3.*KT_MT_SPRING))

//properties of objects
#define MONO_DIAM 1.
#define EPSILON comm_line_epsilon//(10.*KbT) // kt_mt binding energy
#define KT_MT_SPRING comm_line_kt_mt_spring//(30.*KbT/MONO_DIAM/MONO_DIAM)//30.
#define KT_KT_SPRING comm_line_kt_kt_spring//0.//(50.*KbT)
#define KT_KT_SPRING_LENGTH (10.*MONO_DIAM)
#define KT_DIFF_COEFF 1.//e.g.,  D = 1e-3 um^2/s
#define KT_DRAG (KbT / KT_DIFF_COEFF) 

#define RELATIVE_VELOCITY_DRAG (comm_line_relative_velocity_drag*KT_DRAG)//zeta prime
#define SUM_DRAGS (KT_DRAG+RELATIVE_VELOCITY_DRAG)//zeta+zeta prime
#define EFFECTIVE_DRAG (KT_DRAG*KT_DRAG + 2.*KT_DRAG*RELATIVE_VELOCITY_DRAG)// zeta^2 + 2 zeta zeta prime
#define DRAG_OVER_EFFECTIVE (KT_DRAG / EFFECTIVE_DRAG)//zeta / effective
#define RV_DRAG_OVER_EFFECTIVE (RELATIVE_VELOCITY_DRAG / EFFECTIVE_DRAG)// zeta prime / effective
#define SUM_DRAG_OVER_EFFECTIVE (SUM_DRAGS / EFFECTIVE_DRAG)




//properties of simulation and box
#define dt 1e-3//5e-4//1e-4//1e-3//120709: used to be 2e-4
#define TOT_TIME comm_line_time
#define NUMSTEPS comm_line_numsteps
#define NUMSKIP comm_line_numskip

#define LZ comm_line_lz
#define INFINITE_BOX true

//other
#define PULLING_LOAD comm_line_load
#define MAX_FORCE_TEST comm_line_max_force_test
#define MAX_FORCE_TEST_SPRING comm_line_max_force_test_spring
//note 1 unit of force, KbT/MONO_DIAM = .138pN
#define TIME_DELAY comm_line_tdelay//2000.

//testing
//#define TENSION_INDEPENDENT true





#define POSITIVE_SHIFT true
#define NEGATIVE_SHIFT false
#if !TWO_KINETO
#define SHIFT_LENGTH (0.25*LZ)
#define EDGE_PROXIMITY_LENGTH (0.25*LZ)
#else
#define SHIFT_LENGTH (0.15*LZ)
#define EDGE_PROXIMITY_LENGTH (0.15*LZ)
#endif




/////////// things that never change /////////////////////////////////
#define DIMENSION 3
#define KbT 1.0
#define inv_dt (1./dt)

#define unif_rand() (ranf0())

//flags
#define PREVIOUS true
#define CURRENT !PREVIOUS

//numbers
#define PI 3.141592653589793
#define TWOPI 6.283185307179586
