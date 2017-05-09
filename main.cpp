#include "common.h"
#include "protos2.h"//has concentration.h in it
#include "kinetochore.h"

#include "concentration_methods.h"

double free_mono;
double prev_free_mono, total_conc;

concentration null_conc; 
vector<concentration> left_res;
vector<concentration> left_cat;
#if FOUR_POPULATIONS
vector<concentration> left_res_det;
vector<concentration> left_cat_det;
#if TWO_KINETO
vector<concentration> right_res;
vector<concentration> right_cat;
vector<concentration> right_res_det;
vector<concentration> right_cat_det;
#endif
#endif
vector<kinetochore> kt_list;

unsigned long long sim_step;

FILE *ktposfile;
FILE *mtlengthfile;
FILE *rightmtlengthfile;
FILE *tensionfile;
FILE *distribfile;
char ktposfilename[128];
char mtlengthfilename[128];
char rightmtlengthfilename[128];
char tensionname[128];
char distribname[128];


double zero_tension_binding_frac;
double minimum_tension_for_det_mod;

int shift_counter;

//max force test
bool spring_started;
double spring_location;

int comm_line_trial;
double comm_line_time;
long comm_line_numsteps;
double comm_line_free_mono;
double comm_line_depolym;
double comm_line_polym;
double comm_line_cat;
double comm_line_rescue;
double comm_line_p_exp, comm_line_comp_p_exp;
double comm_line_d_exp, comm_line_comp_d_exp;
double comm_line_r_exp, comm_line_comp_r_exp;
double comm_line_c_exp, comm_line_comp_c_exp;
double comm_line_a_exp, comm_line_det_exp, comm_line_det_exp_cat;
double comm_line_att_rate_rescue;
double comm_line_det_rate_rescue;
double comm_line_det_rate_cat_factor;
int comm_line_init_fil_length;
int comm_line_lz, comm_line_numskip;
double comm_line_load;
int comm_line_restart;
double comm_line_kt_mt_spring;
double comm_line_kt_kt_spring;
double comm_line_epsilon;
double comm_line_left_r0;
double comm_line_attachment_range;
double comm_line_tdelay;
bool comm_line_dash_pot;
double comm_line_relative_velocity_drag; 
bool comm_line_sigmoidal_c, comm_line_sigmoidal_p;
double comm_line_sig_cat_coeff, comm_line_sig_polym_coeff, comm_line_sig_depolym_coeff, comm_line_sig_res_coeff;
bool comm_line_max_force_test;
double comm_line_max_force_test_spring;

int main(int argc, char **argv)
{

	int option;

	initialize_command_line_variables();
//only some capital letters left
while( (option = getopt(argc, argv, "t:e:n:s:l:f:i:d:p:c:r:w:x:y:z:m:k:a:b:v:o:u:q:j:g:h:U:M:P:S:I:G:D:C:R:W:X:Y:Z:K:T:")) != -1)
{
 switch(option)
 {
  case 't':
    comm_line_trial = atoi(optarg); break;
  case 'e':
    comm_line_restart = atoi(optarg); break;
  case 'n':
    comm_line_time = atof(optarg);
    comm_line_numsteps = comm_line_time / dt;
    break;
//    comm_line_numsteps = atoi(optarg); break;
  case 's':
    comm_line_numskip = atoi(optarg); break;
  case 'l':
    comm_line_lz = atoi(optarg); break;
    if(comm_line_lz > 2500.)
	comm_line_tdelay = comm_line_lz / 2000. * 100.;
  case 'f':
    comm_line_free_mono = (double) atof(optarg); break;
  case 'i':
    #if TWO_KINETO
    fprintf(stderr, "Initial filament length is set by LZ and KT_KT_SPRING_LENGTH. Command line input ignored.\n");
    #else
    comm_line_init_fil_length = atoi(optarg); 
    #endif
    break;
  case 'd':
    comm_line_depolym = (double) atof(optarg); break;
  case 'p':
    comm_line_polym = (double) atof(optarg); break;
  case 'c':
    comm_line_cat = (double) atof(optarg); break;
  case 'r':
    comm_line_rescue = (double) atof(optarg); break;
  case 'q':
    comm_line_det_exp = (double) atof(optarg); break;
  case 'j':
    comm_line_det_exp_cat = (double) atof(optarg); break;
  case 'v':
    comm_line_a_exp = (double) atof(optarg); break;
  case 'w':
    comm_line_p_exp = (double) atof(optarg); break;
  case 'W':
    comm_line_comp_p_exp = (double) atof(optarg); break;
  case 'x':
    comm_line_d_exp = (double) atof(optarg); break;
  case 'X':
    comm_line_comp_d_exp = (double) atof(optarg); break;
  case 'y':
    comm_line_r_exp = (double) atof(optarg); break;
  case 'z':
    comm_line_c_exp = (double) atof(optarg); break;
  case 'Y': 
    comm_line_comp_r_exp = (double) atof(optarg); break;
  case 'Z':
    comm_line_comp_c_exp = (double) atof(optarg); break;
  case 'm':
    comm_line_load = (double) atof(optarg); break;
  case 'k':
    comm_line_kt_mt_spring = (double) atof(optarg); break;
  case 'h':
    comm_line_kt_kt_spring = (double) atof(optarg); break;
  case 'a':
    comm_line_epsilon = (double) atof(optarg); break;
  case 'b':
    comm_line_left_r0 = (double) atof(optarg); break;
  case 'o':
    comm_line_att_rate_rescue = (double) atof(optarg); break;
  case 'u':
    comm_line_det_rate_rescue = (double) atof(optarg); break;
  case 'U':
    comm_line_det_rate_cat_factor = (double) atof(optarg); break;
  case 'g':
    comm_line_attachment_range = (double) atof(optarg); break;
  case 'D':
    if(atoi(optarg) == 1)
      comm_line_dash_pot = true; 
    break;
  case 'P':
    comm_line_relative_velocity_drag = atof(optarg);
    break;
  case 'S':
    if(atoi(optarg) == 1)
     comm_line_sigmoidal_c = true; 
    break;
  case 'I':
    if(atoi(optarg) == 1)
     comm_line_sigmoidal_p = true;
    break;
  case 'C':
    comm_line_sig_cat_coeff = atof(optarg);
    break;
  case 'R':
    comm_line_sig_res_coeff = atof(optarg);
    break;
  case 'G':
    comm_line_sig_polym_coeff = atof(optarg);
    break;
  case 'M':
    comm_line_sig_depolym_coeff = atof(optarg);
    break;
  case 'T':
    if(atoi(optarg) == 1)
    comm_line_max_force_test = true;
    break;
  case 'K':
    comm_line_max_force_test_spring = atof(optarg);
    break;

  case '?':
    cerr <<  "Unknown option character '" << optopt << "'. ...exiting.\n";
    cerr << "Options are: \n";
    cerr << "-t trial number\n";
    cerr << "-e restart\n";
    cerr << "-n length of sim (in sim tau)\n";//number of time steps\n";
    cerr << "-s number of steps to skip between printing data\n";
    cerr << "-l length of box\n";
    cerr << "-f initial number of free monos\n";
    cerr << "-i initial length of filaments\n";
    cerr << "-d depolymerization rate\n";
    cerr << "-p polymerization probability\n";
    cerr << "-c catastrophe rate\n";
    cerr << "-r rescue rate\n";
    cerr << "-o attachment rate\n";
    cerr << "-u detachment rate\n";
    cerr << "-q detachment exp\n";
    cerr << "-j detachment exp for catastrophe mode\n";
    cerr << "-v attachment exp\n";
    cerr << "-w polym exp\n";
    cerr << "-x depolym exp\n";
    cerr << "-y rescue exp\n";
    cerr << "-z catastrophe exp\n";
    cerr << "-m pulling load\n";
    cerr << "-k kt-mt spring const\n";
    cerr << "-h kt-kt spring const\n";
    cerr << "-a epsilon\n";
    cerr << "-b initial conc of left rescue filaments\n";
    cerr << "-g attachment range\n";
    cerr << "-U detachment rate during shrinking mode\n";
    return 1;

  default:
	abort();
 }
}

fprintf(stderr, "trial t %i restart e %i time n %li numsteps s %i box l %i monos f %g fil length i %i depoly d %g poly p %g cat c %g res r %g det exp q %g att exp v %g poly exp w %g d exp x %g res exp y %g cat exp z %g load m %g att rate o %g det rate u %g\n", comm_line_trial, comm_line_restart, comm_line_numsteps, comm_line_numskip, comm_line_lz, comm_line_free_mono, comm_line_init_fil_length, comm_line_depolym, comm_line_polym, comm_line_cat, comm_line_rescue, comm_line_det_exp, comm_line_a_exp, comm_line_p_exp, comm_line_d_exp, comm_line_r_exp, comm_line_c_exp, comm_line_load, comm_line_att_rate_rescue, comm_line_det_rate_rescue);


	if(!RESTART)
	 initialize();
	else
	 read_restart();

char cpcommand[96];
char mvcommand[96];
sprintf(cpcommand, "cp /var/tmp/new_and_improved_scratch/master_eq/output/*%6.6i output/", TRIALNUMBER);
sprintf(mvcommand, "mv /var/tmp/new_and_improved_scratch/master_eq/output/*%6.6i output/", TRIALNUMBER);

	for(sim_step = 0; sim_step < NUMSTEPS+1; sim_step++)
	{
/*
if(sim_step == 0)
	print_mt_data();
*/

		if(sim_step%(NUMSKIP*20) == 0)
			fprintf(stderr, "Step %i\n", sim_step);

		evolve_mts();//growth and shrinkage of mts


		
//		#if !TENSION_INDEPENDENT
		evolve_kts(sim_step);//pulling by tension, recoil due to MTs (later: kt-kt interactions, diffusion)
//		#endif


		update_system();//update "previous" variables to be "current"


		#if INFINITE_BOX
		if(sim_step%(4*NUMSKIP))
			check_kinetochore_position();
		#endif

		if(sim_step%NUMSKIP == 0)
		{
                 print_mt_data();
//		 if(sim_step%(10*NUMSKIP) == 0)
			print_kt_data();
#if !PERTURB_RESTART
		 if((sim_step != 0)&&(sim_step%(100*NUMSKIP) == 0))
		{
#else
		if((sim_step%(2*NUMSKIP) == 0)&&(sim_step != 0))
		{
		char newrestart[96];
		sprintf(newrestart, "mv restart/restart%6.6i restart/restart%6.6i_%8.8li", TRIALNUMBER, TRIALNUMBER, sim_step-NUMSKIP*2);
		system(newrestart);
#endif
		write_restart();
		 if(sim_step%(500*NUMSKIP) == 0)
		 {
		   system(cpcommand);
		 }
		 }//end if(sim_step... 
		

		}//end if(sim_step%NUMSKIP==0)

/*int wait;
scanf("%i", &wait);
fprintf(stderr, "hello\n");
*/
	}//end for loop over steps

	
fflush(tensionfile);
fclose(tensionfile);
fflush(mtlengthfile);
fclose(mtlengthfile);
fflush(ktposfile);
fclose(ktposfile);
fflush(distribfile);
fclose(distribfile);


write_restart();

system(mvcommand);

fprintf(stderr, "Trial number %i completed.\n", TRIALNUMBER);

return 0;

} 
