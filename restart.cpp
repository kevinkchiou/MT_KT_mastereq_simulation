#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

void write_restart()
{
int ii; 
FILE *restartfile;
char restartname[94];

fprintf(stderr, "writing restart\n");
sprintf(restartname, "restart/restart%6.6i", TRIALNUMBER);
restartfile = fopen(restartname, "w");

for(ii = 0; ii < left_res.size(); ii++)
{
	fprintf(restartfile, "%i %12.10g %12.10g ", ii, left_res[ii].get_prev_conc(), left_cat[ii].get_prev_conc());
	#if FOUR_POPULATIONS
	fprintf(restartfile, "%12.10g %12.10g ", left_res_det[ii].get_prev_conc(), left_cat_det[ii].get_prev_conc());
	#if TWO_KINETO
	fprintf(restartfile, "%12.10g %12.10g %12.10g %12.10g ", right_res[ii].get_prev_conc(), right_cat[ii].get_prev_conc(), right_res_det[ii].get_prev_conc(), right_cat_det[ii].get_prev_conc());
	#endif
	#endif
	fprintf(restartfile, "\n");
}

fprintf(restartfile, "%12.10g\n", prev_free_mono);

for(ii = 0; ii < kt_list.size(); ii++)
	fprintf(restartfile, "%i %12.10g\n", ii, kt_list[ii].get_prev_pos()); 

fflush(restartfile);
fclose(restartfile);

}//write_restart


void read_restart()
{
int ii;
FILE *restartfile;
char restartname[94];
int tempid;
double tempconc1, tempconc2, tempconc3, tempconc4, temppos;
double tempconc5, tempconc6, tempconc7, tempconc8;

sprintf(restartname, "restart/restart%6.6i", TRIALNUMBER);
restartfile = fopen(restartname, "r");


for(ii = 0; ii < LZ+1; ii++)
{
        concentration temp_conc_obj1(0., true);
	concentration temp_conc_obj2(0., false);

        left_res.push_back(temp_conc_obj1);
        left_cat.push_back(temp_conc_obj2);

	#if FOUR_POPULATIONS
	left_res_det.push_back(temp_conc_obj1);
	left_cat_det.push_back(temp_conc_obj2);
	#if TWO_KINETO
	right_res.push_back(temp_conc_obj1);
	right_cat.push_back(temp_conc_obj2);
	right_cat_det.push_back(temp_conc_obj1);
	right_res_det.push_back(temp_conc_obj2);
	#endif
	#endif

	fscanf(restartfile, "%i %lg %lg", &tempid, &tempconc1, &tempconc2);
	#if FOUR_POPULATIONS
	fscanf(restartfile, "%lg %lg", &tempconc3, &tempconc4);
	#if TWO_KINETO
	fscanf(restartfile, "%lg %lg %lg %lg", &tempconc5, &tempconc6, &tempconc7, &tempconc8);
	#endif
	#endif
	left_res[ii].set_conc(tempconc1);
	left_res[ii].set_prev_conc();
	left_cat[ii].set_conc(tempconc2);
	left_cat[ii].set_prev_conc();

	#if FOUR_POPULATIONS
	left_res_det[ii].set_conc(tempconc3);
	left_cat_det[ii].set_conc(tempconc4);
	left_res_det[ii].set_prev_conc();
	left_cat_det[ii].set_prev_conc();
	#if TWO_KINETO
	right_res[ii].set_conc(tempconc5);
        right_cat[ii].set_conc(tempconc6);
        right_res_det[ii].set_conc(tempconc7);
        right_cat_det[ii].set_conc(tempconc8);
        right_res[ii].set_prev_conc();
        right_cat[ii].set_prev_conc();
        right_res_det[ii].set_prev_conc();
        right_cat_det[ii].set_prev_conc();
	#endif
	#endif
}

fscanf(restartfile, "%lg", &tempconc1);
prev_free_mono = tempconc1;
free_mono = prev_free_mono;

for(ii = 0; ii < NUMBER_OF_KTS; ii++)
{
        kinetochore temp_kt;
        kt_list.push_back(temp_kt);

	fscanf(restartfile, "%i %lg\n", &tempid, &temppos);
fprintf(stderr, " temp pos is %g \n", temppos);
	kt_list[ii].set_pos(temppos);
	kt_list[ii].set_prev_pos();
}

initialize_files();

minimum_tension_for_det_mod = 0.;
shift_counter = 0;

fprintf(stderr, "size is %i\n", left_res.size());



#if PERTURB_RESTART
perturb();
write_restart();
#endif


//write_restart();
//exit(1);

}//read_restart



void perturb()
{
if(false)
{
int ii;
double temp1, temp2, temp3;


FILE *perturbfile;
char perturbname[96];
double mag, disp, rhoc, rhor;
int pos;

sprintf(perturbname, "perturbation");
perturbfile = fopen(perturbname, "r");

fscanf(perturbfile, "%lg", &mag);
for(ii = 0; ii < 20; ii++)
{
	fscanf(perturbfile, "%i %lg %lg", &pos, &rhor, &rhoc);
	temp1 = left_res[pos].get_conc();
	left_res[pos].set_conc(temp1+rhor*mag);
        temp1 = left_cat[pos].get_conc();
	left_cat[pos].set_conc(temp1+rhoc*mag);
	left_res[pos].set_prev_conc();
	left_cat[pos].set_prev_conc();
}
fscanf(perturbfile, "%lg", &disp);
temp1 = kt_list[0].get_pos();
kt_list[0].set_pos(temp1 + disp*mag);
kt_list[0].set_prev_pos();

fflush(perturbfile);
fclose(perturbfile);


}
/*for(ii = 0; ii < left_res.size(); ii++)
{
 temp1 = left_res[ii].get_conc();
 temp2 = left_cat[ii].get_conc();

 left_res[ii].set_conc(temp2);
 left_res[ii].set_prev_conc();
 left_cat[ii].set_conc(temp1);
 left_cat[ii].set_prev_conc();
}*/

/*
int ktposition = (int) kt_list[0].get_pos();

temp1 = 0.25*left_cat[ktposition].get_conc();
temp2 = left_res[ktposition].get_conc();
temp3 = left_cat[ktposition].get_conc();

left_res[ktposition].set_conc(temp2+temp1);
left_cat[ktposition].set_conc(temp3-temp1);
left_res[ktposition].set_prev_conc();
left_cat[ktposition].set_prev_conc();
*/
/*
#if TWO_KINETO
kt_list[0].set_pos(kt_list[0].get_pos() + 1.0);
kt_list[1].set_pos(kt_list[1].get_pos() - 1.0);
kt_list[0].set_prev_pos();
kt_list[1].set_prev_pos();
#else
fprintf(stderr, "Hey this is the null perturbation!\n");
#endif
*/
fprintf(stderr, "Hey this is the null perturbation!\n");

}



