#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

//these are various asymmetric initial conditions for filament concentrations

#define STRETCHING 5 

#if (TWO_KINETO)
void symmetric_initial_condition1()
{

left_res[comm_line_init_fil_length-STRETCHING/2].set_conc(0.5*LEFT_R_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING/2].set_conc(0.5*LEFT_C_0_TOT);
left_res[comm_line_init_fil_length-STRETCHING/2].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING/2].set_prev_conc();
right_res[comm_line_init_fil_length-STRETCHING/2].set_conc(0.5*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length-STRETCHING/2].set_conc(0.5*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length-STRETCHING/2].set_prev_conc();
right_cat[comm_line_init_fil_length-STRETCHING/2].set_prev_conc();

}




void initial_condition1()
{
left_res[comm_line_init_fil_length].set_conc(0.5*LEFT_R_0_TOT);
left_cat[comm_line_init_fil_length].set_conc(0.5*LEFT_C_0_TOT);
left_res[comm_line_init_fil_length].set_prev_conc();
left_cat[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length].set_conc(0.25*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length].set_conc(0.25*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length].set_prev_conc();
right_cat[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length-1].set_conc(0.25*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length-1].set_conc(0.25*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length-1].set_prev_conc();
right_cat[comm_line_init_fil_length-1].set_prev_conc();
}



void initial_condition2()
{
left_res[comm_line_init_fil_length].set_conc(0.5*LEFT_R_0_TOT);
left_cat[comm_line_init_fil_length].set_conc(0.5*LEFT_C_0_TOT);
left_res[comm_line_init_fil_length].set_prev_conc();
left_cat[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length].set_conc(0.150*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length].set_conc(0.150*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length].set_prev_conc();
right_cat[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length-1].set_conc(0.135*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length-1].set_conc(0.135*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length-1].set_prev_conc();
right_cat[comm_line_init_fil_length-1].set_prev_conc();
right_res[comm_line_init_fil_length-2].set_conc(0.115*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length-2].set_conc(0.115*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length-2].set_prev_conc();
right_cat[comm_line_init_fil_length-2].set_prev_conc();
right_res[comm_line_init_fil_length-3].set_conc(0.100*RIGHT_R_0_TOT);
right_cat[comm_line_init_fil_length-3].set_conc(0.100*RIGHT_C_0_TOT);
right_res[comm_line_init_fil_length-3].set_prev_conc();
right_cat[comm_line_init_fil_length-3].set_prev_conc();
}



void initial_condition3()
{
left_cat[comm_line_init_fil_length-STRETCHING+3].set_conc(0.1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+3].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING+2].set_conc(0.1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+2].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING+1].set_conc(0.15*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+1].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING].set_conc(0.15*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING-1].set_conc(0.2*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING-1].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING-2].set_conc(0.15*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING-2].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING-3].set_conc(0.1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING-3].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING-4].set_conc(0.05*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING-4].set_prev_conc();
right_res[comm_line_init_fil_length+4].set_conc(0.05*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+4].set_prev_conc();
right_res[comm_line_init_fil_length+3].set_conc(0.1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+3].set_prev_conc();
right_res[comm_line_init_fil_length+2].set_conc(0.15*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+2].set_prev_conc();
right_res[comm_line_init_fil_length+1].set_conc(0.2*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+1].set_prev_conc();
right_res[comm_line_init_fil_length].set_conc(0.15*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length-1].set_conc(0.15*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-1].set_prev_conc();
right_res[comm_line_init_fil_length-2].set_conc(0.1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-2].set_prev_conc();
right_res[comm_line_init_fil_length-3].set_conc(0.1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-3].set_prev_conc();
}


void initial_condition4()//extracted from a pulling and pushing sim, respectively
{
left_cat[comm_line_init_fil_length-STRETCHING].set_conc(7e-3*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING+1].set_conc(1.16e-1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+1].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING+2].set_conc(4.89e-1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+2].set_prev_conc();
left_cat[comm_line_init_fil_length-STRETCHING+3].set_conc(3.88e-1*LEFT_C_0_TOT);
left_cat[comm_line_init_fil_length-STRETCHING+3].set_prev_conc();

right_res[comm_line_init_fil_length-3].set_conc(1.3e-2*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-3].set_prev_conc();
right_res[comm_line_init_fil_length-2].set_conc(1.17e-1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-2].set_prev_conc();
right_res[comm_line_init_fil_length-1].set_conc(3.43e-1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length-1].set_prev_conc();
right_res[comm_line_init_fil_length].set_conc(3.61e-1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length].set_prev_conc();
right_res[comm_line_init_fil_length+1].set_conc(1.46e-1*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+1].set_prev_conc();
right_res[comm_line_init_fil_length+2].set_conc(2e-2*RIGHT_R_0_TOT);
right_res[comm_line_init_fil_length+2].set_prev_conc();

//fprintf(stderr, "hello\n");

}


#endif//2kinteo


void standard_kt_initial_condition()
{
int ii;
for(ii = 0; ii < NUMBER_OF_KTS; ii++)
{
        kinetochore temp_kt;
        kt_list.push_back(temp_kt);
        kt_list[ii].set_pos(comm_line_init_fil_length);
        kt_list[ii].set_prev_pos();
}
}

void symmetric_stretched_kt_initial_condition()
{
        int ii;
        kinetochore temp_kt;
for(ii = 0; ii < NUMBER_OF_KTS; ii++)
{
        kt_list.push_back(temp_kt);
	int temp = 2*((double)ii - 0.5);//if ii = 0 then temp = -1, if ii=1 then temp = +1
        kt_list[ii].set_pos(comm_line_init_fil_length+temp*STRETCHING/2);
        kt_list[ii].set_prev_pos();
}

}



void stretched_kt_initial_condition()
{
	int ii;
	kinetochore temp_kt;
for(ii = 0; ii < NUMBER_OF_KTS; ii++)
{
	kt_list.push_back(temp_kt);
        kt_list[ii].set_pos(comm_line_init_fil_length+(ii-1)*STRETCHING);
        kt_list[ii].set_prev_pos();
}
}


