#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

void evolve_kts(unsigned long long step)
{
#if !TWO_KINETO
if(step * dt > TIME_DELAY)
{
if(MAX_FORCE_TEST)
{
	if(spring_started)
	{
		double spring_force = MAX_FORCE_TEST_SPRING*(spring_location - (kt_list[0].get_prev_pos()+SHIFT_LENGTH*shift_counter));
		fprintf(stderr, "spring force = %g, ktpos = %g, spring location = %g\n", spring_force, kt_list[0].get_prev_pos()+SHIFT_LENGTH*shift_counter, spring_location);
		kt_list[0].load(spring_force);
	}
	else
	{
		spring_location = kt_list[0].get_prev_pos()+SHIFT_LENGTH*shift_counter;
		spring_started = true;
		
		fprintf(stderr, "spring loc = %g\n", spring_location);
	}
}
}
if((step * dt > TIME_DELAY) || RESTART)
	kt_list[0].load(PULLING_LOAD);
#endif

	mt_kt_int();

//	kt_kt_int();
}


void mt_kt_int()
{
int ii;
double force;
double force_2;
double kkforce;
int max_ii = left_res.size();

force = 0.;
force_2 = 0.;
//left kt
for(ii = 0; ii < max_ii; ii++)
{
//force = 0.;
//to be added to kt position
#if FOUR_POPULATIONS
	force -= KT_MT_SPRING*(kt_list[0].get_prev_pos() - (double)ii)*(left_res[ii].get_prev_conc()+left_cat[ii].get_prev_conc());
	if((double) ii > kt_list[0].get_prev_pos())//take the kt position to be inner kt boundary
	   force -= KT_MT_SPRING*(kt_list[0].get_prev_pos() - (double)ii)*(left_res_det[ii].get_prev_conc()+left_cat_det[ii].get_prev_conc());
#if TWO_KINETO
	force_2 -= KT_MT_SPRING*(kt_list[1].get_prev_pos() - (double)ii)*(right_res[ii].get_prev_conc()+right_cat[ii].get_prev_conc());
        if((double) ii > kt_list[1].get_prev_pos())//take the kt position to be inner kt boundary
           force_2 -= KT_MT_SPRING*(kt_list[1].get_prev_pos() - (double)ii)*(right_res_det[ii].get_prev_conc()+right_cat_det[ii].get_prev_conc());
#endif
#else
	force -= KT_MT_SPRING*(kt_list[0].get_prev_pos() - (double)ii)*(left_res[ii].get_frac_att()*left_res[ii].get_prev_conc() + left_cat[ii].get_frac_att()*left_cat[ii].get_prev_conc());
	if(left_res[ii].get_tension() < 0.)
		force -= KT_MT_SPRING*(kt_list[0].get_prev_pos() - (double)ii)*((1.0 - left_res[ii].get_frac_att())*left_res[ii].get_prev_conc() + (1.0 - left_cat[ii].get_frac_att())*left_cat[ii].get_prev_conc());
#endif

}


#if TWO_KINETO
kkforce = KT_KT_SPRING*((double)LZ - kt_list[1].get_prev_pos() - kt_list[0].get_prev_pos() - KT_KT_SPRING_LENGTH);//conversion of right kt position coordinates occurs here

if(DASH_POT)
{
	dash_pot(force, force_2, kkforce);
}
else
{
#endif

kt_list[0].load(force);

#if TWO_KINETO
kt_list[1].load(force_2);

/*
fprintf(stderr, "%g %g %g\n", kkforce, force, force_2);
int wait;
scanf("%i", &wait);
*/
//fprintf(stderr, "%g\n", kkforce);
kt_list[0].load(kkforce);//if kkforce>0, spring stretched too much & kt[0] should move to the right
kt_list[1].load(kkforce);//in that case, kt[1] should move left (which corresponds to moving kt[1] in the direction of longer right side filaments, + in its coordinates)

}
#endif

}

