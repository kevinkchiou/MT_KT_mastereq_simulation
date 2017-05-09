#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

void update_system()
{
int ii;
double total_fil = 0.;
int max_ii = left_res.size();

	for(ii = 0; ii < kt_list.size(); ii++)
	{
		kt_list[ii].set_prev_pos();
	}
	for(ii = 0; ii < max_ii; ii++)
	{
//		left_res[ii].set_prev_tension();
		left_res[ii].set_prev_conc();
//		left_cat[ii].set_prev_tension();
		left_cat[ii].set_prev_conc();
		#if FOUR_POPULATIONS
		left_res_det[ii].set_prev_conc();
		left_cat_det[ii].set_prev_conc();
		
		#if TWO_KINETO
		right_res[ii].set_prev_conc();
                right_cat[ii].set_prev_conc();
                right_res_det[ii].set_prev_conc();
                right_cat_det[ii].set_prev_conc();
		#endif

		#endif
		

		#if !CONST_FREE_MONO
		total_fil += ii*left_res[ii].get_conc();
		total_fil += ii*left_cat[ii].get_conc();
		#if FOUR_POPULATIONS
		total_fil += ii*left_res_det[ii].get_conc();
		total_fil += ii*left_cat_det[ii].get_conc();
		#if TWO_KINETO
		total_fil += ii*right_res[ii].get_conc();
                total_fil += ii*right_cat[ii].get_conc();
                total_fil += ii*right_res_det[ii].get_conc();
                total_fil += ii*right_cat_det[ii].get_conc();
		#endif
		#endif
		#endif
	}
#if !CONST_FREE_MONO
	set_prev_mono_conc(total_fil);
#endif
}


void set_prev_mono_conc(double filaments)
{
	prev_free_mono = total_conc - filaments;
	free_mono = prev_free_mono;//fix numerical error...
//	prev_free_mono = free_mono;
}
