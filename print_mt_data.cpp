#include "common.h"
#include "concentration.h"
#include "kinetochore.h"

void print_mt_data()
{
	print_mean_mt_length();
	print_tension();
	print_distrib_stats();
}

void print_mean_mt_length()
{
	int ii;

	double mean_c_length = 0.;
	double mean_r_length = 0.;
	double mean_rd_length = 0.;
	double mean_cd_length = 0.;
	double mean_all_length = 0.;

	double mean_right_c_length = 0.;
	double mean_right_r_length = 0.;
	double mean_right_rd_length = 0.;
	double mean_right_cd_length = 0.;
	double mean_right_all_length = 0.;	

	double tot_r_conc = 0.;
	double tot_c_conc = 0.;
	double tot_rd_conc = 0.;
	double tot_cd_conc = 0.;

	double tot_right_r_conc = 0.;
	double tot_right_c_conc = 0.;
	double tot_right_rd_conc = 0.;
	double tot_right_cd_conc = 0.;

	int max_ii =left_res.size();
	double tot_att = 0.;
	double tot_right_att = 0.;

	for(ii = 0; ii < max_ii; ii++)
	{
		mean_r_length += ii*left_res[ii].get_prev_conc();
		mean_c_length += ii*left_cat[ii].get_prev_conc();
		#if FOUR_POPULATIONS
		mean_rd_length += ii*left_res_det[ii].get_prev_conc();
		mean_cd_length += ii*left_cat_det[ii].get_prev_conc();

		#if TWO_KINETO
		mean_right_r_length += ii*right_res[ii].get_prev_conc();
		mean_right_c_length += ii*right_cat[ii].get_prev_conc();
		mean_right_rd_length += ii*right_res_det[ii].get_prev_conc();
		mean_right_cd_length += ii*right_cat_det[ii].get_prev_conc();
		#endif
		#endif

		tot_r_conc += left_res[ii].get_prev_conc();
		tot_c_conc += left_cat[ii].get_prev_conc();
		#if FOUR_POPULATIONS
		tot_rd_conc += left_res_det[ii].get_prev_conc();
		tot_cd_conc += left_cat_det[ii].get_prev_conc();
		
		#if TWO_KINETO
		tot_right_r_conc += right_res[ii].get_prev_conc();
		tot_right_c_conc += right_cat[ii].get_prev_conc();
		tot_right_rd_conc += right_res_det[ii].get_prev_conc();
		tot_right_cd_conc += right_cat_det[ii].get_prev_conc();
		#endif
		#endif

		#if !FOUR_POPULATIONS
		tot_att += left_res[ii].get_prev_conc()*left_res[ii].get_frac_att() + left_cat[ii].get_prev_conc()*left_cat[ii].get_frac_att();
		#else//just 4pop
		tot_att += left_res[ii].get_prev_conc() + left_cat[ii].get_prev_conc();
		#if TWO_KINETO
		tot_right_att += right_res[ii].get_prev_conc() + right_cat[ii].get_prev_conc();
		#endif
		#endif
	}




	#if !FOUR_POPULATIONS
        mean_all_length = (mean_r_length + mean_c_length) / (tot_r_conc + tot_c_conc);
	#else
	mean_all_length = (mean_r_length + mean_c_length + mean_rd_length + mean_cd_length) / (tot_r_conc + tot_c_conc + tot_rd_conc + tot_cd_conc);
	#if TWO_KINETO
	mean_right_all_length = (mean_right_r_length + mean_right_c_length + mean_right_rd_length + mean_right_cd_length) / (tot_right_r_conc + tot_right_c_conc + tot_right_rd_conc + tot_right_cd_conc);
	#endif
	#endif
	mean_r_length *= 1./tot_r_conc;
	mean_c_length *= 1./tot_c_conc;
	#if FOUR_POPULATIONS
	mean_rd_length *= 1./tot_rd_conc;
	mean_cd_length *= 1./tot_cd_conc;
	#if TWO_KINETO
	mean_right_r_length *= 1./tot_right_r_conc;
	mean_right_c_length *= 1./tot_right_c_conc;
	mean_right_rd_length *= 1./tot_right_rd_conc;
	mean_right_cd_length *= 1./tot_right_cd_conc;
	#endif
	#endif	

	double shift = SHIFT_LENGTH*shift_counter;
	#if FOUR_POPULATIONS
	fprintf(mtlengthfile, "%g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g\n", sim_step*dt, mean_all_length + shift, mean_r_length + shift, mean_c_length + shift, mean_rd_length + shift, mean_cd_length + shift, tot_r_conc, tot_c_conc, tot_rd_conc, tot_cd_conc, prev_free_mono, tot_att);
	#if TWO_KINETO
	fprintf(rightmtlengthfile, "%g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g\n", sim_step*dt, mean_right_all_length - shift, mean_right_r_length - shift, mean_right_c_length - shift, mean_right_rd_length - shift, mean_right_cd_length - shift, tot_right_r_conc, tot_right_c_conc, tot_right_rd_conc, tot_right_cd_conc, prev_free_mono, tot_right_att);
	fflush(rightmtlengthfile);
	#endif
	#else
        fprintf(mtlengthfile, "%g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g %12.10g\n", sim_step*dt, mean_all_length + shift, mean_r_length + shift, mean_c_length + shift, tot_r_conc, tot_c_conc, prev_free_mono, tot_att);
	#endif
	fflush(mtlengthfile);
}


void print_tension()
{
	int ii;
        int max_ii =left_res.size();

	double tension = 0.;
	double compression = 0.;
	double conc_tension = 0.;
	double conc_compression = 0.;


	for(ii = 0; ii < max_ii; ii++)
	{
		if(ii < kt_list[0].get_pos())
		{
			tension += (ii-kt_list[0].get_pos())*left_res[ii].get_conc();
			tension += (ii-kt_list[0].get_pos())*left_cat[ii].get_conc();
			conc_tension += left_res[ii].get_conc() + left_cat[ii].get_conc();
		}
		else
		{
			compression += (ii-kt_list[0].get_pos())*left_res[ii].get_conc();
			compression += (ii-kt_list[0].get_pos())*left_cat[ii].get_conc();
                        compression += (ii-kt_list[0].get_pos())*left_res_det[ii].get_conc();
                        compression += (ii-kt_list[0].get_pos())*left_cat_det[ii].get_conc();
			conc_compression += left_res[ii].get_conc() + left_cat[ii].get_conc();
			conc_compression += left_res_det[ii].get_conc() + left_cat_det[ii].get_conc();
		}
		
	} 
	
	tension *= 1./conc_tension;
	compression *= 1./conc_compression;

if(finite(tension))
{
	fprintf(tensionfile, "%g %g %g %g %g\n", sim_step*dt, tension, compression, conc_tension, conc_compression);//note that this is all reduced by a factor of spring constant

//        fprintf(stderr, "%g %g %g %g %g\n", sim_step*dt, tension, compression, conc_tension, conc_compression);

	fflush(tensionfile);
}
else
	fprintf(stderr, "tension not finite\n");
}



void print_distrib_stats()
{

        double g_mt_pos_avg = 0.;
        double s_mt_pos_avg = 0.;
        double g_mt_pos_var = 0.;
        double s_mt_pos_var = 0.;
        double g_tot_conc = 0.;
        double s_tot_conc = 0.;

#if TWO_KINETO
        double right_g_mt_pos_avg = 0.;
        double right_s_mt_pos_avg = 0.;
        double right_g_mt_pos_var = 0.;
        double right_s_mt_pos_var = 0.;
        double right_g_tot_conc = 0.;
        double right_s_tot_conc = 0.;
#endif

	int ii; 
	int max_ii = left_res.size();

	for(ii = 0; ii < max_ii; ii++)
	{
		g_mt_pos_avg += (ii - kt_list[0].get_prev_pos())*left_res[ii].get_prev_conc();
		s_mt_pos_avg += (ii - kt_list[0].get_prev_pos())*left_cat[ii].get_prev_conc();
                g_mt_pos_var += (ii - kt_list[0].get_prev_pos())*(ii - kt_list[0].get_prev_pos())*left_res[ii].get_prev_conc();
		s_mt_pos_var += (ii - kt_list[0].get_prev_pos())*(ii - kt_list[0].get_prev_pos())*left_cat[ii].get_prev_conc();
		g_tot_conc += left_res[ii].get_prev_conc();
		s_tot_conc += left_cat[ii].get_prev_conc();
		#if TWO_KINETO
                right_g_mt_pos_avg += (ii - kt_list[1].get_prev_pos())*right_res[ii].get_prev_conc();
                right_s_mt_pos_avg += (ii - kt_list[1].get_prev_pos())*right_cat[ii].get_prev_conc();
                right_g_tot_conc += right_res[ii].get_prev_conc();
                right_s_tot_conc += right_cat[ii].get_prev_conc();
                right_g_mt_pos_var += (ii - kt_list[1].get_prev_pos())*(ii - kt_list[1].get_prev_pos())*right_res[ii].get_prev_conc();
                right_s_mt_pos_var += (ii - kt_list[1].get_prev_pos())*(ii - kt_list[1].get_prev_pos())*right_cat[ii].get_prev_conc();
		#endif
	}	
	
        g_mt_pos_avg *= 1./g_tot_conc;
        s_mt_pos_avg *= 1./s_tot_conc;
        g_mt_pos_var *= 1./g_tot_conc;
	s_mt_pos_var *= 1./s_tot_conc;
        g_mt_pos_var -= g_mt_pos_avg*g_mt_pos_avg;
        s_mt_pos_var -= s_mt_pos_avg*s_mt_pos_avg;

        #if TWO_KINETO
        right_g_mt_pos_avg *= 1./right_g_tot_conc;
        right_s_mt_pos_avg *= 1./right_s_tot_conc;
        right_g_mt_pos_var *= 1./right_g_tot_conc;
        right_s_mt_pos_var *= 1./right_s_tot_conc;
        right_g_mt_pos_var -= right_g_mt_pos_avg*right_g_mt_pos_avg;
        right_s_mt_pos_var -= right_s_mt_pos_avg*right_s_mt_pos_avg;
        #endif


	fprintf(distribfile, "%g %g %g %g %g", sim_step*dt, g_mt_pos_avg, s_mt_pos_avg, g_mt_pos_var, s_mt_pos_var);
	#if TWO_KINETO
	fprintf(distribfile, " %g %g %g %g", right_g_mt_pos_avg, right_s_mt_pos_avg, right_g_mt_pos_var, right_s_mt_pos_var);
	#endif
	fprintf(distribfile, "\n");


fflush(distribfile);
}




