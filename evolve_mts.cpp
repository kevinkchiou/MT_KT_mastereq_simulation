#include "common.h"
#include "protos2.h"//has concentration.h in it
#include "kinetochore.h"

void evolve_mts()
{


	calc_rates();


/*
for(int ii = 0; ii < left_res.size(); ii++)
{
	fprintf(stderr, "ii= %i, kt-ii= %g, k+= %g %g k-= %g %g kr= %g %g kc= %g %g wa= %g %g wd= %g %g\n", ii, (double)ii-kt_list[0].get_prev_pos(), left_res[ii].get_polym_rate(), left_res_det[ii].get_polym_rate(), left_cat[ii].get_depolym_rate(), left_cat_det[ii].get_depolym_rate(), left_cat[ii].get_res_rate(), left_cat_det[ii].get_res_rate(), left_res[ii].get_cat_rate(), left_res_det[ii].get_cat_rate(), left_res_det[ii].get_att_rate(), left_cat_det[ii].get_att_rate(), left_res[ii].get_det_rate(), left_cat[ii].get_det_rate());
}
exit(1);
*/

	kinetics(left_res, left_cat);
	#if FOUR_POPULATIONS

	kinetics(left_res_det, left_cat_det);

	attachment_dynamics(left_res, left_res_det, left_cat, left_cat_det);


	#if TWO_KINETO
	kinetics(right_res, right_cat);
	kinetics(right_res_det, right_cat_det);
	attachment_dynamics(right_res, right_res_det, right_cat, right_cat_det);
	#endif	
	#endif	
	//kinetics(right_res, right_cat);
}



void kinetics(vector<concentration>& growing, vector<concentration>& shrinking)
{
int ii;
double polym_from_shorter_filaments = 0.;//this variable will tell longer filaments how much polymed into them from shorter fil
int max_ii = growing.size()-1;


//fprintf(stderr, "%i\n", 0);

		polym_from_shorter_filaments = growing[0].polymerize(prev_free_mono, free_mono, polym_from_shorter_filaments, false);//polym terms
		//polym boundary bool is for LZ-1 
		growing[0].rescue_switch(shrinking[0]);//rescue & cat terms for filaments that are growing
		shrinking[0].depolymerize(free_mono, shrinking[1], true);//depolym terms
		shrinking[0].catastrophe_switch(growing[0]);//rescue & cat terms for filaments that are shrinking

for(ii = 1; ii < max_ii; ii++)
{
//fprintf(stderr, "%i", ii);

	polym_from_shorter_filaments = growing[ii].polymerize(prev_free_mono, free_mono, polym_from_shorter_filaments, false);
	growing[ii].rescue_switch(shrinking[ii]);
	shrinking[ii].depolymerize(free_mono, shrinking[ii+1], false);
	shrinking[ii].catastrophe_switch(growing[ii]);

/*
if(!finite(growing[ii].get_conc()))
{
	fprintf(stderr, "ii=%i, k+[ii-1]=%g, f[ii-1]=%g, k[ii]=%g, f[ii-1]=%g, base[ii-1]=%g, prev[ii-1]=%g\n", ii, growing[ii-1].get_polym_rate(), growing[ii-1].get_frac_att(), growing[ii].get_polym_rate(), growing[ii].get_frac_att(), growing[ii-1].get_base_polym_rate(), growing[ii-1].get_prev_conc());
	exit(1);
}
	*/
}

//fprintf(stderr, "%i", 2000);
                polym_from_shorter_filaments = growing[max_ii].polymerize(prev_free_mono, free_mono, polym_from_shorter_filaments, true);
                growing[max_ii].rescue_switch(shrinking[max_ii]);
                shrinking[max_ii].depolymerize(free_mono, null_conc, false);
                shrinking[max_ii].catastrophe_switch(growing[max_ii]);


/*fprintf(stderr, "old and new profiles\n");
for(int ii1 = 0; ii1 < 2001; ii1++)
fprintf(stderr, "%i %8.6g %12.10g %12.10g %12.10g %12.10g\n", ii1, kt_list[0].get_prev_pos(), growing[ii1].get_prev_conc(), shrinking[ii1].get_prev_conc(), growing[ii1].get_conc(), shrinking[ii1].get_conc());

exit(1);
*/
}


void attachment_dynamics(vector<concentration>& att1, vector<concentration>& det1, vector<concentration>& att2, vector<concentration>& det2)
{
	int ii;
	int max_ii = att1.size();

	for(ii = 0; ii < max_ii; ii++)
	{
		att1[ii].attach_onoff(det1[ii]);
		det1[ii].detach_onoff(att1[ii]);
		att2[ii].attach_onoff(det2[ii]);
		det2[ii].detach_onoff(att2[ii]);
	}
}




void calc_rates()
{
	int ii; 

//first calculate previous tension
	double pos_diff, pos_diff_2;
	double force;
	int max_ii = left_res.size();
	
	for(ii = 0; ii < max_ii; ii++)
	{
	  pos_diff = kt_list[0].get_prev_pos() - ii*MONO_DIAM;//note use of previous position of kt
	  force = KT_MT_SPRING*pos_diff;//if pos_diff > 0, tension; else compression
	  left_res[ii].set_tension(force);
	  left_cat[ii].set_tension(force);
	  left_res[ii].set_kt_mt_dist(pos_diff);
	  left_cat[ii].set_kt_mt_dist(pos_diff);
	  #if FOUR_POPULATIONS
	  if(pos_diff > 0.)
	  {
	    left_res_det[ii].set_tension(0.);
	    left_cat_det[ii].set_tension(0.);
	  }
	  else
	  {
	    left_res_det[ii].set_tension(force);
	    left_cat_det[ii].set_tension(force);
	  }//filaments still feel compression when detached
	  left_res_det[ii].set_kt_mt_dist(pos_diff);
	  left_cat_det[ii].set_kt_mt_dist(pos_diff);
	  #if TWO_KINETO
	  pos_diff_2 = kt_list[1].get_prev_pos() - ii*MONO_DIAM;
	  force = KT_MT_SPRING*pos_diff_2;
	  right_res[ii].set_tension(force);
	  right_cat[ii].set_tension(force);
	  right_res[ii].set_kt_mt_dist(pos_diff_2);
	  right_cat[ii].set_kt_mt_dist(pos_diff_2);
          if(pos_diff_2 > 0.)
          {
            right_res_det[ii].set_tension(0.);
            right_cat_det[ii].set_tension(0.);
          }//recall right side not measured in absolute position, rather distance from LZ
          else
          {
            right_res_det[ii].set_tension(force);
            right_cat_det[ii].set_tension(force);
          }//filaments still feel compression when detached
          right_res_det[ii].set_kt_mt_dist(pos_diff_2);
          right_cat_det[ii].set_kt_mt_dist(pos_diff_2);
	  #endif
	  #endif


#if !FOUR_POPULATIONS
//then calculate tension-dependent rates and frac_att
		left_res[ii].calculate_frac_att();
		left_cat[ii].set_frac_att(left_res[ii].get_frac_att());
//	fprintf(stderr, "%g %g\n", force, left_res[ii].get_frac_att());
//if fraction attached is different for recovering and shrinking filaments, need to explicitly calculate both
//		left_cat[ii].calc_frac_att();

#endif


#if !FOUR_POPULATIONS
//note the variables polym_rate, etc are only used for the attached fraction of filaments.
//base_polym_rate, etc. are used for unattached filaments.      
//need these lines so exponentially large tension doesn't cause rates to blow up
                if(left_res[ii].get_frac_att() > 1.e-12)
                {
#endif	
//for all sims, calculate polym/depolym/res/cat rates
			left_res[ii].calculate_polym_rate();
			left_res[ii].calculate_cat_rate();
			left_cat[ii].calculate_depolym_rate();
			left_cat[ii].calculate_res_rate();
#if FOUR_POPULATIONS
			if(pos_diff < 0.)
			{
                         left_res_det[ii].calculate_polym_rate();
                         left_res_det[ii].calculate_cat_rate();
                         left_cat_det[ii].calculate_depolym_rate();
                         left_cat_det[ii].calculate_res_rate();
			}
			else
			{
                         left_res_det[ii].set_polym_rate(left_res_det[ii].get_base_polym_rate());
                         left_res_det[ii].set_cat_rate(left_res_det[ii].get_base_res_rate());
                         left_cat_det[ii].set_depolym_rate(left_cat_det[ii].get_base_depolym_rate());
                         left_cat_det[ii].set_res_rate(left_cat_det[ii].get_base_cat_rate());
			}

#endif

#if !FOUR_POPULATIONS
		}
                else
                {
//prevent numerical error
                        left_res[ii].set_polym_rate(0.);
                        left_res[ii].set_cat_rate(0.);
                        left_cat[ii].set_depolym_rate(0.);
                        left_cat[ii].set_res_rate(0.);
                }
#endif

			#if FOUR_POPULATIONS
//attachment/detachment rates in four populations model
			//detached filaments always have base rate

			left_res[ii].calculate_det_rate();
			left_res_det[ii].calculate_att_rate();
			left_cat[ii].calculate_det_rate();
			left_cat_det[ii].calculate_att_rate();

//fprintf(stderr, "i=%i x-i=%g wa=%g\n", ii, pos_diff, left_res_det[ii].get_att_rate());

			#if TWO_KINETO
//calculate everything for right filaments if we have two kinetochores
			right_res[ii].calculate_polym_rate();
			right_res[ii].calculate_cat_rate();
                        right_cat[ii].calculate_depolym_rate();
                        right_cat[ii].calculate_res_rate();
                        right_res[ii].calculate_det_rate();
                        right_res_det[ii].calculate_att_rate();
                        right_cat[ii].calculate_det_rate();
                        right_cat_det[ii].calculate_att_rate();

			if(pos_diff_2 < 0.)
			{
                         right_res_det[ii].calculate_polym_rate();
                         right_res_det[ii].calculate_cat_rate();
                         right_cat_det[ii].calculate_depolym_rate();
                         right_cat_det[ii].calculate_res_rate();
                        }
                        else
                        {
                         right_res_det[ii].set_polym_rate(right_res_det[ii].get_base_polym_rate());
                         right_res_det[ii].set_cat_rate(right_res_det[ii].get_base_res_rate());
                         right_cat_det[ii].set_depolym_rate(right_cat_det[ii].get_base_depolym_rate());
                         right_cat_det[ii].set_res_rate(right_cat_det[ii].get_base_cat_rate());
                        }

			#endif
			#endif
	}
}
