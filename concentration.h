#ifndef CONC_H
#define CONC_H

#include "common.h"


//the concentration object

class concentration{
	protected:
		double base_polym_rate;
		double base_depolym_rate;
/*		double hydrolysis_rate;//in principle, we can have a two-step depolym process
		double cap_rate;
		double uncap_rate; //rate of uncapping*/
		double base_cat_rate;//rate of switching from growth to catastrophe
		double base_res_rate;//rate of switching from catastrophe to rescue
		double base_att_rate;
		double base_det_rate;
		double att_exp_const;
		double det_exp_const;
		double polym_rate;
		double depolym_rate;
		double cat_rate;
		double res_rate;
		double att_rate;
		double det_rate;		

		double att_length;

//		double prev_tension;
		double tension;
		double kt_mt_dist;//special for attachment rate		

		double prev_conc;
		double conc;
		
//		double prev_frac_att;
		double frac_att;

	
	public:
		concentration()
		{
			base_polym_rate = BASE_POLYM_RATE;
			base_depolym_rate = BASE_DEPOLYM_RATE;
			base_res_rate = BASE_RES_RATE;
			base_cat_rate = BASE_CAT_RATE;
			base_att_rate = BASE_ATT_RATE_RESCUE;
			base_det_rate = BASE_DET_RATE_RESCUE;

			att_exp_const = ATT_RATE_EXP_FACTOR_RESCUE;
			det_exp_const = DET_RATE_EXP_FACTOR_RESCUE;

			polym_rate = base_polym_rate;
			depolym_rate = base_depolym_rate;
			res_rate = base_res_rate;
			cat_rate = base_cat_rate;
                        att_rate = base_att_rate;
                        det_rate = base_det_rate;

			att_length = sqrt(2.*EPSILON / KT_MT_SPRING); 

			tension = 0.;
//			prev_tension = tension;
			kt_mt_dist = 0.;

                        conc = 0.;
			prev_conc = conc;

			frac_att = 0.5;			
//			prev_frac_att = frac_att;

		}

//almost the same as default, but this has an initial conc input
		concentration(double init_conc, bool growing)
		{
                        base_polym_rate = BASE_POLYM_RATE;
                        base_depolym_rate = BASE_DEPOLYM_RATE;
                        base_res_rate = BASE_RES_RATE;
                        base_cat_rate = BASE_CAT_RATE;
			if(growing)
			{
                        	base_att_rate = BASE_ATT_RATE_RESCUE;
                        	base_det_rate = BASE_DET_RATE_RESCUE;
	                        att_exp_const = ATT_RATE_EXP_FACTOR_RESCUE;
        	                det_exp_const = DET_RATE_EXP_FACTOR_RESCUE;
			}
			else
			{
				base_att_rate = BASE_ATT_RATE_CAT;
				base_det_rate = BASE_DET_RATE_CAT;
                                att_exp_const = ATT_RATE_EXP_FACTOR_CAT;
                                det_exp_const = DET_RATE_EXP_FACTOR_CAT;
			}
                        polym_rate = base_polym_rate;
                        depolym_rate = base_depolym_rate;
                        res_rate = base_res_rate;
                        cat_rate = base_cat_rate;
                        att_rate = base_att_rate;
                        det_rate = base_det_rate;

                        att_length = sqrt(2.*EPSILON / KT_MT_SPRING);

                        tension = 0.;
//                      prev_tension = tension;
			kt_mt_dist = 0.;

                        conc = 0.;
                        prev_conc = conc;

                        frac_att = 0.5;
//                      prev_frac_att = frac_att;

			prev_conc = init_conc;
			conc = prev_conc;
		}
	
//functions to get protected members
		double get_base_polym_rate() const {return base_polym_rate;}
		double get_base_depolym_rate() const {return base_depolym_rate;}
		double get_base_cat_rate() const {return base_cat_rate;}
		double get_base_res_rate() const {return base_res_rate;}
		double get_base_att_rate() const {return base_att_rate;}
		double get_base_det_rate() const {return base_det_rate;}
		double get_polym_rate() const {return polym_rate;}
		double get_depolym_rate() const {return depolym_rate;}
		double get_cat_rate() const {return cat_rate;}
		double get_res_rate() const {return res_rate;}
//		double get_prev_tension() const {return prev_tension;}
		double get_tension() const {return tension;}
		double get_kt_mt_dist() const {return kt_mt_dist;}
		double get_prev_conc() const {return prev_conc;}
		double get_conc() const {return conc;}
//		double get_prev_frac_att() const {return prev_frac_att;}
		double get_frac_att() const {return frac_att;}
		double get_att_rate() const {return att_rate;}
		double get_det_rate() const {return det_rate;}

//methods to set members
		void set_tension(double force){tension = force;}
		void set_kt_mt_dist(double dist){kt_mt_dist = dist;}
		void set_conc(double conc_in){conc = conc_in; if(conc < 0.) conc = 0.;}
		void set_polym_rate(double rate){polym_rate = rate;}
		void set_depolym_rate(double rate){depolym_rate = rate;}
		void set_cat_rate(double rate){cat_rate = rate;}
		void set_res_rate(double rate){res_rate = rate;}
		void set_frac_att(double frac){frac_att = frac;}
		void set_att_rate(double rate){att_rate = rate;}
		void set_det_rate(double rate){det_rate = rate;}

//set_prev_
		void set_prev_conc()
		{
			if(conc < 0.)
				conc = 0.;
			prev_conc = conc;
		}
	
/**************************************************************/
//Need function(s) to adjust rates based on tension
                void calculate_polym_rate();
                void calculate_depolym_rate();
/*                double calculate_cap_rate();
                double calculate_hydrolysis_rate();
                double calculate_uncap_rate();*/
				void calculate_res_rate();
				void calculate_cat_rate();
				
				void calculate_frac_att();

		void calculate_att_rate();
		void calculate_det_rate();

/*************************************************************/

//functions for kinetics
		double polymerize(double free_conc, double& curr_free_conc, double shorter_fil_polym, bool boundary_bool);//returns amt polymed from this length of filaments
		void depolymerize(double& curr_free_conc, concentration const& longer_fil, bool boundary_bool);
                void rescue_switch(concentration const& cat_fil);//switching terms for filaments undergoing rescue
		void catastrophe_switch(concentration const& res_fil);//switching terms for filaments undergoing catastrophe

		void attach_onoff(concentration const& det_fil);
		void detach_onoff(concentration const& att_fil);
};


extern vector<concentration> left_res;
extern vector<concentration> left_cat;
#if FOUR_POPULATIONS
extern vector<concentration> left_res_det;
extern vector<concentration> left_cat_det;

#if TWO_KINETO
extern vector<concentration> right_res;
extern vector<concentration> right_cat;
extern vector<concentration> right_res_det;
extern vector<concentration> right_cat_det;
#endif
#endif
extern concentration null_conc;





#endif//ifndef 

