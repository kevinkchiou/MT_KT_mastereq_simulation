#include "common.h"
#include "concentration.h"

void concentration::calculate_polym_rate()
{
        double base_rate, exp_factor;
        double tension_dep_rate;
        base_rate = base_polym_rate;
if(tension < 0.)
	exp_factor = COMPRESSION_POLYM_EXP_FACTOR;
else
	exp_factor = POLYM_EXP_FACTOR;

if(!SIGMOIDAL_POLYM_RATE)
{
        tension_dep_rate = base_rate * exp(tension*exp_factor);
}
else
{
        tension_dep_rate = base_rate / (1.0 + SIGMOID_POLYM_COEFF*exp(-1.*tension*exp_factor));
}
//this should only occur for cases of extreme tension/compression...
//this statement avoids computer errors due to infinities
//obviously, can't move more than exist out of of current concentration holder
#if !FOUR_POPULATIONS
	if((tension_dep_rate*frac_att + (1-frac_att)*base_rate)*dt*prev_free_mono > 1.)
		tension_dep_rate = (inv_dt - (1-frac_att)*base_rate) / frac_att / prev_free_mono;
#else
	if(tension_dep_rate*dt*prev_free_mono > 1.)
		tension_dep_rate = inv_dt / prev_free_mono;
#endif
	
        polym_rate = tension_dep_rate;

//fprintf(stderr, "tension %g polym %g\n", tension, polym_rate);
}


void concentration::calculate_depolym_rate()
{
        double base_rate, exp_factor;
        double tension_dep_rate;
        base_rate = base_depolym_rate;

if(tension < 0.)
        exp_factor = COMPRESSION_DEPOLYM_EXP_FACTOR;
else
        exp_factor = DEPOLYM_EXP_FACTOR;
        
if(!SIGMOIDAL_DEPOLYM_RATE)
        tension_dep_rate = base_rate * exp(tension*exp_factor);
else
        tension_dep_rate = base_rate / (1.0 + SIGMOID_DEPOLYM_COEFF*exp(-1.*tension*exp_factor));


#if !FOUR_POPULATIONS
        if((tension_dep_rate*frac_att + (1-frac_att)*base_rate)*dt > 1.)
                tension_dep_rate = (inv_dt - (1-frac_att)*base_rate) / frac_att;
#else
	if(tension_dep_rate*dt > 1.)
		tension_dep_rate = inv_dt;
#endif

        depolym_rate = tension_dep_rate;

//fprintf(stderr, "tension %g depolym %g\n", tension, depolym_rate);
}


void concentration::calculate_res_rate()
{
        double base_rate, exp_factor;
        double tension_dep_rate;
        base_rate = base_res_rate;

if(tension < 0.)
	exp_factor = COMPRESSION_RES_EXP_FACTOR;
else
	exp_factor = RES_EXP_FACTOR;


if(!SIGMOIDAL_RES_RATE)
        tension_dep_rate = base_rate * exp(tension*exp_factor);
else
        tension_dep_rate = base_rate / (1.0 + SIGMOID_RES_COEFF*exp(-1.*tension*exp_factor));

#if !FOUR_POPULATIONS
        if((tension_dep_rate*frac_att + (1-frac_att)*base_rate)*dt > 1.)
                tension_dep_rate = (inv_dt - (1-frac_att)*base_rate) / frac_att;
#else
        if(tension_dep_rate*dt > 1.)
                tension_dep_rate = inv_dt;
#endif

        res_rate = tension_dep_rate;

//fprintf(stderr, "tension %g res %g\n", tension, res_rate);

}


void concentration::calculate_cat_rate()
{
        double base_rate, exp_factor;
        double tension_dep_rate;
        base_rate = base_cat_rate;

if(tension < 0.)
	exp_factor = COMPRESSION_CAT_EXP_FACTOR;
else
	exp_factor = CAT_EXP_FACTOR;


if(!SIGMOIDAL_CAT_RATE)
	tension_dep_rate = base_rate * exp(tension*exp_factor);
else
	tension_dep_rate = base_rate / (1.0 + SIGMOID_CAT_COEFF*exp(-1.*exp_factor*tension));


#if !FOUR_POPULATIONS
        if((tension_dep_rate*frac_att + (1-frac_att)*base_rate)*dt > 1.)
                tension_dep_rate = (inv_dt - (1-frac_att)*base_rate) / frac_att;
#else
        if(tension_dep_rate*dt > 1.)
                tension_dep_rate = inv_dt;
#endif

        cat_rate = tension_dep_rate;

//fprintf(stderr, "tension %g cat %g\n", tension, cat_rate);

}




///only needed for two populations model
void concentration::calculate_frac_att()
{
#if EXPONENTIAL_ATTACHMENT_PROFILE
	if(tension < 0.)
//	if(tension < 4.*KT_MT_SPRING)//pretty sure the two commented lines are opposite what they should be 120307
		frac_att = 1.;
	else
//		frac_att = exp((tension - 4.*KT_MT_SPRING)*FRAC_ATT_EXP_FACTOR);
		frac_att = exp(tension*FRAC_ATT_EXP_FACTOR);
//fprintf(stderr, "%g %g\n", tension, frac_att);
	if(frac_att > 1.)
	{
		fprintf(stderr, "f too large\n");
		frac_att = 1.;
	}
	else if(frac_att < 0.)
	{
		fprintf(stderr, "f less than 0\n");
		frac_att = 0.;
	}
#elif CONST_ATTDET_RATES
/*	if(tension < 0.)
		frac_att = 1.;
	else*/ if(kt_mt_dist >= ATTACHMENT_RANGE)
		frac_att = 0.;
	else
		frac_att = 1./(1.+ base_det_rate / base_att_rate);
#else
	if(tension < 0.)//negative tension is compression
		frac_att = zero_tension_binding_frac;
	else
		frac_att = 1./(1. + exp((0.5*KT_MT_SPRING*kt_mt_dist*kt_mt_dist - EPSILON)/KbT));
#endif
}


void concentration::calculate_att_rate()
{
	double tension_dep_rate=0.;

#if (EXPONENTIAL_ATT_RATE && !LINEAR_ATT_RATE)
	tension_dep_rate = base_att_rate*exp(kt_mt_dist*att_exp_const);
#elif (LINEAR_ATT_RATE && !EXPONENTIAL_ATT_RATE)
	if(kt_mt_dist <= ATTACHMENT_RANGE)
	    tension_dep_rate = base_att_rate - base_att_rate / ATTACHMENT_RANGE * kt_mt_dist;
	else
	    tension_dep_rate = 0.;
#elif (CONST_ATTDET_RATES && !EXPONENTIAL_ATT_RATE && !LINEAR_ATT_RATE)
	if(kt_mt_dist <= ATTACHMENT_RANGE)
            tension_dep_rate = base_att_rate;
        else
            tension_dep_rate = 0.;
#else
	fprintf(stderr, "Need to choose form for attachment rate. Exiting.\n");
	exit(1);
#endif
        if(tension_dep_rate*dt > 1.)
                tension_dep_rate = inv_dt;

        att_rate = tension_dep_rate;
	//fprintf(stderr, "%g \n", att_rate);
}


void concentration::calculate_det_rate()
{
        double tension_dep_rate;

#if (CONST_ATTDET_RATES && !EXPONENTIAL_DET_RATE && !KRAMERS_DET)
#if INF_DETACH
	if(kt_mt_dist <= ATTACHMENT_RANGE)
#endif
            tension_dep_rate = base_det_rate; 
#if INF_DETACH
        else
            tension_dep_rate = inv_dt;
#endif
        if(tension_dep_rate*dt > 1.)
                tension_dep_rate = inv_dt;

#elif (KRAMERS_DET && !EXPONENTIAL_DET_RATE)
	if(tension < minimum_tension_for_det_mod)
	{
		//fprintf(stderr, "mint %g\n", minimum_tension_for_det_mod);
		det_rate = base_det_rate;
	}
	else
	{
	 if(tension < 4.*EPSILON / att_length)
	 {
	  //modification to detachment rate due to tension
	  double tension_mod = 2.*EPSILON/KbT*att_length*tension;
	  tension_mod *= 1./(2.*EPSILON+tension*att_length);
	  tension_mod *= exp(tension*att_length*(1-tension*att_length/(4.*EPSILON))/KbT);

	  if(tension_mod < 1.)
	  {
		minimum_tension_for_det_mod = tension;
		det_rate = base_det_rate;
//fprintf(stderr, "hello %g mint %g\n",  kt_mt_dist, tension);
	  }
	  else
	  {
	 	det_rate = base_det_rate*tension_mod;
		//fprintf(stderr, "hello2\n");
		if(det_rate > inv_dt)
		{
			det_rate = inv_dt;
			//fprintf(stderr, "inv_dt\n");
		}
	  }//else tension mod is >1
	 }//if(tension < 4.*EPSILON / att_length)
	 else
	   det_rate = inv_dt;
	}//else
#elif (EXPONENTIAL_DET_RATE && !KRAMERS_DET)
        tension_dep_rate = base_det_rate*exp(kt_mt_dist*det_exp_const);
	if(tension_dep_rate*dt > 1.)
		tension_dep_rate = inv_dt;
	det_rate = tension_dep_rate;
#else
	fprintf(stderr, "need to choose detachment mode. Exiting\n");
#endif

}




//Note: this scheme for keeping track of "current" monos left won't work....
// trying a scheme where in kinetics() temp = conc[i].polym(temp) where polym returns a double for amt polymed by filaments of a certain length
//should pass conc by const ref -- better performance
double concentration::polymerize(double prev_free_conc, double& curr_free_conc, double shorter_fil_polym, bool boundary_bool)
{
	double on_off = 1.;
	if(boundary_bool)
		on_off = 0.;
#if !FOUR_POPULATIONS	
	double polym_this_length = dt*(
                        prev_free_conc*(-1.*on_off*(frac_att*polym_rate + (1.-frac_att)*base_polym_rate)*prev_conc));
#else
	double polym_this_length = dt*(prev_free_conc*(-1.*on_off*polym_rate*prev_conc));
	//fprintf(stderr, " %12.10g", polym_this_length);
#endif

#if !CONST_FREE_MONO
	//if polym exceeds remaining free monos
	if(curr_free_conc + polym_this_length < 0.)//note polym_this_length <= 0
	{
	 if(curr_free_conc + polym_this_length / (LZ+1) > 0.)//first try a smaller amount, propto grid size
	 {
		polym_this_length *= 1./(LZ+1);
	 }
	 else//then just take the rest of the monos. error should be small
	 {
		polym_this_length = -1.*curr_free_conc;
	 }
	}
#endif

	conc += shorter_fil_polym + polym_this_length;
        curr_free_conc += polym_this_length;
//fprintf(stderr, "thispolym %12.10g prev_mono %12.10g currmono %12.10g\n", polym_this_length, prev_free_conc, curr_free_conc);

        return (-1.*polym_this_length);
}


void concentration::depolymerize(double& curr_free_conc, concentration const& longer_fil, bool boundary_bool)
{
	double on_off = 1.;
	if(boundary_bool)
		on_off = 0.;
#if !FOUR_POPULATIONS
	double depolym_this_length = dt*(
                        -1.*on_off*(frac_att*depolym_rate + (1.-frac_att)*base_depolym_rate)*prev_conc);
	double longer_fil_depolym = dt*(
			 (longer_fil.get_depolym_rate()*longer_fil.get_frac_att()
			   + (1.-longer_fil.get_frac_att())*longer_fil.get_base_depolym_rate())*longer_fil.get_prev_conc());
#else
	double depolym_this_length = dt*(-1.*on_off*depolym_rate*prev_conc);
	double longer_fil_depolym = dt*(longer_fil.get_depolym_rate()*longer_fil.get_prev_conc());

	//fprintf(stderr, " %12.10g %12.10g", depolym_this_length, longer_fil_depolym);

#endif
        conc += depolym_this_length + longer_fil_depolym;

//don't think I need these lines:
#if CONST_FREE_MONO // think it should be a NOT
	curr_free_conc -= depolym_this_length;//depolym_this_length >= 0 // should be a += sign if it's >=0
#endif
//fprintf(stderr, "thisdepolym %12.10g currmono %12.10g\n", depolym_this_length,  curr_free_conc);
}


void concentration::rescue_switch(concentration const& cat_fil)
{
#if !FOUR_POPULATIONS
	double change = dt*(
			(cat_fil.get_frac_att()*cat_fil.get_res_rate() 
			  + (1.-cat_fil.get_frac_att())*cat_fil.get_base_res_rate())*cat_fil.get_prev_conc()
			- (frac_att*cat_rate + (1.-frac_att)*base_cat_rate)*prev_conc
			);
#else
	double change = dt*(cat_fil.get_res_rate()*cat_fil.get_prev_conc() - cat_rate*prev_conc);
//	fprintf(stderr, " %12.10g", change);
#endif
        conc += change;
//fprintf(stderr, "resswitch %g\n", change);
}


void concentration::catastrophe_switch(concentration const& res_fil)
{
#if !FOUR_POPULATIONS
	double change = dt*(
			(res_fil.get_frac_att()*res_fil.get_cat_rate() 
			  + (1.-res_fil.get_frac_att())*res_fil.get_base_cat_rate())*res_fil.get_prev_conc()
			- (frac_att*res_rate + (1.-frac_att)*base_res_rate)*prev_conc
			);
#else
	double change = dt*(res_fil.get_cat_rate()*res_fil.get_prev_conc() - res_rate*prev_conc);
//	fprintf(stderr, " %12.10g\n", change);
#endif
        conc += change;
//fprintf(stderr, "catswitch %g\n", change);
}


void concentration::attach_onoff(concentration const& det_fil)
{
	double change = dt*(det_fil.get_att_rate()*det_fil.get_prev_conc() - det_rate*prev_conc);
	
	conc += change;
}


void concentration::detach_onoff(concentration const& att_fil)
{
	double change = dt*(att_fil.get_det_rate()*att_fil.get_prev_conc() - att_rate*prev_conc);
	
	conc += change;
}

