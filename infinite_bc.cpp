//infinite_bc.cpp

#include "common.h"
#include "kinetochore.h"
#include "concentration.h"


void check_kinetochore_position()
{

bool near_zero_edge = false;
bool near_max_edge = false;
if(kt_list[0].get_pos() < EDGE_PROXIMITY_LENGTH)
	near_zero_edge = true;
#if !TWO_KINETO
else if(kt_list[0].get_pos() > LZ - EDGE_PROXIMITY_LENGTH)
	near_max_edge = true;
#else//recall that kt_list[1] at 0 corresponds to that kinetochore actually being at LZ
else if(kt_list[1].get_pos() < EDGE_PROXIMITY_LENGTH)
	near_max_edge = true;
#endif

if(near_zero_edge)
	shift_simulation(POSITIVE_SHIFT);
else if(near_max_edge)
	shift_simulation(NEGATIVE_SHIFT);
}





void shift_simulation(bool positive_direction)
{
//test
//write_restart();


int direction_multiplier = -1;
if(positive_direction)
{
 direction_multiplier = 1;
}
shift_counter -= direction_multiplier;
int ii;
double temp_double;
int shift_length = (int) SHIFT_LENGTH;

vector<concentration> temp_vector_lr;
vector<concentration> temp_vector_lc;
#if FOUR_POPULATIONS
vector<concentration> temp_vector_lrd;
vector<concentration> temp_vector_lcd;
#if TWO_KINETO
vector<concentration> temp_vector_rr;
vector<concentration> temp_vector_rc;
vector<concentration> temp_vector_rrd;
vector<concentration> temp_vector_rcd;
#endif
#endif

concentration zero_growing_conc(0., true);
concentration zero_shrinking_conc(0., false);


//Shift kinetochore(s)
	temp_double = kt_list[0].get_pos();
	kt_list[0].set_pos(temp_double + direction_multiplier*SHIFT_LENGTH);
	kt_list[0].set_prev_pos();
	#if TWO_KINETO
	temp_double = kt_list[1].get_pos();
	kt_list[1].set_pos(temp_double - direction_multiplier*SHIFT_LENGTH);//note difference of sign
	kt_list[1].set_prev_pos();
	#endif


//now shift all of the filament concentrations

//first copy the vectors of concentrations
for(ii = 0; ii < left_res.size(); ii++)
{
	temp_vector_lr.push_back(left_res[ii]);
        temp_vector_lc.push_back(left_cat[ii]);
#if FOUR_POPULATIONS
        temp_vector_lrd.push_back(left_res_det[ii]);
        temp_vector_lcd.push_back(left_cat_det[ii]);
#if TWO_KINETO
        temp_vector_rr.push_back(right_res[ii]);
        temp_vector_rc.push_back(right_cat[ii]);
        temp_vector_rrd.push_back(right_res_det[ii]);
        temp_vector_rcd.push_back(right_cat_det[ii]);
#endif
#endif
} 

left_res.clear();
left_cat.clear();
#if FOUR_POPULATIONS
left_res_det.clear();
left_cat_det.clear();
#if TWO_KINETO
right_res.clear();
right_cat.clear();
right_res_det.clear();
right_cat_det.clear();
#endif
#endif


if(positive_direction)
{
for(ii = 0; ii < shift_length; ii++)
{
	left_res.push_back(zero_growing_conc);
	left_cat.push_back(zero_shrinking_conc);
	#if FOUR_POPULATIONS
	left_res_det.push_back(zero_growing_conc);
	left_cat_det.push_back(zero_shrinking_conc);
	#if TWO_KINETO
	right_res.push_back(temp_vector_rr[shift_length + ii]);
	right_cat.push_back(temp_vector_rc[shift_length + ii]);
	right_res_det.push_back(temp_vector_rrd[shift_length + ii]);
	right_cat_det.push_back(temp_vector_rcd[shift_length + ii]);
	#endif
	#endif
}

for(ii=shift_length; ii < LZ+1; ii++) 
{
	left_res.push_back(temp_vector_lr[ii-shift_length]);
	left_cat.push_back(temp_vector_lc[ii-shift_length]);
	#if FOUR_POPULATIONS
	left_res_det.push_back(temp_vector_lrd[ii-shift_length]);
	left_cat_det.push_back(temp_vector_lcd[ii-shift_length]);
	#endif
}
#if TWO_KINETO
for(ii=2*shift_length; ii < LZ+1; ii++)
{
	right_res.push_back(temp_vector_rr[ii]);
	right_cat.push_back(temp_vector_rc[ii]);
	right_res_det.push_back(temp_vector_rrd[ii]);
	right_cat_det.push_back(temp_vector_rcd[ii]);
}
for(ii = 0; ii < shift_length; ii++)
{
        right_res.push_back(zero_growing_conc);
        right_cat.push_back(zero_shrinking_conc);
        right_res_det.push_back(zero_growing_conc);
        right_cat_det.push_back(zero_shrinking_conc);
}
#endif
//make sure to set correct tension for new sets of attached filaments  --- actually this is not necessary because calc_rates() runs before
//kinetics and calculates/sets correct tension
//BEWARE!!!!!!!!! right filaments are indexed 0 to LZ but those numbers run from right to left in the box!
}
else//!positive_direction
{
for(ii = 0; ii < shift_length; ii++)
{
	#if FOUR_POPULATIONS
	#if TWO_KINETO
        right_res.push_back(zero_growing_conc);
        right_cat.push_back(zero_shrinking_conc);
        right_res_det.push_back(zero_growing_conc);
        right_cat_det.push_back(zero_shrinking_conc);
	#endif
        left_res_det.push_back(temp_vector_lrd[shift_length + ii]);
        left_cat_det.push_back(temp_vector_lcd[shift_length + ii]);
	#endif
        left_res.push_back(temp_vector_lr[shift_length + ii]);
        left_cat.push_back(temp_vector_lc[shift_length + ii]);
}
#if TWO_KINETO
for(ii=shift_length; ii < LZ+1; ii++)
{
        right_res.push_back(temp_vector_rr[ii-shift_length]);
        right_cat.push_back(temp_vector_rc[ii-shift_length]);
        right_res_det.push_back(temp_vector_rrd[ii-shift_length]);
        right_cat_det.push_back(temp_vector_rcd[ii-shift_length]);
}
#endif
for(ii=2*shift_length; ii < LZ+1; ii++)
{
        left_res.push_back(temp_vector_lr[ii]);
        left_cat.push_back(temp_vector_lc[ii]);
#if FOUR_POPULATIONS
        left_res_det.push_back(temp_vector_lrd[ii]);
        left_cat_det.push_back(temp_vector_lcd[ii]);
#endif
}
for(ii = 0; ii < shift_length; ii++)
{
#if FOUR_POPULATIONS
        left_res_det.push_back(zero_growing_conc);
        left_cat_det.push_back(zero_shrinking_conc);
#endif
        left_res.push_back(zero_growing_conc);
        left_cat.push_back(zero_shrinking_conc);
}
}//end of if/else positive_direction


///This is just a check before the code is ready for actual use
double ctot = 0.;
double ctot_2 = 0.;
for(ii = 0; ii < LZ+1; ii++)
{
	ctot += left_res[ii].get_conc();
	ctot += left_cat[ii].get_conc();
	#if FOUR_POPULATIONS
	ctot += left_res_det[ii].get_conc();
	ctot += left_cat_det[ii].get_conc();
	#if TWO_KINETO
	ctot_2 += right_res[ii].get_conc();
	ctot_2 += right_cat[ii].get_conc();
	ctot_2 += right_res_det[ii].get_conc();
	ctot_2 += right_cat_det[ii].get_conc();
	#endif
	#endif
}
fprintf(stderr, "Shifting. Total concentrations are %g %g\n", ctot, ctot_2);


//test
/*
int wait;
fprintf(stderr, "waiting for an int\n");
scanf("%i", &wait);


write_restart();
exit(1);
*/

}



