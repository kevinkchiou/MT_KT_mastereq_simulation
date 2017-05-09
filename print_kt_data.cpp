#include "common.h"
#include "kinetochore.h"

void print_kt_data()
{

	print_kt_pos();

}


void print_kt_pos()
{
	int ii;

	
	double shift =  shift_counter*SHIFT_LENGTH;
	
	fprintf(ktposfile, "%g", sim_step*dt);
	fprintf(ktposfile, " %12.10g", kt_list[0].get_prev_pos() + shift);
	#if TWO_KINETO
	fprintf(ktposfile, " %12.10g", (double)LZ - kt_list[1].get_prev_pos() + shift);
	#endif
	fprintf(ktposfile, "\n");

	fflush(ktposfile);
}
