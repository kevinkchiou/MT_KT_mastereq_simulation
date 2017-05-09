#include "common.h"
#include "kinetochore.h"

void dash_pot(double left_mts, double right_mts, double kkforce)
{

double force_1 = kkforce*DRAG_OVER_EFFECTIVE - right_mts*RV_DRAG_OVER_EFFECTIVE + left_mts*SUM_DRAG_OVER_EFFECTIVE;
double force_2 = kkforce*DRAG_OVER_EFFECTIVE - left_mts*RV_DRAG_OVER_EFFECTIVE + right_mts*SUM_DRAG_OVER_EFFECTIVE;

/*
fprintf(stderr, "force_1 = %g, kkforce = %g rightdp %g leftnorm %g\n", force_1, kkforce*DRAG_OVER_EFFECTIVE, -1.*right_mts*RV_DRAG_OVER_EFFECTIVE, left_mts*SUM_DRAG_OVER_EFFECTIVE);
fprintf(stderr, "force_2 = %g, kkforce = %g leftdp %g rightnorm %g\n", force_2, kkforce*DRAG_OVER_EFFECTIVE, -1.*left_mts*RV_DRAG_OVER_EFFECTIVE, right_mts*SUM_DRAG_OVER_EFFECTIVE);

int wait;
scanf("%i", &wait);
*/

kt_list[0].move(force_1*dt);
kt_list[1].move(force_2*dt);

}


