#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <vector>

using namespace std;

#define FIRST_TRIAL first
#define LAST_TRIAL last
#define NUMSKIP 400
#define NUM_RESTARTS 50
#define LBOX 2000
#define MAX_DIST 30.
#define MIN_CONC 5.e-5


int main(int argc, char **argv)
{
int option;
int first, last;
while( (option = getopt(argc, argv, "f:l:")) != -1)
{
 switch(option)
 {
   case 'f':
      first = atoi(optarg); break;
   case 'l':
      last = atoi(optarg); break;
   case '?':
    fprintf(stderr, "unknown command line option. Exit\n"); return 1;
 }
}


FILE *infile;
FILE *outfile;
char inname[96];
char outname[96];

int ii,jj,kk;
int x;
int temp_int;
double pr, pc, kt;
double ptot;
double tempr, tempc;
vector<double> pr_list;
vector<double> pc_list;



for(ii = FIRST_TRIAL; ii <= LAST_TRIAL; ii+=1)
{
sprintf(outname, "concs%6.6i", ii);
outfile = fopen(outname, "w");

for(jj = 0; jj < NUM_RESTARTS; jj++)
{
//sprintf(inname, "restart/1331pert_large/restart%6.6i_%8.8i", ii, jj*NUMSKIP);
sprintf(inname, "restart/restart%6.6i_%8.8i", ii, jj*NUMSKIP);
infile = fopen(inname, "r");
ptot = 0.;

for(kk = 0; kk <= LBOX; kk++)
{
fscanf(infile, "%i %lg %lg %lg %lg", &x, &pr, &pc, &tempr, &tempc);

pr_list.push_back(pr);
pc_list.push_back(pc);
ptot += pr + pc;
}//for(kk

fscanf(infile, "%i %i %lg", &temp_int, &temp_int, &kt);

for(kk = 0; kk <= LBOX; kk++)
{
//if(fabs(kt-kk) < MAX_DIST)
if((pr_list[kk]/ptot > MIN_CONC) || (pc_list[kk]/ptot > MIN_CONC))
  fprintf(outfile, "%g %g %g\n", kk-kt, pr_list[kk]/ptot, pc_list[kk]/ptot);
}
pr_list.clear();
pc_list.clear();
fflush(infile);
fclose(infile);
}//for(jj

fflush(outfile);
fclose(outfile);

}//for(ii

return 0;
}



