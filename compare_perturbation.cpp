#include <iostream>
#include <math.h>
#include <vector>
#include <stdio.h>

using namespace std;


#define TRIAL 1331
#define MAX_TIME 10000
#define INTERVAL 20
#define THRESHOLD 82.e-5


int main()
{

int ii, jj, kk;
char ss_conc_name[128];
FILE *ss_conc_file;
char inname[128];
FILE *infile;
double tempx, tempgc, tempsc;
double temp1, temp2, temp3, temp4, temp5;

vector<double> ss_displacement;
vector<double> ss_grow_conc;
vector<double> ss_shrink_conc;
vector<double> ss_kt_pos;

vector<double> pert_pos;
vector<double> pert_grow_conc;
vector<double> pert_shrink_conc;
double current_kt_pos, displacement, ssconc;

double dist_from_ss;

bool end_of_disp_list;

FILE *outfile;
char outname[128];


sprintf(ss_conc_name, "concs%6.6iu", TRIAL);
ss_conc_file = fopen(ss_conc_name, "r");
while(fscanf(ss_conc_file, "%lg %lg %lg", &tempx, &tempgc, &tempsc) != EOF)
{
	ss_displacement.push_back(tempx);
	ss_grow_conc.push_back(82.*tempgc);
	ss_shrink_conc.push_back(82.*tempsc);
}
fflush(ss_conc_file);
fclose(ss_conc_file);


sprintf(outname, "ssdeviation%6.6i", TRIAL);
outfile = fopen(outname, "w");

for(ii = 0; ii <= MAX_TIME; ii+=INTERVAL)
{
sprintf(inname, "restart/%iunpert/restart%6.6i_%8.8i", TRIAL, TRIAL, ii);
infile = fopen(inname, "r");
for(jj = 0; jj <= 2000; jj++)
	fscanf(infile, "%lg %lg %lg %lg %lg", &temp1, &temp2, &temp3, &temp4, &temp5);
fscanf(infile, "%lg", &temp1);
fscanf(infile, "%lg %lg", &temp1, &tempx);
ss_kt_pos.push_back(tempx);
fflush(infile);
fclose(infile);
}


//SS measured.
//////////////////////////////////////////////////////

for(ii = 0; ii <= MAX_TIME; ii+= INTERVAL)
{
sprintf(inname, "restart/%ipert_large/restart%6.6i_%8.8i", TRIAL, TRIAL, ii);
infile = fopen(inname, "r");
for(jj = 0; jj <= 2000; jj++)
{
        fscanf(infile, "%lg %lg %lg %lg %lg", &temp1, &temp2, &temp3, &temp4, &temp5);
	if((temp2 > THRESHOLD)||(temp3>THRESHOLD))
	{
		pert_pos.push_back(temp1);
		pert_grow_conc.push_back(temp2);
		pert_shrink_conc.push_back(temp3);
	}
}

fscanf(infile, "%lg", &temp1);
fscanf(infile, "%lg %lg", &temp1, &tempx);

current_kt_pos = tempx;

jj = 0;
end_of_disp_list = false;
dist_from_ss = 0.;
for(kk = 0; kk < pert_pos.size(); kk++)
{
displacement = pert_pos[kk] - current_kt_pos;	
//fprintf(stderr, "disp %g\n" ,  displacement);
while(displacement > ss_displacement[jj])
{
//fprintf(stderr, "while disp %g ssdisp %g\n" ,  displacement, ss_displacement[jj]);

	jj++;
	if(jj == ss_displacement.size())
	{
	  end_of_disp_list = true; 
	  break;
	}
//fprintf(stderr, "2ndwhile disp %g ssdisp %g\n" ,  displacement, ss_displacement[jj]);
}

if(!end_of_disp_list)
{
if(jj > 0)
  ssconc =  ss_grow_conc[jj]- (ss_displacement[jj] - displacement)/(ss_displacement[jj] - ss_displacement[jj-1])*(ss_grow_conc[jj] - ss_grow_conc[jj-1]);
else//extrapolating here
  ssconc = ss_grow_conc[jj]- (ss_displacement[jj] - displacement)/(ss_displacement[jj+1]-ss_displacement[jj])*(ss_grow_conc[jj+1]-ss_grow_conc[jj]);
/*
if(jj>0)
fprintf(stderr, "ssconc %g from %g and %g, pert: %g, pertpos %g\n", ssconc, ss_grow_conc[jj], ss_grow_conc[jj-1], pert_grow_conc[kk], pert_pos[kk]);

int wait;
scanf("%i", &wait);
*/

}
else
{
ssconc = ss_grow_conc[jj-1] + (ss_grow_conc[jj-1] - ss_grow_conc[jj-2])*(displacement - ss_displacement[jj-1]);
        //Need to extrapolate
}

if(ssconc < 0.)
  ssconc = 0.;


dist_from_ss += (pert_grow_conc[kk] - ssconc)*(pert_grow_conc[kk] - ssconc);


//repeat for shrinking
if(!end_of_disp_list)
{
if(jj > 0)
  ssconc =  ss_shrink_conc[jj]- (ss_displacement[jj] - displacement)/(ss_displacement[jj] - ss_displacement[jj-1])*(ss_shrink_conc[jj] - ss_shrink_conc[jj-1]);
else
  ssconc = ss_shrink_conc[jj]- (ss_displacement[jj] - displacement)/(ss_displacement[jj+1]-ss_displacement[jj])*(ss_grow_conc[jj+1]-ss_grow_conc[jj]);

}
else
{
ssconc = ss_shrink_conc[jj-1] + (ss_shrink_conc[jj-1] - ss_shrink_conc[jj-2])*(displacement - ss_displacement[jj-1]);
        //Need to extrapolate
}


if(ssconc < 0.)
  ssconc = 0.;


dist_from_ss += (pert_shrink_conc[kk] - ssconc)*(pert_shrink_conc[kk] - ssconc);


}//for(kk...


dist_from_ss += (current_kt_pos - ss_kt_pos[ii/INTERVAL])*(current_kt_pos - ss_kt_pos[ii/INTERVAL]); 
dist_from_ss = sqrt(dist_from_ss);

fprintf(outfile, "%i %g\n", ii, dist_from_ss);

pert_grow_conc.clear();
pert_shrink_conc.clear();
pert_pos.clear();

}//for loop over times 

fflush(outfile);
fclose(outfile);



return 0;

}//end main

