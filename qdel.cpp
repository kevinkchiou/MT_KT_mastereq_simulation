#include <iostream>
#include <stdlib.h>

#define MIN_RUN 201
#define MAX_RUN 314

int main()
{
char command1[64];
char command2[64];
int ii;

FILE *infile;
char inname[72];
sprintf(inname, "qdellist");
infile = fopen(inname, "r");

while(fscanf(infile, "%i", &ii) != EOF)
{
sprintf(command2, "qdel %i",ii);
system(command2);
}

fflush(infile);
fclose(infile);


return 0;
} 

