#include<stdio.h>
#include<cstdlib>

#define FIRST 408
#define LAST 418

int main()
{
	char command[96];
	int ii;	


	for(ii = FIRST; ii <= LAST; ii++)
	{
	   sprintf(command, "cp restart/restart%6.6i restart/restart%6.6i_01", ii, ii);
	   system(command);
	}

return 0;

}

