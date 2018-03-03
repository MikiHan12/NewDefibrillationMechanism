#include <stdlib.h>
#include <stdio.h>
// v2: 
double *initialize (const char *argname, int length, double value){
	int i;
	double *arg = (double *)malloc(length*sizeof(double));
	if (!arg) {
		printf ("Error allocating memory for array %s\n",argname);
		exit(1);
	}
	for (i=0;i<length;i++) {
		arg[i]=value;
	}
	return arg;
}

int *int_initialize (const char *argname,int length, int value) {
	int i;
	int *arg = (int *)malloc(length*sizeof(int));
	if (!arg) {
		printf ("Error allocating memory for array %s\n",argname);
		exit(1);
	}
	for (i=0;i<length;i++) {
		arg[i]=value;
	}
	return arg;
}

