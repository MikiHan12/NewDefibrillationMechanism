#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "myutils_v2.h"

void  setup_barkley_model(double *params, int Nx) {
	// Set the params array to the default parameters.
	// Supply an array of length 23 on input.
	// Thresholds:
        double a = 0.6;//0.6
        double b = 0.075;//0.075
        double epsilon = 0.02;
	double Dt = 0.0016; // timestep
	params[1]=a; params[2]=b; params[3]=epsilon; params[4]=Dt; 
}

void barkley_model (double *u, double *v, double *J_other, int N, double *params) {
	
	/* 
	Advances the N instances of the dynamical variables Barkley’s atrial 2V-SIM model one timestep.
	 On input, u and v are arrays of length N containing the old values of these dynamical variables;
	 On output, they contain the new values at the next timestep.
	 On input, the array J_other currently like the gap junction currents and externally injected currents. For this
	 current, current flowing into the intracellular space of the cell is consider positive.  Sign and dimensions
	 of the current are consistent with "du/dt = J_other + internal currents".
	 Default parameters may be defined using the routine above
	 */
	
	double a=params[1], b=params[2], epsilon=params[3], Dt=params[4];
	double old_v;	
	int i;
	for (i=0;i<N;i++) {
		
		old_v=v[i];
                v[i] += Dt*(u[i]-v[i]);
		u[i] += Dt*(u[i]*(1-u[i])*(u[i]-(old_v+b)/a)/epsilon+J_other[i]);
	}
	
}
