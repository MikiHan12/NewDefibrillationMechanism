/*
 *  velocity_calc.cpp
 *  4v_2d_magic_numbers
 *
 *  Created by Niels Otani on 9/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

#include "velocity_calc.h"

void calc_velocity (double *vx, double *vy, double M[4], double *ax, double *ay, double *wavefront_times, int npoints) {
	/* Calculate the velocity of a wavefront based on the times that the wavefront passees through npoints different
	gridpoints.
vx: on output, contains the x-component of the calculated velocity
vy: on output, contains the y-component of the calculated velocity
	M[4]: (input) the elements of the matrix calculated by setup_variance_matrix
	ax(npoints), ay(npoints): (input) contains npoints vectors pointing from the gripoint 
	where the velocity is to be calculated to the npoints gridpoints that will 
	participate in the calculation.  The first vector should point to the gridpoint
	where the velocity is calculated itself; thus should have ax(0)=0 and ay(0)=0 on input.
	wavefront_times(npoints): (input) times at which the wavefront passes through the gridpoints defined 
	by ax and ay relative to the gridpoint at which the calculation is to be performed.  Note that the
	wavefront time at this gridpoint should be in wavefront_time(0) on input.
npoints: (input) the number of points that will participate in the calculation.
	This routine works by minimizing the sum of the squares of (time_difference-a*grad(t_wavefront))
	*/
	double A = M[0], B = M[1], C = M[2], D = M[3];
	double b1 = 0.0, b2 = 0.0;
	int ipoints;
	for (ipoints=0;ipoints<npoints;ipoints++) {
		b1 += (wavefront_times[ipoints]-wavefront_times[0])*ax[ipoints];
		b2 += (wavefront_times[ipoints]-wavefront_times[0])*ay[ipoints];
	}
	double det = A*D-B*C;
	double dtdx = (D*b1-B*b2)/det;
	double dtdy = (-C*b1 + A*b2)/det;
	double gradtsq = dtdx*dtdx + dtdy*dtdy;
	*vx = dtdx/gradtsq;
	*vy = dtdy/gradtsq;
	return;
}

void setup_variance_matrix (double M[4], double *ax, double *ay, int npoints) {
	/* Setup the matrix used to calculate the velocity vector.
		M[4]: (output) contains the elements of the 2x2 matrix.  Supply a matrix of length 4 on input.
		ax(npoints), ay(npoints): (input) contains npoints vectors pointing from the gripoint 
			where the velocity is to be calculated to the npoints gridpoints that will 
			participate in the calculation.  The first vector should point to the gridpoint
			where the velocity is calculated itself; thus should have ax(0)=0 and ay(0)=0 on input.
		npoints: (input) the number of points that will participate in the calculation
		This routine need only be called once before calculations are performed, since it is a
		constant matrix.
	*/
	int k, ipoints;
	for (k=0;k<4;k++) {
		M[k] = 0;
	}
	for (ipoints=0;ipoints<npoints;ipoints++) {
		M[0] += ax[ipoints]*ax[ipoints];
		M[1] += ax[ipoints]*ay[ipoints];
		M[3] += ay[ipoints]*ay[ipoints];
	}
	M[2] = M[1];
	return;
}

int ready_to_calc_velocity (double current_time, double *wavefront_times, int npoints, double time_diff_threshold) {
	/* Check to see whether all the wavefront times are in place to do the velocity calculation.
		The calculation will be done when the last wavefront time involved in the calculation has just been updated
		(so that one of the wavefront_times is the current_time) and all the wavefront_times are less than
		time_diff_threshold behind the current_time.
		current_time: (input) the current time in the simulation
		wavefront_times(npoints): (input) a vector of length npoints containing the last time a wavefront passed
			each gripoint participating in the calculation.
		npoints: (input) the number of gridpoints participating in the calculation.
		time_diff_threshold: (input) all gridpoint wavefront_times must be within this value of the current_time;
			otherwise, the wavefront times will be assumed to not be from the same wavefront, and the calculation
			will not be ok'ed by this routine.
		return value: (output) = 1 if it is time to do the velocity calculation; 0 otherwise.
	*/
	int ipoints;
	// Check to see one of the points has just been updated.
	// (This prevents the same calculation from being performed over ana over)
	for (ipoints=0;ipoints<npoints;ipoints++) {
		if (wavefront_times[ipoints]==current_time)
			break;
	}
	if (ipoints==npoints)
		return (0);
	// Check to see that all the points have recently been updated:
	for (ipoints=0;ipoints<npoints;ipoints++) {
		if ( !(current_time-wavefront_times[ipoints] < time_diff_threshold ) )
			return(2);
	}
	return (1);
}
		


		
	
	

		
		
		
	


