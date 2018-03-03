/*
 *  velocity_calc.h
 *  4v_2d_magic_numbers
 *
 *  Created by Niels Otani on 9/7/06.
 *  Copyright 2006 __MyCompanyName__. All rights reserved.
 *
 */

void calc_velocity (double *vx, double *vy, double M[4], double *ax, double *ay, double *wavefront_times, int npoints);
void setup_variance_matrix (double M[4], double *ax, double *ay, int npoints); 
int ready_to_calc_velocity (double current_time, double *wavefront_times, int npoints, double time_diff_threshold);

