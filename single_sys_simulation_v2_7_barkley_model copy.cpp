// Modified from flavios_tailored_model/alisa_modified/4_varonff_nfo_mod.cpp.
// Modification started Sept 5, 2006  --NFO

// Dec 7, 2006:  Modified from "first attempt" to include the atrial_4v_sim_mod2.cpp
// file contained in build_ion_channel_model_steep_exp_RA_v2.  
// This code reads in the file, coeffs.bin, which contains coefficients for the w-gate advance.
// It is created by improve_LA_RA_transition_v2.m.  Don't forget to set Nx equal to its value
// before running that code.

// June 16, 2010: Modified from atrial_4v_sim_steep.cpp to do a single system instead of two.
// This code is intended to be used to test my magic numbers defibrillation idea.  As such, I 
// need a spiral wave code that has a steep APD restitution function that creates slowly traveling
// wavebacks.  Hopefully this code will do this.

// April 12, 2011: single_system_simulation_v2: Add APD vs. next DI and cycle length diagnostics.

// April 19, 2011: v2_0_field_stim: Extend code to allow field stimuli.

// April 20, 2011: v2_2_field_stim: Add run continuation capability. 

// April 26, 2011: v2_3_field_stim: Allow batching with respect to stimulus times.

// July 19, 2011: v2_5_field_stim: Allow batching with respect to the first (and only) stimulus time.

// July 19, 2011: v2_6_field_stim: Add a central obstacle.

// July 20, 2011: v2_6_1_field_stim: Add the second stimulus again.

// Nov 16, 2015: v2_6_1s_field_stim: Simplify and comment..

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include "velocity_calc.h"
#include "barkley_model.h"
#include "myutils_v2.h"
#include "single_sys_simulation_v2_7_barkley.h"

// Main program:
int main (int argc, char * const argv[]) { 
	
	const int N_stim = 2;  // Number of electric field stimuli to apply.
	double Field_stim_time[N_stim]; // when to apply them.
	
	int i_stim;
	
	// Allow the possibility of doing several runs, each with the same first electric field
	// stimulus time, but each with a different second stimulus time:
	Field_stim_time[0] = 1580;
	for (Field_stim_time[1] = 1790; Field_stim_time[1]<=1790; Field_stim_time[1] += 10.0) {
		printf("Running field stim at %7.1f ms\n",Field_stim_time[1]);
		simulation(Field_stim_time,N_stim);
	}
	
	return 0;
	
}

// The simulation itself:

int simulation (double *Field_stim_time, int N_stim) { 
	/* This is a two-dimensional simulation of our attempts to apply electric field stimuli to stop a spiral wave that is rotating around a circular obstacle.  The full simulation consists of two stages: In the first stage, two current stimuli are applied to start a spiral wave rotating around the circular obstacle.  In the second stage, we use electric field stimuli to try to stop the spiral wave.
		
	   This code allows you to save the simulation after the spiral wave has formed, and then reload it again and again into the code, so you can continue from that point, so you don't have to run the spiral wave formation stage over and over. To start the simulation from the beginning, set the variable continuation = 0. Following the formation of the spiral wave, the simulation will save a copy of the state of the simulation to a file.  You can then start a new simulation at that point by setting continuation = 1, and defining the variable continuation_filename with the name of the continuation file.
	 
		This code send all its results to files.  The files that are of the form uxxxxxx, where xxxxxx is a number, contains the membrane potential information for the timestep xxxxxx. By default, this is the only information that is outputted.  By setting other flags, other information about the APDs, DIs and velocities are available.
	 */
	
	// Simulation parameters:
	const int Nx = 144, Ny = 144; 	//450,450  System computational grid size
	const double t_end = Field_stim_time[1] + 600.0;  // Time of end of the simulation
	const double Dx = 0.1667, Dy = 0.1667, Dz = 1.0;  //  Grid spacing in x and y (in cm?)
	const double D = 1; // Diffusion coefficienthttps://mycourses.rit.edu/d2l/lp/documents/aCo2IEDcUcEil1If0UHT3259205/1/6
	const int ntplot = 50, ntplot_start = 0; // plotting interval (in timesteps), and timestep to turn on plotting
	//const int nthist = 40;  // recording interval for the history files. nthist=40 corresponds to 2ms when Dt = 0.05.
	const int Napd = 240, Ndi = 300, Nvabs = 100; // Dimension of additional diagnostic plots, when used, in pixels
	const double obs_radius	= 0.0; // Obstacle radius (cm)
    const double x0 = 0, y0 = -7.0; //1.5,-4.5// //the center of the second obstacle
    const double obs2_radius = 0.99;//1.5//10//Second obstacle radius (cm)
    const double I = 75;//75
    const int I_E_duration = 184; // default = 23
    const int stimulus_it = 9350;//good:9350// when to apply the obstacle stimulus
	
	double params[5]; // Array that will hold the parameters for the simulation.
	setup_barkley_model(params,Nx); // Call the routine that defines the array "params"
	double Dt = params[4]; // get the timestep size from the "params" array.
	int Nt = 20000,Nt_NoObs=0;
    
    //lround(t_end/Dt) + 1; // Number of timesteps to run
	
	double stim_duration = Dt*10; // ms
	double stim_amplitude = 2.5*(-10); //-15// (Minimum activation current is -0.11)
	
	/* If the simulation is started from the beginning, two local stimuli are required to start 
	 the spiral wave. These are not the same as the electric field stimuli, used later to try to
	 stop the spiral wave. */
	// const double CL2 = Dt*2650;// good 
	const double CL2 = Dt*3150;//3150lower spiral wave
	
	const double DI_rel_ref = 20.0; // Time (ms) that must have elapsed since last depolarization for field stimulus to create depolarization
	
	int  i, j,k;
	FILE *fp;
	char filename[256];
	
	// Pointers to arrays to hold the dynamical variables:
	double *u, *v, *u_old;
	
	// Pointers to arrays holding additional diagnostic information.
	// Comment out if not needed:
	/*
	double *last_apd, *last_di, *cycle_length, *last_depolarization, *last_repolarization;
	int *apd_di, *apd_nextdi, *vabs_di;
	double *vabs, *vbackabs;
	 */
	
	int it;
	int i_stim;
	
	// Run name:
	char runname[256];
	sprintf(runname,"%ld",lround(Field_stim_time[0]));
	for (i_stim=1;i_stim<N_stim;i_stim++) {
		sprintf(runname,"%s_%li",runname,lround(Field_stim_time[i_stim]));
	}
	
	// Continuation controls (see description above):
	const int continuation = 0;
	const char *continuation_filename = "continuation_file_1580_1790_4699";
	
	// Data folder:
	sprintf(filename,"data_%s",runname);
	mkdir(filename,0000755);
		
	// Initial conditions for the dynamical variables:
	// NOTE: these arrays are one-dimensional, even though they represent dynamical variable
	// data on a 2D grid that is Nx by Ny. u(i,j) = u[i+j*Nx]; that is u(i,j) is stored 
	// in the (i+j*Nx)th element of the array.  
	// Note that i runs from 0 to Nx-1 and j runs from 0 to Ny-1.
	u = initialize ("u", Nx*Ny, 0);
	v = initialize ("v", Nx*Ny, 0);
		
	u_old = initialize ("u_old", Nx*Ny, 0);		
	
	
	
	it = 0;
	
	if (continuation==1) { // If 1, overwrite initial conditions from a continuation file:
		
		fp = fopen(continuation_filename,"rb");
		if (!fp) {printf("Couldn't read continuation file %s\n",continuation_filename); exit(1); }
		printf("Reading continuation file\n");
		
		fread(u,sizeof(double),Nx*Ny,fp);
		fread(v,sizeof(double),Nx*Ny,fp);
		
		
		fread(u_old,sizeof(double),Nx*Ny,fp);
				
		fread(&it,sizeof(int),1,fp);
		fclose(fp);
		
	}
	else { // If continuation is not 1, start from the beginning.
		// Set last timestep equal to the timestep marking the end of the spiral wave setup period:
		Nt = 20000;
        Nt_NoObs=0;// Comment out if a full run (spiral wave setup; then electric field pulses) is desired.
	}
	
	// Define arrays that will contain..
	// Currents flowing through the gap junctions..
	double *x_current = initialize ("x_current",(Nx-1)*Ny,0.); // in the x-direction,
	double *y_current = initialize ("y_current", Nx*(Ny-1),0.); // and in the y-direction.
	// Current flowing out of nodes into gap junctions:
	double *I_gap = initialize ("I_gap", Nx*Ny, 0); 
	// Applied stimulus current (mainly for starting the spiral waves):
	double *J_ext = initialize ("J_ext", Nx*Ny, 0);
	// Current to give to the routine single-cell model:
	double *J_other = initialize ("J_other",Nx*Ny,0);
	// Gap junction condunctances in the x- and y-directions:
	double *G_x = initialize ("G_x",(Nx-1)*Ny,0);
	double *G_y = initialize ("G_y",Nx*(Ny-1),0);
    double *G_x_NoObs = initialize ("G_x_NoObs",(Nx-1)*Ny,0);
	double *G_y_NoObs = initialize ("G_y_NoObs",Nx*(Ny-1),0);
    
    double *I_E = initialize ("I_E",Nx*Ny,0) ;
    	
	// Define gap junction conductances.
	// The obstacle is also defined here, as a region in which the conductances are 0.
	
	// First define the conductances between cells in the x-direction:
	for (i=0;i<Nx-1;i++) {
		for (j=0;j<Ny;j++) {
			double xx1, yy1, xx2;
			xx1 = Dx*(i-0.5*Nx);
            xx2 = Dx*(i+1-0.5*Nx);
			yy1 = Dy*(j-0.5*Ny);
			
            if ( (xx1-x0)*(xx1-x0)+(yy1-y0)*(yy1-y0) >= obs2_radius*obs2_radius &&
                 (xx2-x0)*(xx2-x0)+(yy1-y0)*(yy1-y0) >= obs2_radius*obs2_radius ){
                G_x[i+j*(Nx-1)] = D/(Dx*Dx);
                G_x_NoObs[i+j*(Nx-1)] = D/(Dx*Dx);// Normal conductances outside two obstacles radius.
            }
            else{
                G_x[i+j*(Nx-1)] = 0.0; // Zero conductances inside the second obstacle.
                G_x_NoObs[i+j*(Nx-1)] = D/(Dx*Dx);
            }

		}}
	
	// Then define the conductances in the y-direction:
	for (i=0;i<Nx;i++) {
		for (j=0;j<Ny-1;j++) {
			double xx1, yy1, yy2;
			xx1 = Dx*(i-0.5*Nx);
			yy1 = Dy*(j-0.5*Ny);
			yy2 = Dy*(j+1-0.5*Ny);
			
            if ( (xx1-x0)*(xx1-x0)+(yy1-y0)*(yy1-y0) >= obs2_radius*obs2_radius &&
                 (xx1-x0)*(xx1-x0)+(yy2-y0)*(yy2-y0) >= obs2_radius*obs2_radius ) {
                G_y[i+j*Nx] = D/(Dy*Dy);
                G_y_NoObs[i+j*Nx] = D/(Dy*Dy);// Normal conductances outside two obstacles radius.
            }
            else{
                G_y[i+j*Nx] = 0.0;
                G_y_NoObs[i+j*Nx] = D/(Dy*Dy);// Zero conductances inside the second obstacle.
            }

		}}
	
	
	// MAIN TIMESTEP LOOP:
	
	while (it<Nt) {
		
		// Calculate currrents flowing through the gap junctions:
		for (i=0;i<Nx-1;i++) {
			for (j=0;j<Ny;j++) {
                if (it<Nt_NoObs){
                    x_current[i+j*(Nx-1)] = G_x_NoObs[i+j*(Nx-1)] * ( u_old[i+j*Nx] - u_old[(i+1)+j*Nx] );
                }
                else{
				    x_current[i+j*(Nx-1)] = G_x[i+j*(Nx-1)] * ( u_old[i+j*Nx] - u_old[(i+1)+j*Nx] );
                }
			}}
		for (i=0;i<Nx;i++) {
			for (j=0;j<Ny-1;j++) {
                if (it<Nt_NoObs){
                    y_current[i+j*Nx] = G_y_NoObs[i+j*Nx] * ( u_old[i+j*Nx] -u_old[i+(j+1)*Nx] );
                }
                else{
				    y_current[i+j*Nx] = G_y[i+j*Nx] * ( u_old[i+j*Nx] -u_old[i+(j+1)*Nx] );
                }
			}}
		
		// Calculate I_gap, the currents flowing out of nodes into surrounding gap junctions:
		for (i=0;i<Nx;i++) {
			for (j=0;j<Ny;j++) {
				I_gap[i+j*Nx] = 0.0;
			}}
		for (i=0;i<Nx-1;i++) {
			for (j=0;j<Ny;j++) {
				I_gap[i+j*Nx] += x_current[i+j*(Nx-1)];
			}}
		for (i=1;i<Nx;i++) {
			for (j=0;j<Ny;j++) {
				I_gap[i+j*Nx] -= x_current[(i-1)+j*(Nx-1)];
			}}
		for (i=0;i<Nx;i++) {
			for (j=0;j<Ny-1;j++) {
				I_gap[i+j*Nx] += y_current[i+j*Nx];
			}}
		for (i=0;i<Nx;i++) {
			for (j=1;j<Ny;j++) {
				I_gap[i+j*Nx] -= y_current[i+(j-1)*Nx];
			}}
		
		// Calculate external stimulus, as needed:
		if ( it==2) { // Turn on the first simulus current on timestep 2:
			for (i=0;i<Nx;i++) {
				for (j=0;j<Ny;j++) {
					if (j>Ny-20) // Only in a strip 5 cells wide on the bottom of the system
						J_ext[i+j*Nx] = stim_amplitude;
					else
						J_ext[i+j*Nx] = 0.0;
				}
			}
		}
		if ( it==2+lround(stim_duration/Dt) ) { // Turn off the first stimulus:
			for (i=0;i<Nx;i++) {
				for (j=0;j<Ny;j++)
					J_ext[i+j*Nx] = 0.0;
			}
		}
		if ( it==2+lround(CL2/Dt)) { // Turn on the second stimulus:
			for (i=0;i<Nx;i++) {
				for (j=0;j<Ny;j++) {
					if (i< Nx/3)//3*Nx/4)-1 // Only in the left half of the system
						J_ext[i+j*Nx] = stim_amplitude;
					else
						J_ext[i+j*Nx] = 0.0;
				}
			}
		}
		if (it == 2+lround(CL2/Dt)+lround(stim_duration/Dt) ) { // Turn off the second stimulus
			for (i=0;i<Nx;i++)
				for (j=0;j<Ny;j++)
					J_ext[i+j*Nx] = 0.0;
		}
        
        if (it==stimulus_it){
            for (i=0;i<Nx;i++) {
				for (j=0;j<Ny;j++) {
					
						J_ext[i+j*Nx] += I_E[i+j*Nx] ;
					
				}
			}
     	}
            
		 if (it==stimulus_it+I_E_duration){
            for (i=0;i<Nx;i++) {
				for (j=0;j<Ny;j++) {
					
						J_ext[i+j*Nx] -= I_E[i+j*Nx] ;
					
				}
			}
		}
        
        
		for (i=0;i<Nx;i++) {
			for (j=0;j<Ny;j++) {
				int index = i+j*Nx;
				J_other[index] = - J_ext[index] - I_gap[index];
			}
		}
		
        
        
        
        //Add electrical field on the boundary of the obstacles on x-direction.
        //Decide the points that will be added Current I onto. 
        //int k_Pos=0;
        //int k_Neg=0;
        
        for (j=0;j<Ny;j++) {
		for (i=0;i<=Nx-1;i++) {
            double  xx, xx1, yy, yy1;
            
			xx = Dx*(i-0.5*Nx);
			yy = Dy*(j-0.5*Ny);
            
            /*// For NW-SE E-field:
            xx1 = Dx*((i+1)-0.5*Nx);
            yy1 = Dx*((j-1)-0.5*Ny);*/ 
            
            /*// For SW-NE-field:
            xx1 = Dx*((i+1)-0.5*Nx);
            yy1 = Dx*((j+1)-0.5*Ny); */
          
            
           // For horizontal E-field
            xx1 = Dx*((i+1)-0.5*Nx);
            yy1 = Dx*((j)-0.5*Ny);  
            
            //Add the point outside of the obstcle into the array if its next point is inside.
			if ((xx-x0)*(xx-x0)+(yy-y0)*(yy-y0) >= obs2_radius*obs2_radius) 
                if ( (xx1-x0)*(xx1-x0)+(yy1-y0)*(yy1-y0) < obs2_radius*obs2_radius ){
                   
                    I_E[i+j*Nx] = I;//*(sqrt(2)/2);//(i,j))
                    
                }
               
            //Add the point outside of the obstcle into the array if its previous point is inside.
            if ((xx-x0)*(xx-x0)+(yy-y0)*(yy-y0) < obs2_radius*obs2_radius )
                if ((xx1-x0)*(xx1-x0)+(yy1-y0)*(yy1-y0) >= obs2_radius*obs2_radius ){
                    //I_plus_x[k_Pos]=i+1;
                    //I_plus_y[k_Pos]=j;
                
                    // For NW-SE E-field:
                   // I_E[(i+1)+(j-1)*Nx] = -I;//*(sqrt(2)/2);
                    
                    // For SW-NE E-field:
                    //I_E[(i+1)+(j+1)*Nx] = -I;
                
                    // For horizontal E-field:
                    I_E[(i+1)+(j)*Nx] = -I;

                }
        }}
         
        
        
		// Update the main dynamical variables:
		// Notes:
		// (1) If you intend to use a different ion channel model, it is best to 
		// write a new function, basing it on the existing function; then call it
		// by its new name here. (If you use a different header file to declare it,
		// change that too, above, where the .h files are declared.) This allows you 
		// to build your own set of models that you can drop in here as needed, without
		// having to modify anything else.
		// (2) Since these are single_cell models, stimulus currents and gap junction
		// currents look like current from the outside.  These currenta are therefore
		// passed to the function through J_other, the array that is assumed to contain
		// currents introduced to the cell from the outside.
		barkley_model (u, v, J_other, Nx*Ny, params);  
		
		
		// Inform the user of the timestep every so often:
		if (it/ntplot*ntplot==it)
			printf("Timestep = %i, time = %f\n",it,it*Dt);
		
		// Write arrays into binary files every ntplot timesteps:
		if ( (it/ntplot*ntplot==it) && (it>=ntplot_start) ) {
			// printf("Outputting\n");
			sprintf(filename,"data_%s/u%07i",runname,it); fp = fopen(filename,"wb"); 
			if (!fp) {printf("Couldn't find data directory\n"); exit(1);}
			fwrite(u,sizeof(double),Nx*Ny,fp); fclose(fp);
            sprintf(filename,"data_%s/v%07i",runname,it); fp = fopen(filename,"wb"); 
			if (!fp) {printf("Couldn't find data directory\n"); exit(1);}
			fwrite(v,sizeof(double),Nx*Ny,fp); fclose(fp);
		}
		
		old_new(u_old,u,Nx,Ny);
		
		it++;
		
		// Write a continuation file at the end of the simulation.
		// This file may be used to restart the simulation at the point where this simulation ends.
		
		if (it==Nt-1) {
			
			printf("Writing continuation file\n");
			sprintf(filename,"continuation_file_%s_%i",runname,it);
			fp = fopen(filename,"wb");

			fwrite(u,sizeof(double),Nx*Ny,fp);
			fwrite(v,sizeof(double),Nx*Ny,fp);
			
			fwrite(u_old,sizeof(double),Nx*Ny,fp);
			
			
			fwrite(&it,sizeof(int),1,fp);
			
			fclose(fp);
		}
	}
	
	/*
	 for (i=0;i<24;i++) {
	 fclose(fphist[i]);
	 }
	 */
	
	return 0;
	
}

void old_new (double old[], double update[], int Nx, int Ny){
	int i, j;
	for (i=0;i<Nx;i++) {
		for (j=0;j<Ny;j++) {
			old[i+j*Nx]=update[i+j*Nx];
        }
	}
}



	
