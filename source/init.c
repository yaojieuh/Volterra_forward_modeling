#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include <time.h>

#include "init.h"
void init_acou_homog(FILE *filefprintf, double dx, double dz, int nx, int nz,double c0, double *vel, double *vpe){
 
 int i;
 long index;
 int ndims = nz*nx;

 fprintf(filefprintf, "--!\tAcoustic media initialization homgeneous ...\n");
 fprintf(filefprintf, "--!\t                                \n");

  for(i=0; i<ndims; i++) {
	vel[i] = 1500.0;
	vpe[i]=1-c0*c0/vel[i]/vel[i];
  }

}

void init_acou_layer(FILE *filefprintf, double dx, double dz, int nx, int nz, double c0, double *vel, double *vpe){
 
 int i, j;
 long index;
 double x,z;
 
 int nzl1= 20, nzl2=40;



 fprintf(filefprintf, "--!\tAcoustic media initialization three layers ...\n");
 fprintf(filefprintf, "--!\t                                \n");
 FILE* file2;  
 char fname2[100];

  	sprintf(fname2,"Velocitymodel.dat");
  	file2= fopen( fname2, "w" );
	  for(j=0; j<nx; j++){
		x=j*dx;
	      	for(i=0; i<nz; i++) {
			z=(i)*dz;
			index= j*nz+i;
			if (i<nzl1){
				vel[index]=1500;
				
			}
		 	else if (i<nzl2){
                  		vel[index]=2000;
                 
	      		}
	     		else{
                  		vel[index]=1500;
                  	 }
			vpe[index]=1-c0*c0/vel[index]/vel[index];

	 		fprintf(file2, "%lf %lf %lf  %lf\n", x, z,vel[index],vpe[index]);
	 	 }
		fprintf(file2,"  \n");
	}
	fclose(file2);
}





void init_source_ricker_fwps(FILE *filefprintf, int nt, double dt,int nw, double dw, double *source, double *fren, double *sw){

  int i=0; 
  double at=0.0, arg=0.0;
  double sig=1.5, gam=8.0, tau=1.0;
  double pi=3.14159265358979323846;
  double fmax=10;
  FILE* file2;  
  char fname2[100];
  fprintf(filefprintf, "--!\tSource wavelet initialization  ...\n");
  fprintf(filefprintf, "--!\t                                   \n");

  sprintf(fname2,"Ricker_fwps_%.2lf_Hz.dat",fmax);
  file2= fopen( fname2, "w" );

  sig *= fmax;
  for(i=0; i<nt; i++ ) {

         at = i * dt;
         arg = -1.0*sqrt(2.0/pi)*sig*gam;
         arg = arg*(sig-2.0*sig*gam*(sig*at-tau)*(sig*at-tau));
         source[i] = 0.05*arg*exp(-gam*(sig*at-tau)*(sig*at-tau));
         fprintf(file2, " %lf %lf \n", i*dt, source[i]);

     }
		              
   fclose(file2);
   sprintf(fname2,"Rickerfren_fwps_%.2lf_Hz.dat",fmax);
   file2 = fopen(fname2,"w");
   double  aw=32.0, ac;
   int iw;
   for(iw=0;iw<nw;iw++){
          aw = iw*dw + 0.001;
	  fren[iw]=aw;
          at = 0.0;
          arg = aw*at;

          ac = 0.5; 
	  sw[2*iw]  =      ac*dt*cos(arg)*source[0];
	  sw[2*iw+1] = -1.0*ac*dt*sin(arg)*source[0];

          ac = 1.0;
	  for(i=1;i<nt-1;i++){
                at = i*dt;
                arg = aw*at;
		sw[2*iw] += ac*dt*cos(arg)*source[i];
		sw[2*iw+1] -= ac*dt*sin(arg)*source[i];
	   }

                at = (nt-1)*dt;
                ac = 0.5;
                arg = aw*at;
		sw[2*iw]  += ac*dt*cos(arg)*source[nt-1];
		sw[2*iw+1] -= ac*dt*sin(arg)*source[nt-1];
    
         fprintf(file2," %lf %lf %lf \n", fren[iw], sw[2*iw], sw[2*iw+1]);
	 }
	fclose(file2);

}

