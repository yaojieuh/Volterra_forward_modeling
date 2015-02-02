//!########################################################################
//!                                                                       #
//! Copyright (C) University .                     #
//! This file is a modeling for 2d acoustic media  #
//!                                                                       #
//! 
//! Feb. 2014                                                      # 
//!                                                                       #
//!########################################################################
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include "init.h"
#include "P0an.h"
//#include "Reflectiondata.h"

int main( int argc, char *argv[] )
{
  // grids dimensions 
  int nz    = 51; 
  int nx    = 101;
  int ndims = nx*nz;       // total grid points number



  // physical grid size
  double dz = 5.0;     // modeling grid spacial step in z
  double dx = 5.0;    // modeling grid spacial step in x

 
  FILE* file1;
  char fname1[100];
  int ret;
  sprintf(fname1,"test_2D_acoustic.out");
  file1 = fopen(fname1,"w");


  fprintf(file1, "--!\t                                     \n");
  fprintf(file1, "--!\t2D acoustic Volterra modeling research code\n");
  fprintf(file1, "--!\t                                     \n");

  ret = fflush(file1);

  // velocity field
  double c0=1500; //reference velocity
  double* vel  = malloc( ndims*sizeof(double) );
  double* vpe  = malloc( ndims*sizeof(double) );
  //init_acou_homog( file1, dim_w, vel);
  init_acou_layer(file1,  dx,dz,nx,nz, c0, vel,vpe);
                                                                                                                                                                            

 
   // time and frequency discretization parameter
  int   nt = 1000;        // number of steps
  double dt = 0.001;      // time step
  int nw=101;
  double dw=2;
  double* sourcet = malloc(nt*sizeof(double) );
  double* sourcefren  = malloc( (2*nw)*sizeof(double) );
  double* fren  = malloc( nw*sizeof(double) );
  init_source_ricker_fwps( file1, nt, dt,nw,dw, sourcet,fren, sourcefren); 


  int* ps = malloc(2*sizeof(int));
  double* pscor = malloc(2*sizeof(double));
  ps[1] = 0;
  ps[0] = (nx-1)/2;
  pscor[1] = (ps[0])*dz;
  pscor[0] = (ps[1])*dx;   

   int nz2    = 2*nz+1; 
   int nx2    = 2*nx+1;
   double* G0r = malloc(nz2*nx2*nw*sizeof(double) );
   double* G0i = malloc(nz2*nx2*nw*sizeof(double) );
   double* VG0r = malloc(nz2*nx2*nw*sizeof(double) );
   greenfunctiontable(file1,nx2, nz2, dx,dz, nw,  fren, c0,G0r, G0i,VG0r);

   /*
   double* P0r = malloc(nz*nx*nw*sizeof(double) );
   double* P0i = malloc(nz*nx*nw*sizeof(double) );
   P0num(file1,  nw, nx, nz, dx,dz,sourcefren, ps,  G0r,G0i, P0r,P0i);

   double* P1r = malloc(nz*nx*nw*sizeof(double) );
   double* P1i = malloc(nz*nx*nw*sizeof(double) );
   
   P1num(file1,  nw, nx,  nz,dx, dz,c0, ps, sourcefren, fren, VG0r,  vpe, P0r,P0i, P1r,P1i);

   nt=401;

   double* Pret = malloc(nz*nx*nt*sizeof(double) );
   
   Pwtot(nw,nx, nz,  dw, dx, dz, dt, nt,fren,P1r, P1i, Pret);
   // Pwtot(nw,nx,  nz,  dw, dx, dz, dt, nt,   fren,P0r, P0i, Pret);
 

   double* precor = malloc(2*nx*sizeof(double));
   int* prec = malloc(2*nx*sizeof(int));
   

        double pi = 3.14159265358979323846;
   	int nk=nx,ik,iw,iz2;
   	int nk2=(nk-1)/2;
   	double dk=2*pi/(nx-1)/dx;
	double r0,k,w,coef=0;
	double ac,arg,k1,k2,q,z1,kr,ag;
	double P0tempr, P0tempi;
	double* P0Vr = malloc(nk*nz*sizeof(double));
   	double* P0Vi = malloc(nk*nz*sizeof(double));
	double* P0kr = malloc(nk*sizeof(double));
   	double* P0ki = malloc(nk*sizeof(double));
   	double* P0rewsr = malloc(nx*nz*sizeof(double));
   	double* P0rewsi = malloc(nx*nz*sizeof(double));	
        double* P0numr = malloc(nw*nx*sizeof(double));
        double* P0numi = malloc(nw*nx*sizeof(double));
	double* P0numsr = malloc(nw*nx*nz*sizeof(double));
        double* P0numsi = malloc(nw*nx*nz*sizeof(double));
	double* P0numlr = malloc(nw*nx*nz*sizeof(double));
        double* P0numli = malloc(nw*nx*nz*sizeof(double));
   	double* P0rewtr = malloc(nw*nx*nz*sizeof(double));
  	double* P0rewti = malloc(nw*nx*nz*sizeof(double));
  FILE* file2;  
  char fname2[100];
        sprintf(fname2,"green.dat");
        file2 = fopen(fname2,"w");
	
	for (iz=1;iz<nz-300;iz++){		
		z=iz*dz;
		fprintf(file1,"Depth  %lf   \n",z);
		for(ix=0;ix<nx;ix++){
    			prec[2*ix] =  ps[0]+iz;
    			prec[2*ix+1] =  ix;
    			precor[2*ix] = (prec[2*ix])*dz;
    			precor[2*ix+1] = (prec[2*ix+1])*dx;
   
   		}	
		P0num(nx, nw,  c0,sourcefren, pscor, precor, P0numr, P0numi);
		for(ix=0;ix<nx;ix++){
			for(iw=0;iw<nw;iw++){
				P0numsr[iw*nx*nz+iz*nx+ix]=P0numr[ix*nw+iw];
				P0numsi[iw*nx*nz+iz*nx+ix]=P0numi[ix*nw+iw];
			}	
		}
	}	

   	for (iw=0;iw<nw;iw++){
	        fprintf(file1,"Frequency  %lf   \n",iw*dw);
		w=iw*dw+0.001;
		k=w/c0;
		
		for (iz=0;iz<nz;iz++){	
			for (ik=0;ik<nk;ik++){		
				P0Vr[iz*nk+ik]=0;
				P0Vi[iz*nk+ik]=0;
			}
			for (ix=0;ix<nx;ix++){		
				P0rewsr[iz*nx+ix]=0;
				P0rewsi[iz*nx+ix]=0;
			}
			
		}		
		
		for (iz=1;iz<nz-300;iz++){
			z=iz*dz;
			for (ik=0;ik<nk;ik++){
				kr=(ik-nk2)*dk;
				P0Vr[(iz-1)*nk+ik]=0;
				P0Vi[(iz-1)*nk+ik]=0;
				for (ix=0;ix<nx;ix++){
					ag=kr*precor[2*ix+1];
					P0Vr[(iz-1)*nk+ik]+=(cos(ag)*P0rewsr[(iz-1)*nx+ix]+sin(ag)*P0rewsi[(iz-1)*nx+ix])*v[(iz-1)*nx+ix];
					P0Vi[(iz-1)*nk+ik]+=(cos(ag)*P0rewsi[(iz-1)*nx+ix]-sin(ag)*P0rewsr[(iz-1)*nx+ix])*v[(iz-1)*nx+ix];
				}
				P0Vr[(iz-1)*nk+ik]*=dx/2/pi;
				P0Vi[(iz-1)*nk+ik]*=dx/2/pi;
							
			}
			for(ik=0;ik<nk;ik++){
				kr=(ik-nk2)*dk;
				k2=k*k-kr*kr;
				P0kr[ik]=0;
				P0ki[ik]=0;
				if (k2>0.5*k*k){
				//if (k2>=0){
					q=sqrt(k2);
					ag=q*z;		
					P0tempr=0;
					P0tempi=0;
					for (iz2=0;iz2<iz;iz2++){
						z1=iz2*dz;
						k1=k*k/q;
						P0tempr=sin(q*(z-z1))*P0Vr[iz2*nk+ik]*k1*dz;
						P0tempi=sin(q*(z-z1))*P0Vi[iz2*nk+ik]*k1*dz;
						P0kr[ik]+=P0tempr;
						P0ki[ik]+=P0tempi;	
					}
						
				}
			
					
			}
			
			for (ix=0;ix<nx;ix++){
				x=precor[2*ix+1];			
				P0rewsr[iz*nx+ix]=0;
				P0rewsi[iz*nx+ix]=0;
				for(ik=0;ik<nk;ik++){
					kr=(ik-nk2)*dk;
					k2=k*k-kr*kr;						
					if (k2>0){
						ag=kr*x;
						P0rewsr[iz*nx+ix]+=cos(ag)*P0kr[ik]-sin(ag)*P0ki[ik];
						P0rewsi[iz*nx+ix]+=cos(ag)*P0ki[ik]+sin(ag)*P0kr[ik];
					}
				}
				P0rewsr[iz*nx+ix]*=dk;
				P0rewsi[iz*nx+ix]*=dk;
				//P0rewsr[iz*nx+ix]=coef*(sourcefren[3*iw+1]*sin(k*r0-pi/4)+sourcefren[3*iw+2]*cos(k*r0-pi/4));
				//P0rewsi[iz*nx+ix]=coef*(-sourcefren[3*iw+1]*cos(k*r0-pi/4)+sourcefren[3*iw+2]*sin(k*r0-pi/4));
				P0rewsr[iz*nx+ix]+=P0numsr[iw*nx*nz+iz*nx+ix];
				P0rewsi[iz*nx+ix]+=P0numsi[iw*nx*nz+iz*nx+ix];
				P0rewtr[iw*nx*nz+iz*nx+ix]=	P0rewsr[iz*nx+ix];
				P0rewti[iw*nx*nz+iz*nx+ix]=	P0rewsi[iz*nx+ix];
				//fprintf(file2," %lf  %lf  %lf \n",  x,iz*dz,P0rewsr[iz*nx+ix]);
							
			}
			
			//fprintf(file2,"  \n");
		}
		
	}
    
	double* Pswr1 = malloc(nw*nx*sizeof(double));
        double* Pswi1 = malloc(nw*nx*sizeof(double));
        double* Prets = malloc(nt*nx*sizeof(double));
        for (ix=0;ix<nx;ix++){
                for (iw=0;iw<nw;iw++){
                        //Pswr1[ix*nw+iw]=P0rewtr[iw*nx*nz+90*nx+ix];
                        //Pswi1[ix*nw+iw]=P0rewti[iw*nx*nz+90*nx+ix];
			  Pswr1[ix*nw+iw]=P0numsr[iw*nx*nz+90*nx+ix];
                          Pswi1[ix*nw+iw]=P0numsi[iw*nx*nz+90*nx+ix];
			fprintf(file2," %lf  %lf  %lf \n",  ix*dx,iw*dw,Pswr1[ix*nw+iw]);
                }
		fprintf(file2,"  \n");
        }
	

        //P0num(nx, nw,  c0,sourcefren, pscor, precor, Pswr1, Pswi1);
        Pwtot(nx,  nw, dw, dt, nt, precor, Pswr1, Pswi1, Prets);
*/
}
