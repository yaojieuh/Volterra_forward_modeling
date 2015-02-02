#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void greenfunction(double *r1, double *r2, double k,  double* G0)
{
	double r, coe;
	double pi = 3.14159265358979323846; 
	r=sqrt((r1[0]-r2[0])*(r1[0]-r2[0])+(r1[1]-r2[1])*(r1[1]-r2[1]));
	coe=sqrt(1.0/(pi*k*r*8));
	G0[0]=coe*sin(k*r-pi/4);
	G0[1]=coe*cos(k*r-pi/4);
	
}
void prodcomp( double *z1, double *z2, double *zres ){

     // real part
     zres[0] = z1[0]*z2[0]-z1[1]*z2[1]; 
     // imaginary part  
     zres[1] = z1[1]*z2[0]+z1[0]*z2[1]; 

}


void P0num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz, double *sourcefren, int *ps,  double *G0r,double *G0i, double *P0r,double *P0i){
	int iw,ix,iz;
        int xindex,zindex;
   	double G0[2],sw[2],po[2];
  	double x,z;
 	int nz2    = 2*nz+1; 
   	int nx2    = 2*nx+1;

	fprintf(filefprintf, "--!\tP0 calculation  ...\n");
  	fprintf(filefprintf, "--!\t   		\n");

	FILE* file2;  
  	char fname2[100];
  	
        for(iw=0;iw<nw;iw++){
		//omega=fren[iw];
		sw[0]=sourcefren[2*iw];
		sw[1]=sourcefren[2*iw+1];
  		for(ix=0;ix<nx;ix++){
			xindex=ix-ps[0]+nx;
			for(iz=0;iz<nz;iz++){
				zindex=iz-ps[1]+nz;
				G0[0]=G0r[iw*(nx2*nz2)+xindex*nz2+zindex];
				G0[1]=G0i[iw*(nx2*nz2)+xindex*nz2+zindex];
				prodcomp(sw, G0, po );
				P0r[iw*(nx*nz)+ix*nz+iz]=po[0];
				P0i[iw*(nx*nz)+ix*nz+iz]=po[1];
			}
        		
  		}
		if((iw%50)==0){
		sprintf(fname2,"P0num_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=(iz-(nz-1)/2)*dz;		
                          	fprintf(file2," %f %f %.12lf %.12lf \n", x, z, P0r[iw*(nx*nz)+ix*nz+iz],P0i[iw*(nx*nz)+ix*nz+iz]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
		
	}
}






void greenfunction2(double x, double z, double k,  double* G00)
{
	double pi = 3.14159265358979323846; 
	double Exp[2],exp2[2],exp3[2];
	double kz2,kn;
	int n,i;
	G00[0]=0;
	G00[1]=0;
	double L=2000;
	for (i=0;i<151;i++){
		n=i-75;
		kz2=k*k-(2*pi*n/L)*(2*pi*n/L);
		Exp[0]=cos(2*pi*n*x/L);
		Exp[1]=sin(2*pi*n*x/L);
		if (kz2>0){
			kn=sqrt(kz2);			
			exp2[0]=cos(kn*fabs(z));
			exp2[1]=sin(kn*fabs(z));
			prodcomp(Exp, exp2, exp3 );
			G00[0]+=exp3[1]/kn/2/L;
			G00[1]+=exp3[0]/kn/2/L;//sign issue!!!
		}
		else{
			kn=sqrt(-kz2);
			exp2[0]=exp(-kn*fabs(z));
			exp2[1]=0;
			prodcomp(Exp, exp2, exp3 );
			G00[0]-=exp3[0]/kn/2/L;
			G00[1]-=exp3[1]/kn/2/L;
		}
		
	}
	
}

void greenfunction3(double x, double z, double k,  double* G0)
{
	double r, coe;
	double pi = 3.14159265358979323846; 
	r=sqrt(x*x+z*z);
	if (r==0){
		G0[0]=0;
		G0[1]=-0.25;
	}else{
	coe=sqrt(1.0/(pi*k*r*8));
	G0[0]=coe*sin(k*r-pi/4);
	G0[1]=coe*cos(k*r-pi/4);
	}
	
}

void greenfunctionV(double x, double z, double k,  double* G00)
{
	double pi = 3.14159265358979323846; 
	double Exp[2],exp2[2],exp3[2];
	double kz2,kn;
	int n,i;
	G00[0]=0;
	G00[1]=0;
	double L=2000;
	for (i=0;i<151;i++){
		n=i-75;
		kz2=k*k-(2*pi*n/L)*(2*pi*n/L);
		Exp[0]=cos(2*pi*n*x/L);
		Exp[1]=sin(2*pi*n*x/L);
		if (kz2>0){
			kn=sqrt(kz2);			
			exp2[0]=sin(kn*fabs(z));
			exp2[1]=0;
			prodcomp(Exp, exp2, exp3 );
			G00[0]+=exp3[0]/kn/L;
			G00[1]+=exp3[1]/kn/L;
		}
		
		
	}
	
}
void greenfunctiontable(FILE *filefprintf, int nx2, int nz2, double dx, double dz, int nw, double* fren, double c0,double* G0r,double* G0i,
double* VG0r,double* VG0i)
{
	int ix,iz,iw;
	double omega, k;
        double x,z;
        double G00[2];

	fprintf(filefprintf, "--!\tGreen's Function calculation  ...\n");
  	fprintf(filefprintf, "--!\t                                   \n");
	
	FILE* file2;  
  	char fname2[100];
        for(iw=0;iw<nw;iw++){
		fprintf(filefprintf, "--!\t  Frequency_%d     \n",iw);
		omega=fren[iw];
                k=omega/c0;
	        for (ix=0;ix<nx2;ix++){
			x=(ix-(nx2-1)/2)*dx;
			for (iz=0;iz<nz2;iz++){
				z=(iz-(nz2-1)/2)*dz;
				//greenfunction2(x, z, k, G00);
				greenfunction3(x, z, k, G00);
				G0r[iw*(nx2*nz2)+ix*nz2+iz]=G00[0];
				G0i[iw*(nx2*nz2)+ix*nz2+iz]=G00[1];
				if(z>0){
					greenfunctionV(x, z, k, G00);
					//VG0r[iw*(nx2*nz2)+ix*nz2+iz]=G0r[iw*(nx2*nz2)+ix*nz2+iz]-G00[0];
					//VG0r[iw*(nx2*nz2)+ix*nz2+iz]=G0i[iw*(nx2*nz2)+ix*nz2+iz]-G00[1];
					VG0r[iw*(nx2*nz2)+ix*nz2+iz]=G00[0];
					VG0i[iw*(nx2*nz2)+ix*nz2+iz]=G00[1];
				}else{
					VG0r[iw*(nx2*nz2)+ix*nz2+iz]=0;
					VG0i[iw*(nx2*nz2)+ix*nz2+iz]=0;
				}
			}
		}

		
		if((iw%50)==0){
		sprintf(fname2,"Green_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx2;ix++){
			x=(ix-(nx2-1)/2)*dx;
			for (iz=0;iz<nz2;iz++){
				z=(iz-(nz2-1)/2)*dz;		
                          	fprintf(file2," %f %f %.12lf %.12lf %.12lf  %.12lf\n", x, z, G0r[iw*(nx2*nz2)+ix*nz2+iz],G0i[iw*(nx2*nz2)+ix*nz2+iz],VG0r[iw*(nx2*nz2)+ix*nz2+iz],VG0i[iw*(nx2*nz2)+ix*nz2+iz]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}	
	
}


void P1num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz,double c0, int *ps, double *sourcefren,double* fren, double *VG0r, double *VG0i, double* vpe, double *P0r,double *P0i, double *P1r,double *P1i){
	int iw,ix,iz,ix2,iz2;
        int xindex,zindex;
   	double VG0[2],Prtemp,Pitemp;
	double omega,k;
  	double x,z;
 	int nz2    = 2*nz+1; 
   	int nx2    = 2*nx+1;

	fprintf(filefprintf, "--!\tP1 calculation  ...\n");
  	fprintf(filefprintf, "--!\t   		\n");

	FILE* file2;  
  	char fname2[100];
  	
        for(iw=0;iw<nw;iw++){
		omega=fren[iw];
                k=omega/c0;
		
		for(iz=0;iz<nz;iz++){			
  			for(ix=0;ix<nx;ix++){
				P1r[iw*(nx*nz)+ix*nz+iz]=P0r[iw*(nx*nz)+ix*nz+iz];
				P1i[iw*(nx*nz)+ix*nz+iz]=P0i[iw*(nx*nz)+ix*nz+iz];				
				Prtemp=0;
				Pitemp=0;
				for(iz2=0;iz2<iz;iz2++){
					zindex=iz-iz2+nz;
					for(ix2=0;ix2<nx;ix2++){
						xindex=ix-ix2+nx;
						VG0[0]=VG0r[iw*(nx2*nz2)+xindex*nz2+zindex];
						VG0[1]=VG0i[iw*(nx2*nz2)+xindex*nz2+zindex];
						Prtemp+=vpe[ix2*nz+iz2]*(VG0[0]*P1r[iw*(nx*nz)+ix2*nz+iz2]-VG0[1]*P1i[iw*(nx*nz)+ix2*nz+iz2]);
						Pitemp+=vpe[ix2*nz+iz2]*(VG0[0]*P1i[iw*(nx*nz)+ix2*nz+iz2]+VG0[1]*P1r[iw*(nx*nz)+ix2*nz+iz2]);
					}
				}
				Prtemp*=dx*dz*k*k;
				Pitemp*=dx*dz*k*k;
				P1r[iw*(nx*nz)+ix*nz+iz]+=Prtemp;
				P1i[iw*(nx*nz)+ix*nz+iz]+=Pitemp;	
						

			}
		}

		if((iw%50)==0){
		sprintf(fname2,"P1_iw_%d.dat", iw);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf %.12lf \n", x, z, P1r[iw*(nx*nz)+ix*nz+iz],P1i[iw*(nx*nz)+ix*nz+iz]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
}



void Pwtot(int nw, int nx,  int nz,  double dw, double dx,double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret){
	int ix,iw,it,iz;
	double pi = 3.14159265358979323846; 
	double t,x,z, arg=0;
        for(iz=0;iz<nz;iz++){
	for(ix=0;ix<nx;ix++){
		x=(ix-(nx-1)/2)*dx;
		for(it=0;it<nt;it++){
			t=it*dt;
			Pret[iz*nt*nx+ix*nt+it]=0;
			for (iw=0;iw<nw;iw++){
				arg = fren[iw]*t;		
				Pret[iz*nt*nx+ix*nt+it]+=cos(arg)*P1r[iw*(nx*nz)+ix*nz+iz]-sin(arg)*P1i[iw*(nx*nz)+ix*nz+iz];
			}
			Pret[iz*nt*nx+ix*nt+it]*=dw*2/2/pi;
			
		}
		
		
	}
	}
	FILE* file2;  
  	char fname2[100];
	for(it=0;it<nt;it++){
		if((it%100)==0){
		sprintf(fname2,"snapshot_%d.dat", it);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (iz=0;iz<nz;iz++){
				z=iz*dz;		
                          	fprintf(file2," %f %f %.12lf \n", x, z, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}

	for(iz=0;iz<nz;iz++){
		if((iz%10)==0){
		sprintf(fname2,"trace_%d.dat", iz);
		file2 = fopen(fname2,"w");
 		for (ix=0;ix<nx;ix++){
			x=(ix-(nx-1)/2)*dx;
			for (it=0;it<nt;it++){
				t=it*dt;		
                          	fprintf(file2," %f %f %.12lf \n", x, t, Pret[iz*nt*nx+ix*nt+it]);
			}
                          fprintf(file2," \n");
		}

       		 fclose(file2);
		}
	}
	
}
