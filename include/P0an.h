void greenfunction(double *r1, double *r2, double k,  double* G0);
void greenfunction2(double x, double z, double k,  double* G00);
void greenfunctiontable(FILE *filefprintf, int nx2, int nz2, double dx, double dz, int nw, double* fren, double c0,double* G0r,double* G0i,
double* VG0r,double* VG0i);
void prodcomp( double *z1, double *z2, double *zres );
void greenfunction3(double x, double z, double k,  double* G0);
void P0num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz, double *sourcefren, int *ps,  double *G0r,double *G0i, double *P0r,double *P0i);
void P1num(FILE *filefprintf, int nw, int nx, int nz,double dx,double dz,double c0, int *ps, double *sourcefren,double* fren, double *VG0r, double *VG0i,double* vpe, double *P0r,double *P0i, double *P1r,double *P1i);
void Pwtot(int nw, int nx,  int nz,   double dw,double dx, double dz, double dt, int nt,  double* fren, double *P1r, double *P1i, double *Pret);
void greenfunctionV(double x, double z, double k,  double* G00);
