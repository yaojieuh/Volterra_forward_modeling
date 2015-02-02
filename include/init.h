void init_acou_homog(FILE *filefprintf, double dx, double dz, int nx, int nz, double c0, double *vel, double *vpe);
void init_acou_layer(FILE *filefprintf, double dx, double dz, int nx, int nz, double c0, double *vel, double *vpe);
void init_source_ricker_fwps(FILE *filefprintf, int nt, double dt,int nw, double dw, double *source, double *fren, double *sw);
