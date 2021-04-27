#include <stdio.h>
#include <stdlib.h>
#include <time.h>
void __fastscapeapi_MOD_fastscape_init(int *ierr);
void __fastscapeapi_MOD_fastscape_set_nx_ny(int* nx,int* ny,int *ierr);
void __fastscapeapi_MOD_fastscape_set_xl_yl(double* xl,double* yl,int *ierr);
void __fastscapeapi_MOD_fastscape_set_dt(double* dt,int *ierr);
void __fastscapeapi_MOD_fastscape_set_erosional_parameters(double* kf1_arr,double* kf2,double* m,double* n,double* kd1_arr,double* kd2,double* g1,double* g2, double* p_flow_dir_exp,int *ierr);
//void fastscape_set_precipitation_rate(double* preci_rate);
void __fastscapeapi_MOD_fastscape_set_marine_parameters(double* sealevel,double* poro1,double* poro2,double* z1,double* z2, double* ratio,double* L,double* kds1,double* kds2,int *ierr);
void __fastscapeapi_MOD_fastscape_set_bc(int* bc,int *ierr);
void __fastscapeapi_MOD_fastscape_init_h(double *h,int *ierr);
void __fastscapeapi_MOD_fastscape_setup(int *ierr);
void __fastscapeapi_MOD_fastscape_copy_h(double *h,int *ierr);
void __fastscapeapi_MOD_fastscape_set_u(double *u,int *ierr);
void __fastscapeapi_MOD_fastscape_execute_step(int *ierr);
void __fastscapeapi_MOD_fastscape_get_step(int* istep,int *ierr);
void __fastscapeapi_MOD_fastscape_debug(int *ierr);
void __fastscapeapi_MOD_fastscape_copy_total_erosion(double* etot,int *ierr);
void __fastscapeapi_MOD_fastscape_copy_erosion_rate(double* erate,int *ierr);
void __fastscapeapi_MOD_fastscape_copy_basement(double* b,int *ierr);
void __fastscapeapi_MOD_fastscape_vtk(double* etot, double* vexp,int *ierr);
void __fastscapeapi_MOD_fastscape_view(int *ierr);

int chkerr(int ierr);

int main(void){
    srand(time(NULL));
    int  nx,ny,istep,nstep,nn,nfreq,bc,nr_fields;
    // double* u[],ux[],uy[],h[],b[],etot[],erate[],a[],chi[],catchment[],sedflux[],sedflux_shore[];
    double kf1,kf2,m,n,kd1,kd2,g1,g2,preci_rate,p_flow_dir_exp;
    double xl,yl,dx,dy,dt;
    double sealevel, poro1, poro2, ratio, L, kds1, kds2, z1, z2;
    double vexp;
    int i,j,ij;
    int ierr;

    bc = 1010;
    nx=201;
    ny=201;
    nn=nx*ny;
    xl=200.e3;
    yl=200.e3;
    dx=xl/(nx-1);
    dy=yl/(ny-1);
    dt=1.e3;
    kf1=1.e-5;
    kf2=2.e-5;
    // !kf2=kf1;
    m=0.4e0;
    n=1.e0;
    kd1=1.e-2;
    kd2=5.e-2;
    kd1 = 0.e0;
    kd2 = 0.e0;
    // !kd2=kd1;
    g1=1.e0;
    g2=1.e0;
    p_flow_dir_exp = -2.e0;


    preci_rate = 1.e0; // precipitation rate
    sealevel = 0.e0;
    poro1 = 0.0e0;
    poro2 = 0.0e0;
    ratio = 0.5e0;
    L = 0.5e2;
    kds1 = 2.e2;
    kds2 = 1.e2;
    z1 = 1.e3;
    z2 = 1.e3;

    // allocate memory for working arrays
    //double *u creates pointer to start of array. (double*) typecasts the allocated array
    // malloc returns void pointer. Thus we need typecasting. malloc does not assign a value
    // using calloc(nn,sizeof(double)) would initialise array with 0 values.
    double *u = (double*)malloc( nn * sizeof( double ) );
    double *ux= (double*)malloc( nn * sizeof( double ) );
    double *uy= (double*)malloc( nn * sizeof( double ) );
    double *h= (double*)malloc( nn * sizeof( double ) );
    double *b= (double*)malloc( nn * sizeof( double ) );
    double *etot= (double*)malloc( nn * sizeof( double ) );
    double *erate= (double*)malloc( nn * sizeof( double ) );
    double *a= (double*)malloc( nn * sizeof( double ) );
    double *chi= (double*)malloc( nn * sizeof( double ) );
    double *catchment= (double*)malloc( nn * sizeof( double ) );
    double *sedflux= (double*)malloc( nn * sizeof( double ) );
    double *sedflux_shore= (double*)malloc( nn * sizeof( double ) );
    double *kf1_arr = (double*)malloc( nn * sizeof( double ) );
    double *kd1_arr =(double*)malloc( nn * sizeof( double ) );
    nr_fields=2;
    // double field[nr_fields][nn];
    // double **field = (double**)malloc(sizeof(double*)*nr_fields + sizeof(double*)*nr_fields*nn);
    // double **field = (double **)malloc(nr_fields * sizeof(double*));
    // for (i=0; i<nr_fields; i++){
    //      field[i] = (double *)malloc(nn * sizeof(double));
    //  }
     // printf("%d\n",field[0] );
     // printf("%d\n",field[1] );
     __fastscapeapi_MOD_fastscape_init(&ierr);chkerr(ierr);
     __fastscapeapi_MOD_fastscape_set_nx_ny(&nx,&ny,&ierr);chkerr(ierr);
     __fastscapeapi_MOD_fastscape_setup(&ierr);chkerr(ierr);
     // if (ierr != 0) {
     //   printf("Internal FastScape error occurred: error code = %d\n",ierr);
     //   printf("C driver taking action - returning with exit code 1\n");
     //   return(1);
     // }

     __fastscapeapi_MOD_fastscape_set_xl_yl(&xl,&yl,&ierr);chkerr(ierr);
     __fastscapeapi_MOD_fastscape_set_dt(&dt,&ierr);chkerr(ierr);
     __fastscapeapi_MOD_fastscape_set_marine_parameters(&sealevel, &poro1, &poro2, &z1, &z2, &ratio, &L, &kds1, &kds2, &ierr);chkerr(ierr);
     __fastscapeapi_MOD_fastscape_set_bc(&bc,&ierr);chkerr(ierr);
    for (i=0;i <nn; i++){
        h[i] = rand() % 10 + 0;
        u[i] = 0.0;
        kf1_arr[i] = kf1;
        kd1_arr[i] = kd1;
    }
    for(j=0;j<ny;j++){
        for(i=0;i<ny;i++){
            ij=(j)*nx+i;
            if (j<ny/2){
                h[ij]=h[ij]-200.e0;
            }
            if (j>ny/2){
                h[ij]=h[ij]+1000.e0;
            }
        }
    }
    __fastscapeapi_MOD_fastscape_set_erosional_parameters(kf1_arr,&kf2,&m,&n,kd1_arr,&kd2,&g1,&g2,&p_flow_dir_exp,&ierr);chkerr(ierr);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }

    for(j=0;j<ny;j++){
        for(i=0;i<ny;i++){
            ij=(j)*nx+i;
            u[ij]=0.0;
        }
    }

    __fastscapeapi_MOD_fastscape_init_h(h,&ierr);chkerr(ierr);
    __fastscapeapi_MOD_fastscape_set_u(u,&ierr);chkerr(ierr);
    __fastscapeapi_MOD_fastscape_copy_h(h,&ierr);chkerr(ierr);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }
    nstep=500;
    nfreq=100;
    istep=nstep;
    __fastscapeapi_MOD_fastscape_view(&ierr);chkerr(ierr);

    while (istep<=nstep){
        __fastscapeapi_MOD_fastscape_execute_step(&ierr);chkerr(ierr);
        __fastscapeapi_MOD_fastscape_get_step(&istep,&ierr);chkerr(ierr);
        __fastscapeapi_MOD_fastscape_debug(&ierr);chkerr(ierr);
        if (istep%nfreq==0){
            printf("%d\n",istep );
            __fastscapeapi_MOD_fastscape_copy_total_erosion(etot,&ierr);chkerr(ierr);
            __fastscapeapi_MOD_fastscape_copy_erosion_rate(erate,&ierr);chkerr(ierr);
            __fastscapeapi_MOD_fastscape_copy_h(h,&ierr);chkerr(ierr);
            __fastscapeapi_MOD_fastscape_copy_basement(b,&ierr);chkerr(ierr);
            // for (i=0;i <nn; i++){
            //     printf("Element #%d: %f\n",i,h[i]);
            // }
            vexp = 2.0;
            __fastscapeapi_MOD_fastscape_vtk(etot,&vexp,&ierr);chkerr(ierr);
            vexp = -2.0;
            __fastscapeapi_MOD_fastscape_vtk(erate,&vexp,&ierr);chkerr(ierr);
        }
    }
    // fastscape_copy_h(h);
    // for (i=0;i <nn; i++){
    //     printf("Element #%d: %f\n",i,h[i]);
    // }
    //  free memory again
    free(u);
    free(ux);
    free(uy);
    free(h);
    free(b);
    free(etot);
    free(erate);
    free(a);
    free(chi);
    free(catchment);
    free(sedflux);
    free(sedflux_shore);
    free(kf1_arr);
    free(kd1_arr);

}

int chkerr(int ierr){
  if (ierr != 0) {
    printf("Internal FastScape error occurred: error code = %d\n",ierr);
    printf("C driver taking action - returning with exit code 1\n");
    exit(1);
  }
}
