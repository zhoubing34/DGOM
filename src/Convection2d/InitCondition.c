#include "ConvectionDriver2d.h"

/**
 * @brief
 * Set initial condition
 *
 * @details
 * Set initial scalar field distribution and the constant flow field
 * 1. scalar distribution
 * 2. constant velocity field
 *
 */
double InitCondition(PhysDomain2d * phys, PhysDomain2d *flowRate){
    double sigma = 125*1e3/(33*33);
    double xc, yc, t;
    double w, dt = 1e6;

    MultiReg2d *mesh = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;
    const int K = mesh->K;
    int k, n, f, sk=0;


    /* initial position */
    xc = 0.0; yc = 0.6;

    /* initial scalar field */
    for(k=0;k<K;++k){
        for(n=0;n<shape->Np;++n){
            t = -sigma * ( ( mesh->x[k][n] - xc )*( mesh->x[k][n] - xc )
                           + ( mesh->y[k][n] - yc )*( mesh->y[k][n] - yc ) );
            phys->f_Q[sk++] = (float) exp(t);
        }
    }

    /* flow rate field */
    w = 5*M_PI/6;
    sk=0;
    for (k=0; k<K; ++k){
        for (n=0;n<shape->Np; ++n){
            flowRate->f_Q[sk++] = (float)(-w * mesh->y[k][n]); // flow rate at x-coordinate
            flowRate->f_Q[sk++] = (float)(w * mesh->x[k][n]);  // flow rate at y-coordinate
        }
    }

    /* time step */
    double ri, dscal, cfl = 0.5;
    int Nfp = shape->Nfp;
    int Nfaces = shape->Nfaces;

    for(k=0;k<mesh->K;++k){
        for(f=0;f<Nfaces;++f){
            for(n=0;n<Nfp;n++) {
                int surfid = k*shape->Np + shape->Fmask[f][n];
                double J = mesh->J[surfid];
                int index = n * 6 + f * Nfp * 6 + k * Nfp * Nfaces * 6;
                double sJ = phys->surfinfo[index + 3]; /* edge length */
                dt = min(dt, J/sJ);
            }
        }
    }

    double gdt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    dt = cfl*gdt/((shape->N+1)*(shape->N+1));

//    sk = 0;
//    for(k=0;k<mesh->K;++k){
//        ri = 1.e6; dscal=1e6;
//        for(f=0;f<Nfaces;++f){
//            for(n=0;n<Nfp;n++) {
//                int index = n * 6 + f * Nfp * 6 + k * Nfp * Nfaces * 6;
//                double sJ = phys->surfinfo[index + 3]; /* edge length */
//                ri = min(ri, sJ * 2 / shape->Nfp);
//
//            }
//        }
//
//        for(n=0;n<shape->Np;n++){
//            double u = flowRate->f_Q[sk++];
//            double v = flowRate->f_Q[sk++];
//            double flows = sqrt(u*u + v*v);
//            dscal = min(mesh->ciradius[k]/flows, dscal);
//
//        }
//
//        dt = min(dt, 2.0/3.0*cfl*dscal*ri);
//    }

    return dt;
}

/**
 * @brief
 * Print the mesh information of each process in files
 *
 * @details
 *
 * @author
 * li12242, Tianjin University, li12242@tju.edu.cn
 *
 */
#define DSET_NAME_LEN 1024
void PrintPhys ( PhysDomain2d *phys, char *name ){
    int n, m, fld, rank, nprocs, ret, sk, sf;
    char filename[DSET_NAME_LEN];

    MultiReg2d *mesh = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    rank = mesh->procid;
    nprocs = mesh->nprocs;

    ret = snprintf(filename, DSET_NAME_LEN,"%s.%d-%d.txt",name,rank,nprocs);
    if (ret >= DSET_NAME_LEN) {
        fprintf(stderr, "name too long \n");
        exit(-1);
    }

    FILE *fig = fopen(filename, "w");

    fprintf(fig, "field data: \n");
    fprintf(fig, "\n Nfields = %d\n", phys->Nfields);
    fprintf(fig, "\n parNtotalout = %d\n", phys->parNtotalout);

    for( fld = 0; fld<phys->Nfields; fld++) {
        fprintf(fig, "\n field:%d \n", fld);
        sk = fld;
        for (n = 0; n < mesh->K; n++) {
            for (m = 0; m < shape->Np; m++) {
                fprintf(fig, "%f,\t", phys->f_Q[sk]);
                sk += 2;
            }
            fprintf(fig, "\n");
        }
    }

#if 0

    sk = 0;
    fprintf(fig, "\n Node vgeo = \n");
    for (n = 0; n < mesh->K; n ++){
        fprintf(fig, "Ele %d:\n", n);
        for( m = 0; m < shape->Np; m++)
            fprintf(fig, "drdx = %f, drdy = %f, dsdx = %f, dsdy = %f\n",
                    mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++], mesh->vgeo[sk++]);
        fprintf(fig, "\n");
    }
#endif
    fclose(fig);

}
