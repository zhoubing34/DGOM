#include "SWEDriver2d.h"

/* private function */
void SWEFlux(SWESolver *solver, real h,   real qx,   real qy,
             real *Eh, real *Eqx, real *Eqy, real *Gh, real *Gqx, real *Gqy);

void SWENumFlux2d(SWESolver *solver, real nx, real ny,
                  real hM, real hP, real qxM, real qxP, real qyM, real qyP,
                  real *Fhs, real *Fqxs, real *Fqys){
    /* rotation */
    real qnM,qnP,qvM,qvP;
    qnM =  nx*qxM + ny*qyM;
    qnP =  nx*qxP + ny*qyP;
    qvM = -ny*qxM + nx*qyM;
    qvP = -ny*qxP + nx*qyP;
    real EhM,EqxM,EqyM,EhP,EqxP,EqyP;
    real GhM,GqxM,GqyM,GhP,GqxP,GqyP;

    SWEFlux(solver, hM, qnM, qvM, &EhM, &EqxM, &EqyM, &GhM, &GqxM, &GqyM);
    SWEFlux(solver, hP, qnP, qvP, &EhP, &EqxP, &EqyP, &GhP, &GqxP, &GqyP);

    /* calculation of wave speed */
    real sM, sP;
    real gra  = (real) solver->gra;
    real hmin = (real) solver->hcrit;

    if( (hM>hmin) & (hP>hmin) ){
        real us, cs;
        real unM=qnM/hM,unP=qnP/hP;
        us = (real)(0.5*(unM + unP) + sqrt(gra*hM) - sqrt(gra*hP));
        cs = (real)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        sM = (real)min(unM-sqrt(gra*hM), us-cs);
        sP = (real)min(unM+sqrt(gra*hM), us+cs);
    }else if ( (hM>hmin) & (hP<=hmin) ){
        real unM=qnM/hM;
        sM = (real)(unM -  sqrt(gra*hM) );
        sP = (real)(unM +2*sqrt(gra*hM) );
    }else if ( (hM<=hmin) & (hP>hmin) ){
        real unP=qnP/hP;
        sM = (real)(unP -2*sqrt(gra*hP) );
        sP = (real)(unP +  sqrt(gra*hP) );
    }else{ /* both dry element */
        sM = 0; sP = 0;
    }

    /* HLL function */
    real Fhn,Fqxn,Fqyn;
    if (sM>0){
        Fhn = EhM; Fqxn = EqxM; Fqyn = EqyM;
    }else if((sM<=0) & (sP>=0)){
        Fhn  = sP*EhM  - sM*EhP  + sM*sP*(hP  - hM );
        Fqxn = sP*EqxM - sM*EqxP + sM*sP*(qxP - qxM);
        Fqyn = sP*EqyM - sM*EqyP + sM*sP*(qyP - qyM);
    }else if(sP<0){
        Fhn = EhP; Fqxn = EqxP; Fqyn = EqyP;
    }else{
        printf("Wrong status of numerical wave speed!\n");
        printf("hM = %f, hP = %f\n", hM, hP);
        printf("sM = %f, sP = %f\n", sM, sP);
        exit(-2);
    }

    /* inverse rotation */
    *Fhs  = Fhn;
    *Fqxs = nx*Fqxn - ny*Fqyn;
    *Fqys = ny*Fqxn + nx*Fqyn;
}

/**
 * @brief
 * Calculation of the source term.
 *
 * @details
 * The source term
 * \f[S = \begin{bmatrix} 0 \cr
 * -gh\frac{\partial z}{\partial x} \cr
 * -gh\frac{\partial z}{\partial y} \end{bmatrix}\f]
 *
 *
 * @param[in] phys
 * @param[in] solver
 * @param[in] vgeo
 * @param[in] Qk
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Sour   | real*  | source term
 *
 */
void SWESource(PhysDomain2d *phys, SWESolver *solver,
               int k, real *vgeo, real *Qk, real *Sour){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;
    real gra = (real)solver->gra;

    real *f_Dr    = shape->f_Dr;
    real *f_Ds    = shape->f_Ds;
    double *bot   = solver->bot[k];
    double hmin   = solver->hcrit;

    int n,m,geoid=0,isdry=0;
    real Sh=0, Sqx=0, Sqy=0;
    int Nfields  = phys->Nfields;


    for(n=0;n<shape->Np;n++){
        const real h = Qk[n*Nfields];
        if (h<hmin) {isdry = 1;}
    }

    if(!isdry) { /* wet element */
        for (n = 0; n < shape->Np; ++n) {

            const real *ptDr = f_Dr + n * shape->Np;
            const real *ptDs = f_Ds + n * shape->Np;

            const real drdx = vgeo[geoid++], drdy = vgeo[geoid++];
            const real dsdx = vgeo[geoid++], dsdy = vgeo[geoid++];

            const real h = Qk[n * Nfields];
            if (h < hmin) { isdry = 1; }

            int sk = 0;
            for (m = 0; m < shape->Np; ++m) {
                const real dr = ptDr[m];
                const real ds = ptDs[m];
                const real dx = drdx * dr + dsdx * ds;
                const real dy = drdy * dr + dsdy * ds;

                const real z = (real) bot[sk++];

                Sqx += -gra * h * dx * z;
                Sqy += -gra * h * dy * z;
            }

            int ind = n * Nfields;
            Sour[ind++] = Sh;
            Sour[ind++] = Sqx;
            Sour[ind++] = Sqy;
        }
    }else{  /* dry element */
        for (n = 0; n < shape->Np; ++n) {
            int ind = n * Nfields;
            Sour[ind++] = 0;
            Sour[ind++] = 0;
            Sour[ind++] = 0;
        }
    }
}

/**
 * @brief
 * Calculation of flux terms.
 *
 * @details
 * The flux term
 * \f[ E = \begin{bmatrix}q_x \cr \frac{q_x^2}{h} + \frac{1}{2}gh^2 \cr
 * \frac{q_xq_y}{h}\end{bmatrix} \quad G = \begin{bmatrix}q_y \cr \frac{q_xq_y}{h} \cr
 * \frac{q_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} \f]
 *
 *
 * @param[in] phys
 * @param[in] solver
 * @param[in] Np        number of points
 * @param[in] Qk        variables
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Eflx     | real*    | flux in the x direction
 * Gflx     | real*    | flux in the y direction
 *
 */
void SWEFlux2d(PhysDomain2d *phys, SWESolver *solver, real *Qk, real *Eflx, real *Gflx){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    real h,qx,qy;
    int i,sk=0,xk=0,yk=0;

    for(i=0;i<shape->Np;i++){
        h  = Qk[sk++];
        qx = Qk[sk++];
        qy = Qk[sk++];

        real Eh,Eqx,Eqy;
        real Gh,Gqx,Gqy;

        SWEFlux(solver, h, qx, qy, &Eh, &Eqx, &Eqy, &Gh, &Gqx, &Gqy);

        Eflx[xk++] = Eh;
        Eflx[xk++] = Eqx;
        Eflx[xk++] = Eqy;

        Gflx[yk++] = Gh;
        Gflx[yk++] = Gqx;
        Gflx[yk++] = Gqy;
    }

}

void SWEFlux(SWESolver *solver,
             real h,   real qx,   real qy,
             real *Eh, real *Eqx, real *Eqy,
             real *Gh, real *Gqx, real *Gqy){

    real hmin = (real)solver->hcrit;
    real gra  = (real)solver->gra;

    if(h>hmin){
        *Eh  = qx;
        *Eqx = (real)(qx*qx/h + 0.5*gra*h*h);
        *Eqy = qx*qy/h;
        *Gh  = qy;
        *Gqx = qx*qy/h;
        *Gqy = (real)(qy*qy/h + 0.5*gra*h*h);
    }else{
        *Eh  = 0; *Eqx = 0; *Eqy = 0;
        *Gh  = 0; *Gqx = 0; *Gqy = 0;
    }
}

/**
 * @brief
 * Predict the wave speed and give local delta time
 *
 * @details
 * The wave speed is given as \f[ c = \sqrt{gh}+\left| \mathbf{u} \right| \f],
 * and the delta time is derived \f[ dt = dx/dt \f], while \f$ dx \f$ is the
 * characteristic length.
 *
 * @param[in] phys
 * @param[in] solver
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * dt   | double      | delta time
 *
 */
double SWEPredictDt(PhysDomain2d *phys, SWESolver *solver, double CFL){
    double dt   = 1.0e4;
    double c    = 0;
    double gra  = solver->gra;
    double hmin = solver->hcrit;
    double dx   = solver->dx;
    double h,u,v,gdt;
    int k,i,isdry,ind;

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    int Nfields = phys->Nfields;
    for(k=0;k<mesh->K;k++){
        isdry = 0; /* initialization */
        c     = 0;
        for(i=0;i<shape->Np;i++){
            ind = (k*shape->Np +i)*Nfields;
            h   = phys->f_Q[ind++];

            if(h<hmin) { /* dry element, do nothing */
                isdry = 1;
            }else{       /* wet element, predict wave speed */
                u = phys->f_Q[ind++]; u /= h;
                v = phys->f_Q[ind++]; v /= h;
                c = max(c, sqrt(gra*h) + sqrt(u*u+v*v));
            }
        }

        if(!isdry){ /* only calculate the wet element */
            dt = min(dt, dx/c);
        }
    }

    /* gather all the delta time in all process */
    dt = CFL*dt;
    MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    return gdt;
}

