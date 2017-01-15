#include "swe_driver2d.h"

/* Private function */
void SWE_NodalHLL2d(SWE_Solver2d *solver, real hM, real hP,
                    real qnM, real qnP, real qvM, real qvP,
                    real *Fhn, real *Fqxn, real *Fqyn);

/**
 * @brief
 * Calculate the numerical flux for each node.
 *
 * @details
 * Calculation of the numerical flux, which is the approximation of dual flux term
 * \f$ \mathbf{n} \cdot \mathbf{F}(U^-, u^+) \f$.
 *
 * @param[SWE_Solver2d*] solver
 * @param[real] nx outward vector of x component
 * @param[real] ny outward vector of y component
 * @param[real] hM
 * @param[real] hP
 * @param[real] qxM
 * @param[real] qxP
 * @param[real] qyM
 * @param[real] qyP
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Fhs   | real*    | numerical flux of h function
 * Fqxs  | real*    | numerical flux of qx function
 * Fqys  | real*    | numerical flux of qy function
 *
 */
void SWE_NodalNumFlux2d(SWE_Solver2d *solver, real nx, real ny,
                        real hM, real hP, real qxM, real qxP, real qyM, real qyP,
                        real *Fhs, real *Fqxs, real *Fqys){
    /* rotation of conserve variable */
    real qnM,qnP,qvM,qvP;
    qnM =  nx*qxM + ny*qyM;
    qnP =  nx*qxP + ny*qyP;
    qvM = -ny*qxM + nx*qyM;
    qvP = -ny*qxP + nx*qyP;

    /* HLL numerical flux */
    real Fhn,Fqxn,Fqyn;
    SWE_NodalHLL2d(solver, hM, hP, qnM, qnP, qvM, qvP, &Fhn, &Fqxn, &Fqyn);

    /* inverse rotation */
    *Fhs  = Fhn;
    *Fqxs = nx*Fqxn - ny*Fqyn;
    *Fqys = ny*Fqxn + nx*Fqyn;
}

/**
 * @brief
 * Calculate the HLL flux function for each point.
 *
 * @details
 * HLL flux function is the approximation of dual flux term
 *
 * @param[in] solver SWE solver structure
 * @param[in] hM water depth at local node
 * @param[in] hP water depth at adjacent cell
 * @param[in] qnM
 * @param[in] qnP
 * @param[in] qvM
 * @param[in] qvP
 * @param[inout] Fhn   Numerical flux term for function h
 * @param[inout] Fqxn  Numerical flux term for function qx
 * @param[inout] Fqyn  Numerical flux term for function qy
 *
 *
 */
void SWE_NodalHLL2d(SWE_Solver2d *solver, real hM, real hP,
                    real qnM, real qnP, real qvM, real qvP,
                    real *Fhn, real *Fqxn, real *Fqyn){

    real EhM,EqnM,EqvM,EhP,EqnP,EqvP;
    real GhM,GqnM,GqvM,GhP,GqnP,GqvP;

    SWE_NodalFlux2d(solver, hM, qnM, qvM, &EhM, &EqnM, &EqvM, &GhM, &GqnM, &GqvM);
    SWE_NodalFlux2d(solver, hP, qnP, qvP, &EhP, &EqnP, &EqvP, &GhP, &GqnP, &GqvP);

    /* calculation of wave speed */
    real sM, sP;
    real gra  = (real) solver->gra;
    real hmin = (real) solver->hcrit;
    real us, cs;
    real unM,unP;

    if( (hM>hmin) & (hP>hmin) ){
        unM=qnM/hM;
        unP=qnP/hP;
        us = (real)(0.5*(unM + unP)   + sqrt(gra*hM)   - sqrt(gra*hP));
        cs = (real)(0.5*(sqrt(gra*hM) + sqrt(gra*hP) ) + 0.25*(unM - unP));

        sM = (real)min(unM-sqrt(gra*hM), us-cs);
        sP = (real)max(unP+sqrt(gra*hP), us+cs);
    }else if ( (hM>hmin) & (hP<=hmin) ){
        unM=qnM/hM;
        sM = (real)(unM -  sqrt(gra*hM) );
        sP = (real)(unM +2*sqrt(gra*hM) );
    }else if ( (hM<=hmin) & (hP>hmin) ){
        unP=qnP/hP;
        sM = (real)(unP -2*sqrt(gra*hP) );
        sP = (real)(unP +  sqrt(gra*hP) );
    }else{ /* both dry element */
        sM = 0; sP = 0;
    }

    /* HLL function */
    if ( (sM>=0) & (sP>0) ){
        *Fhn = EhM; *Fqxn = EqnM; *Fqyn = EqvM;
    }else if((sM<0) & (sP>0)){
        *Fhn  = (sP*EhM  - sM*EhP  + sM*sP*(hP  - hM ))/(sP - sM);
        *Fqxn = (sP*EqnM - sM*EqnP + sM*sP*(qnP - qnM))/(sP - sM);
        *Fqyn = (sP*EqvM - sM*EqvP + sM*sP*(qvP - qvM))/(sP - sM);
    }else if( (sM<0)&(sP<=0) ){
        *Fhn = EhP; *Fqxn = EqnP; *Fqyn = EqvP;
    }else if( (sM==0) & (sP==0) ){
        *Fhn = 0; *Fqxn = 0; *Fqyn = 0;
    }else{
        printf("Wrong status of numerical wave speed!\n");
        printf("hM  = %f, hP  = %f, hmin = %f\n", hM, hP, hmin);
        printf("qnM = %f, qnP = %f\n", qnM, qnP);
        printf("sM  = %f, sP  = %f\n", sM, sP);
        printf("us  = %f, cs  = %f\n", us, cs);
        exit(-2);
    }
}

/**
 * @brief
 * Calculation of the source term for element.
 *
 * @details
 * The source term
 * \f[S = \begin{bmatrix} 0 \cr
 * -gh\frac{\partial z}{\partial x} \cr
 * -gh\frac{\partial z}{\partial y} \end{bmatrix}\f]
 *
 * @param[in] phys
 * @param[in] solver
 * @param[in] vgeo  volume geometric coefficients of element: drdx, drdy, dsdx, dsdy
 * @param[in] Qk    variables in element
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * soureTerm| real*    | source term
 *
 */
void SWE_ElementalSource2d(PhysDomain2d *phys, SWE_Solver2d *solver,
                           int k, real *vgeo, real *Qk, real *soureTerm){

    MultiReg2d   *mesh  = phys->mesh;
    StdRegions2d *shape = mesh->stdcell;

    const real gra     = (real)solver->gra;
    const real *f_Dr   = shape->f_Dr;
    const real *f_Ds   = shape->f_Ds;
    const double *bot  = solver->bot[k];
    const double hcrit = solver->hcrit;
    const int Np       = shape->Np;
    const int Nfields  = phys->Nfields;

    int n,m,sk,geoid=0,ind=0;

    for (n=0;n<Np;++n) {

        const real *ptDr = f_Dr + n*Np;
        const real *ptDs = f_Ds + n*Np;
        const real drdx  = vgeo[geoid++], drdy = vgeo[geoid++];
        const real dsdx  = vgeo[geoid++], dsdy = vgeo[geoid++];
        const real h = Qk[n*Nfields];

        real Sh = 0, Sqx = 0, Sqy = 0;

        if(h>hcrit){ // for wet nodes
            for (sk = 0, m = 0; m < Np; ++m) {
                const real dr = ptDr[m];
                const real ds = ptDs[m];
                const real dx = drdx * dr + dsdx * ds;
                const real dy = drdy * dr + dsdy * ds;

                const real z = (real) bot[sk++];

                Sqx += dx * z;
                Sqy += dy * z;
            }
        }

        soureTerm[ind++] = Sh;
        soureTerm[ind++] = -gra*h*Sqx;
        soureTerm[ind++] = -gra*h*Sqy;
    }
    return;
}

/**
 * @brief
 * Calculation of flux terms for element.
 *
 * @details
 * The flux term
 * \f[ E = \begin{bmatrix}q_x \cr \frac{q_x^2}{h} + \frac{1}{2}gh^2 \cr
 * \frac{q_xq_y}{h}\end{bmatrix} \quad G = \begin{bmatrix}q_y \cr \frac{q_xq_y}{h} \cr
 * \frac{q_y^2}{h} + \frac{1}{2}gh^2 \end{bmatrix} \f]
 *
 * @param[in] phys
 * @param[in] solver
 * @param[in] Qk        variables in element
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Eflx     | real*    | flux in the x direction
 * Gflx     | real*    | flux in the y direction
 *
 */
void SWE_ElementalFlux2d(PhysDomain2d *phys, SWE_Solver2d *solver, real *Qk, real *Eflx, real *Gflx){

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

        SWE_NodalFlux2d(solver, h, qx, qy, &Eh, &Eqx, &Eqy, &Gh, &Gqx, &Gqy);

        Eflx[xk++] = Eh;
        Eflx[xk++] = Eqx;
        Eflx[xk++] = Eqy;

        Gflx[yk++] = Gh;
        Gflx[yk++] = Gqx;
        Gflx[yk++] = Gqy;
    }
    return;
}

/**
 * @brief
 * Calculate the flux term of each node.
 *
 * @details
 * The function returns the flux term on each node with
 * \f[ E = \begin{bmatrix} E_h \cr E_qx \cr E_qy \end{bmatrix} =
 * \begin{bmatrix} q_x \cr \frac{q_x^2}{h}+\frac{1}{2}gh^2 \cr \frac{q_xq_y}{h} \end{bmatrix}, \quad
 * G = \begin{bmatrix} G_h \cr G_qx \cr G_qy \end{bmatrix} =
 * \begin{bmatrix} q_y \cr \frac{q_xq_y}{h} \cr \frac{q_y^2}{h}+\frac{1}{2}gh^2 \end{bmatrix}\f]
 *
 * @param[SWE_Solver2d] solver
 * @param[real] h
 * @param[real] qx
 * @param[real] qy
 *
 * @return
 * return values:
 * name     | type     | description of value
 * -------- |----------|----------------------
 * Eh   | real*   | flux term of h  on x coordinate
 * Eqx  | real*   | flux term of qx on x coordinate
 * Eqy  | real*   | flux term of qy on x coordinate
 * Gh   | real*   | flux term of h  on y coordinate
 * Gqx  | real*   | flux term of qx on y coordinate
 * Gqy  | real*   | flux term of qy on y coordinate
 *
 */
void SWE_NodalFlux2d(SWE_Solver2d *solver,
                     real h, real qx, real qy,
                     real *Eh, real *Eqx, real *Eqy,
                     real *Gh, real *Gqx, real *Gqy){

    real hcrit = (real)solver->hcrit;
    real gra  = (real)solver->gra;

    if(h>hcrit){
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
    return;
}
