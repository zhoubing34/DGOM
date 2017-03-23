//
// Created by li12242 on 12/27/16.
//

#include "dg_phys_strong_surf_opt.h"

#define DEBUG 0

/**
 * @brief calculate the surface integral of strong form
 * @details
 * \f[ M^{-1} \cdot M_{es} (E \cdot nx + G \cdot ny - Fhs) \f]
 * where the \f[ Fhs \f] is the numerical flux, \f[ F = (E, G) \f] is the flux term.
 *
 * @param[in,out] phys pointer to dg_phys structure;
 * @param[in] slipwall_func boundary condition for slip-wall;
 * @param[in] nonslipwall_condition boundary condition for non-slip-wall;
 * @param[in] nodal_flux function to calculate the flux term;
 * @param[in] numerical_flux numerical flux function;
 *
 */
void dg_phys_strong_surf_opt2d(dg_phys *phys,
                               Wall_Condition slipwall_func,  // slip wall condition
                               Wall_Condition non_slipwall_func, // non-slip wall condition
                               Nodal_Flux_Fun nodal_flux, // flux term
                               Numerical_Flux numerical_flux) // numerical flux function
{
    dg_edge *edge = dg_phys_edge(phys);
    dg_cell *cell = dg_phys_cell(phys);
    const int Nfield = dg_phys_Nfield(phys);
    const int Nedge = dg_edge_Nedge(edge);
    const int Nfaces = dg_cell_Nfaces(cell);
    const int Nfptotal = dg_cell_Nfptotal(cell);
    const int Np = dg_cell_Np(cell);

    dg_real *f_Q = dg_phys_f_Q(phys);
    dg_real *f_inQ = dg_phys_f_recvQ(phys); // phys->f_recvQ;
    dg_real *f_ext = dg_phys_f_extQ(phys); // phys->f_extQ;
    dg_real *f_LIFT  = dg_cell_f_LIFT(cell); //phys->cell->f_LIFT;
    dg_real *f_rhsQ = dg_phys_f_rhsQ(phys); //phys->f_rhsQ;

    register int f,m,n,fld,surfid=0,nodeid=0;

    for(f=0;f<Nedge;f++){
        const int k1 = edge->surfinfo[surfid++];
        const int k2 = edge->surfinfo[surfid++];
        const int f1 = edge->surfinfo[surfid++];
        const int f2 = edge->surfinfo[surfid++];
        const int ftype = edge->surfinfo[surfid++];
        const int Nfp = dg_cell_Nfp(cell)[f1];

        dg_real flux_M[Nfp*Nfaces*Nfield];
        dg_real flux_P[Nfp*Nfaces*Nfield];
        int fp_M[Nfp], fp_P[Nfp];
        for(m=0;m<Nfp;m++){
            const int idM = (int)edge->nodeinfo[nodeid++];
            const int idP = (int)edge->nodeinfo[nodeid++];
            fp_M[m] = (int)edge->nodeinfo[nodeid++];
            fp_P[m] = (int)edge->nodeinfo[nodeid++];
            const dg_real nx = edge->nodeinfo[nodeid++];
            const dg_real ny = edge->nodeinfo[nodeid++];
            const dg_real fsc = edge->nodeinfo[nodeid++];

            dg_real f_M[Nfield], f_P[Nfield], Fhs[Nfield];
            // local face2d values
            for(fld=0;fld<Nfield;fld++) {f_M[fld] = f_Q[idM*Nfield+fld];}

            // adjacent nodal values
            switch (ftype){
                case FACE_INNER:
                    for(fld=0;fld<Nfield;fld++){ f_P[fld] = f_Q[idP*Nfield+fld]; }
                    break;
                case FACE_PARALL:
                    for(fld=0;fld<Nfield;fld++){ f_P[fld] = f_inQ[idP*Nfield+fld]; }
                    break;
                case FACE_SLIPWALL:
                    slipwall_func(nx, ny, f_M, f_P);
                    break;
                case FACE_NSLIPWALL:
                    non_slipwall_func(nx, ny, f_M, f_P);
                    break;
                default:
                    for(fld=0;fld<Nfield;fld++){ f_P[fld] = f_ext[idP*Nfield+fld]; }
                    break;
            }

            dg_real E_M[Nfield], G_M[Nfield];
            dg_real E_P[Nfield], G_P[Nfield];

            nodal_flux(f_M, E_M, G_M);
            nodal_flux(f_P, E_P, G_P);
            numerical_flux(nx,ny,f_M,f_P,Fhs);

            dg_real *t = flux_M + m*Nfield;
            dg_real *s = flux_P + m*Nfield;
            for(fld=0;fld<Nfield;fld++){
                t[fld] = fsc*(  nx*E_M[fld] + ny*G_M[fld] - Fhs[fld] );
                s[fld] = fsc*( -nx*E_P[fld] - ny*G_P[fld] + Fhs[fld] );
            }

        }

        for(n=0;n<Np;n++){
            const dg_real *ptLIFT = f_LIFT + n*Nfptotal;
            dg_real *f_rhsM = f_rhsQ + Nfield*(n+k1*Np);

            for(m=0;m<Nfp;m++){
                const int col1 = fp_M[m];
                const dg_real L = ptLIFT[col1];
                const dg_real *t = flux_M+m*Nfield;
                for(fld=0;fld<Nfield;fld++)
                    f_rhsM[fld] += L*t[fld];
            }
            if (ftype == FACE_INNER){
                dg_real *f_rhsP = f_rhsQ + Nfield*(n+k2*Np);
                for(m=0;m<Nfp;m++){
                    const int col2 = fp_P[m];
                    const dg_real L = ptLIFT[col2];
                    const dg_real *s = flux_P+m*Nfield;
                    for(fld=0;fld<Nfield;fld++) {f_rhsP[fld] += L*s[fld];}
                }
            }
        }

    }

    return;
}