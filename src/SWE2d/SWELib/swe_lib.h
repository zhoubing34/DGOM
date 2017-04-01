//
// Created by li12242 on 17/3/28.
//

#ifndef DGOM_SWE2D_LIB_H
#define DGOM_SWE2D_LIB_H

#include "PhysField/dg_phys.h"
#include "Utility/nc_library.h"
#include "Utility/utility.h"
#include "Utility/arg_section.h"
#include "PhysField/dg_phys_strong_vol_opt.h"
#include "PhysField/dg_phys_strong_surf_opt.h"

/** input file marco */
#define HEAD_LINE   "[------------]"
#define SECTION_NUM 6 ///< number section in input file;
#define NFIELD 4 ///< number of physical field;

typedef enum {
    SWE_HELP = 0,
    SWE_CREATE_INPUT = 1,
    SWE_RUN = 2,
    SWE_UNKNOWN = 3,
}SWE_Run_Type;

typedef enum{
    ZERO_GRADIENT = 4,
    CLAMPED_ALL = 5,
    CLAMPED_H = 6,
    CLAMPED_FLOW = 7,
    FLATHER_H = 8,
    FLATHER_FLOW = 9,
} SWE_OBC;

typedef struct{
    char filename[MAX_NAME_LENGTH]; /// input file;
    dg_phys *phys; ///< physical field;
    double ftime; ///< final time;
    double dt; ///< user-set time interval;
    double outDt; ///< output time interval;
    double cfl; ///< CFL number;
    /* physical parameters */
    double gra; ///< gravity acceleration;
    double hcrit; ///< minimum water depth;
    double *m; ///< manning roughness coefficient for friction term;
    char outfilename[MAX_NAME_LENGTH]; ///< name of output file;
    NC_File *outfile; ///< NetCDF output file;
}SWE_Solver;

/* swe_output.c */
void swe_output();
void swe_save_var(int tstep, double time);
/* swe_run.c */
void swe_run();
/* swe_flux.c */
int swe_flux_term(dg_real *var, dg_real *Eflux, dg_real *Gflux);
int swe_hll_flux(dg_real nx, dg_real ny, dg_real *varM, dg_real *varP, dg_real *Fhs);
/* swe_rhs.c */
void swe_rhs(dg_phys *phys, dg_real frka, dg_real frkb, dg_real fdt);
/* swe_obc.c */
int swe_obc(dg_real nx, dg_real ny, dg_real *f_M, dg_real *f_ext, int obc_ind, dg_real *f_P);

#endif //DGOM_SWE2D_LIB_H
