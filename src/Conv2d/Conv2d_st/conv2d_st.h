#ifndef CONVECTION2D_ST_H
#define CONVECTION2D_ST_H

#include "Utility/nc_library.h"
#include "Utility/arg_section.h"
#include "MultiArea/Grid/dg_grid_reader.h"
#include "../ConvLib/conv_lib2d.h"

#define SEC_NUM 6
#define HEADEND   "[============]"
#define HEADLINE  "[------------]"

typedef enum{
    Conv_Rotation = 0,
    Conv_AdvDiff = 1,
    Conv_Adv = 2,
} Conv_St_Case;

/* conv_input */
void conv_input_st(int argc, char *agrv[]);
void conv_arg_section_free(arg_section **section_p);
arg_section** conv_read_inputfile(char *filename);

/* conv_init */
void conv_init_st();

/* conv_extsol */
void conv_normerr();

#endif