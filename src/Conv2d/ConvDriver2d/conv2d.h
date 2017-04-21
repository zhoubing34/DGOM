#ifndef CONVECTION2D_H
#define CONVECTION2D_H

#include "Utility/nc_library.h"
#include "Utility/arg_section.h"
#include "MultiArea/Grid/dg_grid_reader.h"
#include "../ConvLib/conv_lib2d.h"

#define SEC_NUM 6
#define HEADEND   "[============]"
#define HEADLINE  "[------------]"

/* conv_input */
void conv_input(int argc, char *agrv[]);
void conv_arg_section_free(arg_section **section_p);
arg_section** conv_read_inputfile(char *);
/* conv_init */
void conv_init();

#endif