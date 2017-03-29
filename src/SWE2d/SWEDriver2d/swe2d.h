//
// Created by li12242 on 16/7/30.
//

#ifndef DGOM_SWEDRIVER2D_H
#define DGOM_SWEDRIVER2D_H

#include "../SWELib/swe_lib.h"
/* swe_input.c */
void swe_input(int argc, char **argv);
arg_section** swe_read_section();
void swe_free_section(arg_section **section_p);
/* swe_init.c */
void swe_init();

#endif //DGOM_SWEDRIVER2D_H
