//
// Created by li12242 on 17/1/12.
//

#include "swe_input.h"
#include "swe_mesh.h"
#include "Utility/arg_section.h"

#define DEBUG 0

/* ============ useful local variables ============ */
#define HEAD_TITLE "[============]"
#define HEAD_LINE   "[------------]"
#define FILE_NAME "swe_parameter.inc"
#define SECTION_NUM 6
/* ================================================ */

static arg_section** swe_create_section(){
    int section_num = SECTION_NUM;
    arg_section** section_p = (arg_section**) calloc((size_t)section_num, sizeof(arg_section*));

    int ind = 0;
    /// 0. section: case info
    char caseinfo[] = HEAD_TITLE "DGOM: 2d swe solver\n"
            HEAD_LINE "case indicator (1 parameter)\n"
            HEAD_LINE "    1. case indicator  |-- 1. dam break wet\n"
            HEAD_LINE "                       |-- 2. dam break dry\n"
            HEAD_LINE "                       |-- 3. parabolic bowl\n"
            HEAD_LINE "                       |-- 4. user specific\n";
    int var_num = 1;
    section_p[ind++] = section_create(caseinfo, var_num);
    /// 1. section: std cell info
    char scinfo[] = HEAD_TITLE "standard element info (2 parameters)\n"
            HEAD_LINE "    1. element type  |--- 0. triangle\n"
            HEAD_LINE "                     |--- 1. quadrilateral\n"
            HEAD_LINE "    2. order of polynomial\n";
    var_num = 2;
    section_p[ind++] = section_create(scinfo, var_num);
    /// 2. section: mesh info
    char meshinfo[] = HEAD_TITLE "mesh info (2 parameters)\n"
            HEAD_LINE "    1. num of elements in x direction\n"
            HEAD_LINE "    2. num of elements in y direction\n"
            HEAD_LINE "    3. user specific case name\n";
    var_num = 3;
    section_p[ind++] = section_create(meshinfo, var_num);
    /// 3. section: time info
    char timeinfo[] = HEAD_TITLE "time info (2 parameters)\n"
            HEAD_LINE "    1. CFL number\n"
            HEAD_LINE "    2. final time\n";
    var_num = 2;
    section_p[ind++] = section_create(timeinfo, var_num);
    /// 4. section: physical info
    char physinfo[] = HEAD_TITLE "physical parameter (3 parameters)\n"
            HEAD_LINE "    1. gravity acceleration\n"
            HEAD_LINE "    2. hcrit -- minimum water depth of wet cell\n"
            HEAD_LINE "    3. manning coefficient for friction\n";
    var_num = 3;
    section_p[ind++] = section_create(physinfo, var_num);
    /// 5. section: LDG info
    char LDGinfo[] = HEAD_TITLE "LDG parameter for viscosity term (1 parameters)\n"
            HEAD_LINE "    1. C12\n";
    var_num = 1;
    section_p[ind++] = section_create(LDGinfo, var_num);

    return section_p;
}

static void swe_free_section(arg_section **section_p){
    int n;
    for(n=0;n<SECTION_NUM;n++)
        section_free(section_p[n]);

    free(section_p);
    return;
}

void swe_create_input(){
    arg_section **section_p = swe_create_section();

    int n;
    for(n=0;n<SECTION_NUM;n++)
        section_print(section_p[n]);

    FILE *fp = fopen(FILE_NAME, "w");
    for(n=0;n<SECTION_NUM;n++){
        section_write_file(section_p[n], fp);
    }
    fclose(fp);
    swe_free_section(section_p);
}

arg_section** swe_read_input(){

    arg_section **arg = swe_create_section();
    FILE *fp = fopen(FILE_NAME, "r");
    int n;
    for(n=0;n<SECTION_NUM;n++){
        section_read_file(arg[n], fp);
    }

    fclose(fp);
    //print_section(arg);
    return arg;
}

swe_solver* swe_create_solver(){

    arg_section **arg = swe_read_input();
    swe_solver *solver = (swe_solver*) calloc(1, sizeof(swe_solver));

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    /// 0. section: case info
    arg_section *sec_p = arg[0];
    sscanf(sec_p->arg_str[0], "%d\n", &(solver->caseid));

    if(!procid){
        if( (solver->caseid==swe_dambreakdry) || (solver->caseid==swe_dambreakwet)
            || (solver->caseid==swe_parabolicbowl) || (solver->caseid==swe_userset) ){
            printf(HEAD_TITLE " DGOM: 2d swe solver\n"
                           HEAD_LINE " case = %d\n", solver->caseid);
        }else{
            fprintf(stderr, "%s (%d): Unknown case id %d\n",
                    __FILE__, __LINE__, solver->caseid);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    /// 1. section: cell info
    sec_p = arg[1];
    sscanf(sec_p->arg_str[0], "%d\n", &(solver->celltype));
    sscanf(sec_p->arg_str[1], "%d\n", &(solver->N));

    if(!procid){
        if( solver->celltype==TRIANGLE ){
            printf(HEAD_LINE " cell type: triangle\n");
        }else if(solver->celltype==QUADRIL){
            printf(HEAD_LINE " cell type: quadrilateral\n");
        }else{
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FILE__, __LINE__, solver->celltype);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    /// 2. section: mesh info
    sec_p = arg[2];
    int Mx, My;
    double xmin, xmax, ymin, ymax;
    sscanf(sec_p->arg_str[0], "%d\n", &(Mx));
    sscanf(sec_p->arg_str[1], "%d\n", &(My));
    solver->casename = calloc(strlen(sec_p->arg_str[2]), sizeof(char));
    strcpy(solver->casename, sec_p->arg_str[2]);
    char *casename = solver->casename;
    switch (solver->caseid){
        case swe_dambreakwet:
            xmin=0, xmax=1000, ymin=-100, ymax=100;
            if(!procid){ printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
            break;
        case swe_dambreakdry:
            xmin=0, xmax=1000, ymin=-100, ymax=100;
            if(!procid){ printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
            break;
        case swe_parabolicbowl:
            xmin = -4000, xmax = 4000, ymin = -4000, ymax = 4000;
            if(!procid){  printf(HEAD_LINE " cell num: Mx=%d, My=%d\n", Mx, My); }
            solver->phys = swe_uniform_mesh(solver, Mx, My, xmin, xmax, ymin, ymax);
            break;
        case swe_userset:
            if(!procid){  printf(HEAD_LINE " case name: %s\n", casename); }
            solver->phys = swe_file_mesh(solver, casename);
            break;
        default:
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
    /// 3. section: time info
    sec_p = arg[3];
    sscanf(sec_p->arg_str[0], "%lf\n", &(solver->cfl));
    sscanf(sec_p->arg_str[1], "%lf\n", &(solver->ftime));
    if(!procid) printf(HEAD_LINE " cfl = %f\n", solver->cfl);
    if(!procid) printf(HEAD_LINE " final time = %f\n", solver->ftime);
    /// 4. section: phys info
    sec_p = arg[4];
    sscanf(sec_p->arg_str[0], "%lf\n", &(solver->gra));
    sscanf(sec_p->arg_str[1], "%lf\n", &(solver->hcrit));
    sscanf(sec_p->arg_str[2], "%lf\n", &(solver->roughness));

    if(!procid) printf(HEAD_LINE " gra = %f\n", solver->gra);
    if(!procid) printf(HEAD_LINE " hcrit = %f\n", solver->hcrit);
    if(!procid) printf(HEAD_LINE " manning factor = %f\n", solver->roughness);

    /// 5. section: LDG info
    sec_p = arg[5];
    sscanf(sec_p->arg_str[0], "%lf\n", &(solver->c12));
    if(!procid) printf(HEAD_LINE " c12 = %f\n", solver->c12);

    swe_free_section(arg);
    return solver;
}
