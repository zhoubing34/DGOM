//
// Created by li12242 on 17/1/12.
//

#include "swe_input.h"
#include "swe_mesh.h"
#include "swe_driver2d.h"

#define DEBUG 0

/* ============ useful local variables ============ */
#define HEAD_TITLE "[============]"
#define HEAD_LINE   "[------------]"
#define ARG_LEN 100
#define FILE_NAME "swe_parameter.inc"
#define SECTION_NUM 6
/* ================================================ */

typedef struct{
    int info_linenum; ///< number of lines
    char *info_str; ///< strings of info
    int arg_num; ///< argument number
    char **arg_str; ///< argument str
} swe_arg_section;


static char** arg_create(int N){
    char **arg = (char**) calloc((size_t)N, sizeof(char*));
    arg[0] = (char*) calloc( (size_t)N*ARG_LEN, sizeof(char) );
    int i;
    for(i=1;i<N;++i){
        arg[i] = arg[i-1]+ ARG_LEN;
    }

    return arg;
}

static void arg_free(char **arg){
    free(arg[0]);
    free(arg);
}

static void print_section(swe_arg_section **arg){
#if DEBUG
    printf("into print_section\n");
#endif
    int n,i;
    for(n=0;n<SECTION_NUM;n++){
        swe_arg_section *sec_p = arg[n];
        printf("%s\n", sec_p->info_str);
        for(i=0;i<sec_p->arg_num;i++){
            printf("%s", sec_p->arg_str[i]);
        }
    }
#if DEBUG
    printf("out print_section\n");
#endif
    return;
}

static swe_arg_section** section_create(){
    swe_arg_section** section_p = (swe_arg_section**) calloc(SECTION_NUM, sizeof(swe_arg_section*));
    int n;
    for(n=0;n<SECTION_NUM;n++){
        section_p[n] = (swe_arg_section*) calloc(SECTION_NUM, sizeof(swe_arg_section));
    }

    // 0. section: case info
    int id = 0, var_num = 1;
    char caseinfo[] = HEAD_TITLE "DGOM: 2d swe solver\n"
            HEAD_LINE "case indicator (1 parameter)\n"
            HEAD_LINE "    1. case indicator  |-- 1. dam break wet\n"
            HEAD_LINE "                       |-- 2. dam break dry\n"
            HEAD_LINE "                       |-- 3. parabolic bowl\n"
            HEAD_LINE "                       |-- 4. user specific";
    section_p[id]->info_linenum = 6;
    section_p[id]->info_str = (char*) calloc(strlen(caseinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, caseinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    // 1. section: std cell info
    id++, var_num = 2;
    char scinfo[] = HEAD_TITLE "standard element info (2 parameters)\n"
            HEAD_LINE "    1. element type  |--- 0. triangle\n"
            HEAD_LINE "                     |--- 1. quadrilateral\n"
            HEAD_LINE "    2. order of polynomial";
    section_p[id]->info_linenum = 4;
    section_p[id]->info_str = (char*) calloc(strlen(scinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, scinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    // 2. section: mesh info
    id++, var_num = 3;
    char meshinfo[] = HEAD_TITLE "mesh info (2 parameters)\n"
            HEAD_LINE "    1. num of elements in x direction\n"
            HEAD_LINE "    2. num of elements in y direction\n"
            HEAD_LINE "    3. user specific mesh name";
    section_p[id]->info_linenum = 4;
    section_p[id]->info_str = (char*) calloc(strlen(meshinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, meshinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    // 3. section: time info
    id++, var_num = 2;
    char timeinfo[] = HEAD_TITLE "time info (2 parameters)\n"
            HEAD_LINE "    1. CFL number\n"
            HEAD_LINE "    2. final time";
    section_p[id]->info_linenum = 3;
    section_p[id]->info_str = (char*) calloc(strlen(timeinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, timeinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    // 4. section: physical info
    id++, var_num = 4;
    char physinfo[] = HEAD_TITLE "physical parameter (3 parameters)\n"
            HEAD_LINE "    1. gravity acceleration\n"
            HEAD_LINE "    2. hcrit -- minimum water depth of wet cell\n"
            HEAD_LINE "    3. manning coefficient for friction\n"
            HEAD_LINE "    4. file name of topography elevation";
    section_p[id]->info_linenum = 5;
    section_p[id]->info_str = (char*) calloc(strlen(physinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, physinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    // 5. section: LDG info
    id++, var_num = 1;
    char LDGinfo[] = HEAD_TITLE "LDG parameter for viscosity term (1 parameters)\n"
            HEAD_LINE "    1. C12";
    section_p[id]->info_linenum = 2;
    section_p[id]->info_str = (char*) calloc(strlen(LDGinfo)+1, sizeof(char));
    strcpy(section_p[id]->info_str, LDGinfo);
    section_p[id]->arg_num = var_num;
    section_p[id]->arg_str = arg_create(var_num);

    return section_p;
}


static void section_free(swe_arg_section** section_p){
    int n;
    for(n=0;n<SECTION_NUM;n++){
        free(section_p[n]->info_str);
        arg_free(section_p[n]->arg_str);
        free(section_p[n]);
    }
    free(section_p);
}

void swe_create_input(){
    swe_arg_section **arg_section = section_create();
    print_section(arg_section);

    FILE *fp = fopen(FILE_NAME, "w");

    swe_arg_section *sec_p=NULL;
    int n,m;
    for(n=0;n<SECTION_NUM;n++){
        sec_p = arg_section[n];
        fprintf(fp, "%s\n", sec_p->info_str);
        for(m=0;m<sec_p->arg_num;m++){
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
    section_free(arg_section);
}

swe_arg_section** swe_read_input(){
#if DEBUG
    printf("into swe_read_input\n");
#endif
    swe_arg_section **arg = section_create();
    swe_arg_section *sec_p = NULL;

    FILE *fp = fopen(FILE_NAME, "r");
    char buffer[ARG_LEN];
    int n,i;
    for(n=0;n<SECTION_NUM;n++){
        sec_p = arg[n];
#if DEBUG
        printf("arg[%d], line num=%d, arg num=%d\n", n, sec_p->info_linenum, sec_p->arg_num);
#endif
        for(i=0;i<sec_p->info_linenum;i++){ fgets(buffer, ARG_LEN, fp); }
        for(i=0;i<sec_p->arg_num;i++){
            fgets(buffer, ARG_LEN, fp);
            strcpy(sec_p->arg_str[i], buffer);
        }
    }

    fclose(fp);
    //print_section(arg);
    return arg;
}

swe_solver* swe_create_solver(){

    swe_arg_section **arg = swe_read_input();
    swe_solver *solver = (swe_solver*) calloc(1, sizeof(swe_solver));

    int procid;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);

    /// 0. section: case info
    swe_arg_section *sec_p = arg[0];
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
    char *meshfile = sec_p->arg_str[2];
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
            if(!procid){  printf(HEAD_LINE " mesh file: %s\n", meshfile); }
            solver->phys = swe_file_mesh(solver, meshfile);
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
    char *botfile=NULL;
    sscanf(sec_p->arg_str[0], "%lf\n", &(solver->gra));
    sscanf(sec_p->arg_str[1], "%lf\n", &(solver->hcrit));
    sscanf(sec_p->arg_str[2], "%lf\n", &(solver->roughness));

    if(!procid) printf(HEAD_LINE " gra = %f\n", solver->gra);
    if(!procid) printf(HEAD_LINE " hcrit = %f\n", solver->hcrit);
    if(!procid) printf(HEAD_LINE " manning factor = %f\n", solver->roughness);


    switch (solver->caseid){
        case swe_dambreakwet:
            solver->bot = swe_flat_topography(solver); break;
        case swe_dambreakdry:
            solver->bot = swe_flat_topography(solver); break;
        case swe_parabolicbowl:
            solver->bot = swe_parabolic_topography(solver); break;
        case swe_userset:
            sscanf(sec_p->arg_str[3], "%s\n", botfile);
            printf(HEAD_LINE " topography file = %s\n", botfile);
            solver->bot = swe_read_topography(solver, botfile);
            break;
    }
    /// 5. section: LDG info
    sec_p = arg[5];
    sscanf(sec_p->arg_str[0], "%lf\n", &(solver->c12));
    if(!procid) printf(HEAD_LINE " c12 = %f\n", solver->c12);
    section_free(arg);
    return solver;
}
