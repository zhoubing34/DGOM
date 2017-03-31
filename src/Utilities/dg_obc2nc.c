//
// Created by li12242 on 17/2/21.
//
#include "Utility/utility.h"
#include "Utility/nc_library.h"
#include "Utility/unit_test.h"

NC_File * create_ncfile(char *file_name);
void write_obc2nc(char *file_name, NC_File *obcfile);

/**
 * @brief
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    NC_File *obcfile = create_ncfile(argv[1]);
    write_obc2nc(argv[1], obcfile);

    nc_file_close(obcfile);
    nc_file_free(obcfile);

    MPI_Finalize();
    return 0;
}

/**
 * @brief
 * @param file_name
 * @param obcfile
 */
void write_obc2nc(char *file_name, NC_File *obcfile){
    FILE *fp;
    dg_fopen(fp, file_name, "Unable to open obc file");
    char buffer[MAX_NAME_LENGTH];
    fgets(buffer, MAX_NAME_LENGTH, fp);
    fgets(buffer, MAX_NAME_LENGTH, fp);
    int Nv = obcfile->dim_vec_p[1]->len;
    int Nfield = obcfile->dim_vec_p[2]->len;

    double *p_Q = vector_double_create(Nv*Nfield);
    double time;

    MPI_Offset start_t, count_t;
    start_t = 0; count_t = 1;
    MPI_Offset start_v[3] = {start_t, 0, 0};
    MPI_Offset count_v[3] = {1, Nv, Nfield};
    int i,fld,tmp, fend=0;

    int fid = obcfile->ncid;
    int time_id = obcfile->var_vec_p[1]->id;
    int obc_id = obcfile->var_vec_p[2]->id;
    while(1){
        for(i=0;i<Nv;i++){
            if(fscanf(fp, "%d", &tmp) != 1){ fend = 1; break; }
            fscanf(fp, "%lf", &time);
            for(fld=0;fld<Nfield;fld++){
                fscanf(fp, "%lf", p_Q+i*Nfield+fld);
            }
        }
        if(fend){ break; } // read file to the end, exit;
        /* put time variable */
        ncmpi_put_vara_double_all(fid, time_id, &start_t, &count_t, &time);
        /* put obc variable */
        ncmpi_put_vara_double_all(fid, obc_id, start_v, count_v, p_Q);
        start_t++;
        start_v[0] = start_t;
    }
    printf(HEADLINE "finish reading obc data\n");
    printf(HEADLINE "Nt = %lld\n", start_t);
    printf(HEADEND "\n\n");
    vector_double_free(p_Q);
    fclose(fp);
    return;
}
/**
 * @brief
 * @param file_name
 * @return
 */
NC_File* create_ncfile(char *file_name){

    char *ncfile_name = (char*) calloc(MAX_NAME_LENGTH, sizeof(char));
    strcpy(ncfile_name, file_name);
    strcat(ncfile_name, ".nc");

    /* read obc info */
    FILE *fp;
    dg_fopen(fp, file_name, "Unable to open obc file");
    printf(HEADSTART "\n");
    printf(HEADLINE "open file: %s\n", file_name);
    int Nv, Nfield;
    if( fscanf(fp, "%d", &Nv)!=1 ){ // reading Nv
        fprintf(stderr, "%s (%d): Error in reading Nv\n", __FUNCTION__, __LINE__);
        exit(-1);
    };
    if( fscanf(fp, "%d", &Nfield)!=1 ){ // reading Nfield
        fprintf(stderr, "%s (%d): Error in reading Nfield\n", __FUNCTION__, __LINE__);
        exit(-1);
    };
    printf(HEADLINE "reading dimensions:\n");
    printf(HEADLINE "Nv = %d\n", Nv);
    printf(HEADLINE "Nfield = %d\n", Nfield);
    /* reading vertex list */
    int *vertlist = vector_int_create(Nv);
    int i;
    for(i=0;i<Nv;i++){
        if(fscanf(fp, "%d", vertlist+i)!=1){
            fprintf(stderr, "%s (%d)\nError occurs in reading vertex %d\n",
                    __FUNCTION__, __LINE__, i);
            exit(-1);
        }
        vertlist[i] -= 1; // change to C type
    }

    /* creating dimensions */
    NC_Dim *nfld = nc_dim_create("Nfield", Nfield);
    NC_Dim *nv = nc_dim_create("Nv", Nv);
    NC_Dim *nt = nc_dim_create("Nt", 0);
    /* creating variables */
    int ndim = 3;
    NC_Dim **dimarray = (NC_Dim**) calloc(ndim, sizeof(NC_Dim*));
    dimarray[0] = nt;
    dimarray[1] = nv;
    dimarray[2] = nfld; /* the inner loop dimension comes the last */
    NC_Var *time = nc_var_create("time", 1, dimarray, NC_DOUBLE);
    NC_Var *vert = nc_var_create("vert", 1, dimarray+1, NC_INT);
    NC_Var *data = nc_var_create("obc", 3, dimarray, NC_DOUBLE);
    /* creating nc files */
    int nvar = 3;
    NC_Var **vararray = (NC_Var**) calloc(nvar, sizeof(NC_Var*));
    vararray[0] = vert;
    vararray[1] = time;
    vararray[2] = data;
    NC_File *obcfile = nc_file_create(file_name, ndim, dimarray, nvar, vararray);
    /* rename the file */
    free(obcfile->name);
    obcfile->name = ncfile_name;
    /* create nc file */
    nc_file_define(obcfile);
    printf(HEADLINE "generate NetCDF file: %s\n", ncfile_name);

    /* assignment */
    ncmpi_put_var_int_all(obcfile->ncid, obcfile->var_vec_p[0]->id, vertlist);

    vector_int_free(vertlist);
    fclose(fp);
    return obcfile;
}