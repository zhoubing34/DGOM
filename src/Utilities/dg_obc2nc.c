//
// Created by li12242 on 17/2/21.
//
#include "Utility/utility.h"
#include "Utility/nc_library.h"

nc_file * create_ncfile(char *file_name);
void write_obc2nc(char *file_name, nc_file *obcfile);

int main(int argc, char **argv){

    int Nobc = argc-1;
    int i;
    MPI_Init(&argc, &argv);
    nc_file **obcfile = (nc_file**) calloc(Nobc, sizeof(nc_file*));
    for(i=0;i<Nobc;i++){
        obcfile[i] = create_ncfile(argv[i+1]);
        write_obc2nc(argv[i+1], obcfile[i]);
        nc_file_close(obcfile[i]);
        nc_file_free(obcfile[i]);
    }
    free(obcfile);

    MPI_Finalize();
    return 0;
}

void write_obc2nc(char *file_name, nc_file *obcfile){
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

    int fid = obcfile->id;
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

    vector_double_free(p_Q);
    fclose(fp);
    return;
}

nc_file * create_ncfile(char *file_name){

    char *ncfile_name = (char*) calloc(MAX_NAME_LENGTH, sizeof(char));
    strcpy(ncfile_name, file_name);
    strcat(ncfile_name, ".nc");

    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* read obc info */
    FILE *fp;
    dg_fopen(fp, file_name, "Unable to open obc file");
    int Nv, Nfield;
    if( fscanf(fp, "%d", &Nv)!=1 ){ // reading Nv
        fprintf(stderr, "%s (%d)\nerror in reading Nv\n", __FUNCTION__, __LINE__);
        exit(-1);
    };
    if( fscanf(fp, "%d", &Nfield)!=1 ){ // reading Nfield
        fprintf(stderr, "%s (%d)\nerror in reading Nfield\n", __FUNCTION__, __LINE__);
        exit(-1);
    };
    /* reading vertex list */
    int *vertlist = vector_int_create(Nv);
    int i;
    for(i=0;i<Nv;i++){
        if(fscanf(fp, "%d", vertlist+i)!=1){
            fprintf(stderr, "%s (%d)\nerror in reading vertex %d\n",
                    __FUNCTION__, __LINE__, i);
            exit(-1);
        }
        vertlist[i] -= 1; // change to C type
    }

    /* creating dimensions */
    nc_dim *nfld = nc_dim_create("Nfield", Nfield);
    nc_dim *nv = nc_dim_create("Nv", Nv);
    nc_dim *nt = nc_dim_create("Nt", 0);
    /* creating variables */
    int ndim = 3;
    nc_dim **dimarray = (nc_dim**) calloc(ndim, sizeof(nc_dim*));
    dimarray[0] = nt;
    dimarray[1] = nv;
    dimarray[2] = nfld; /* the inner loop dimension comes the last */
    nc_var *time = nc_var_create("time", 1, dimarray, NC_DOUBLE);
    nc_var *vert = nc_var_create("vert", 1, dimarray+1, NC_INT);
    nc_var *data = nc_var_create("obc", 3, dimarray, NC_DOUBLE);
    /* creating nc files */
    int nvar = 3;
    nc_var **vararray = (nc_var**) calloc(nvar, sizeof(nc_var*));
    vararray[0] = vert;
    vararray[1] = time;
    vararray[2] = data;
    nc_file *obcfile = nc_file_create(file_name, procid, nprocs, ndim, dimarray, nvar, vararray);
    /* rename the file */
    free(obcfile->name);
    obcfile->name = ncfile_name;
    /* create nc file */
    nc_file_init(obcfile);

    /* assignment */
    ncmpi_put_var_int_all(obcfile->id, obcfile->var_vec_p[0]->id, vertlist);

    vector_int_free(vertlist);
    fclose(fp);
    return obcfile;
}