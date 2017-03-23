#include "dg_grid_BS.h"
#include "dg_grid_reader.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

/**
 * @brief Generation of uniform triangle grid.
 * @details
 * Generate uniform triangle mesh with specific coordinate range and number of elements.
 * The flag 'type' defines how the square divides into two triangles.
 * @param[in] cell standard cell object;
 * @param[in] Mx,My number of cells on x and y coordinate;
 * @param[in] xmin,xmax  range of x coordinate;
 * @param[in] ymin,ymax  range of y coordinate;
 * @param[in] type  flag for triangle dividing, 0 for '\' and 1 for '/';
 *
 * @return grid geometry grid object.
 */
dg_grid* dg_grid_uniform_tri(dg_cell *cell, int Mx, int My,
                             double xmin, double xmax,
                             double ymin, double ymax, int type) {
    /* check stand element */
    if(dg_cell_celltype(cell) !=  TRIANGLE){
        fprintf(stderr, "%s (%s):%d\nthe input element type %d is not triangle!\n",
               __FUNCTION__, __FILE__, __LINE__, dg_cell_celltype(cell));
        exit(-1);
    }

    /* parameter */
    const int K = Mx*My*2;
    const int Nv = (Mx+1)*(My+1);

    /* allocation */
    int **EToV = matrix_int_create(K, 3);
    double VX[Nv], VY[Nv];

    /* vertex coordinate */
    /* the order of vertex is from left to right, and from bottom to the top */
    int ind; /* vertex index */
    int dim1, dim2;
    for (dim1 = 0; dim1<My+1; dim1++){
        for (dim2 = 0; dim2<Mx+1; dim2++){
            ind = dim1*(Mx+1) + dim2;
            VX[ind] = (double)dim2/Mx * (xmax - xmin) + xmin;
            VY[ind] = (double)dim1/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int top[2], bot[2];
    for (dim1 = 0; dim1 < My; dim1++){
        for (dim2 = 0; dim2 < Mx; dim2++){
            /* Assignment of vertex index, which start from 0 */
            bot[0] = dim1*(Mx + 1)+dim2;     /* bottom left  */
            bot[1] = dim1*(Mx + 1)+dim2 + 1; /* bottom right */
            top[0] = bot[0] + (Mx + 1); /* top left  */
            top[1] = bot[1] + (Mx + 1); /* top right */

            /* Assignment of EToV */
            if ( type==1 ){ /* for '/' division */
                EToV[dim1 * Mx * 2 + dim2][0] = bot[0];
                EToV[dim1 * Mx * 2 + dim2][1] = top[1];
                EToV[dim1 * Mx * 2 + dim2][2] = top[0];

                EToV[dim1 * Mx * 2 + Mx + dim2][0] = bot[0];
                EToV[dim1 * Mx * 2 + Mx + dim2][1] = bot[1];
                EToV[dim1 * Mx * 2 + Mx + dim2][2] = top[1];
            }else if( type==0 ){ /* for '\' division */
                EToV[dim1 * Mx * 2 + dim2][0] = bot[1];
                EToV[dim1 * Mx * 2 + dim2][1] = top[1];
                EToV[dim1 * Mx * 2 + dim2][2] = top[0];

                EToV[dim1 * Mx * 2 + Mx + dim2][0] = bot[0];
                EToV[dim1 * Mx * 2 + Mx + dim2][1] = bot[1];
                EToV[dim1 * Mx * 2 + Mx + dim2][2] = top[0];
            }
        }
    }

    /* initialize */
    dg_grid* grid = dg_grid_create(cell, K, Nv, VX, VY, NULL, EToV);

    /* free memory */
    matrix_int_free(EToV);
    return grid;
}

/**
 * @brief Generation of uniform quadrilateral mesh.
 * @details
 * Generate uniform quadrilateral mesh with specific coordinate range and number of elements.
 * @param[in] cell dg_cell structure;
 * @param[in] Mx,My number of cells on x and y coordinate;
 * @param[in] xmin,xmax range of x coordinate;
 * @param[in] ymin,ymax range of y coordinate;
 * @return grid dg_grid structure.
 */
dg_grid* dg_grid_uniform_quad(dg_cell *cell, int Mx, int My,
                              double xmin, double xmax,
                              double ymin, double ymax)
{
    /* check stand element */
    if(dg_cell_celltype(cell) != QUADRIL){
        fprintf(stderr, "%s (%s):%d\nthe input element type %d is not quadrilateral!\n",
                __FUNCTION__, __FILE__, __LINE__, dg_cell_celltype(cell));
        exit(-1);
    }

    /* parameter */
    const int K = Mx*My;
    const int Nv = (Mx+1)*(My+1);

    /* allocation */
    int **EToV = matrix_int_create(K, 4);
    double VX[Nv], VY[Nv];

    /* vertex coordinate */
    int ind; /* vertex index */
    int dim1, dim2;
    for (dim1=0; dim1<My+1; dim1++){
        for (dim2=0; dim2<Mx+1; dim2++){
            /* the order of vertex is from left to right, and from bottom to the top */
            ind = dim1*(Mx+1) + dim2;
            VX[ind] = (double)dim2/Mx * (xmax - xmin) + xmin;
            VY[ind] = (double)dim1/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int top[2], bot[2];
    for (dim1 = 0; dim1 < My; dim1++){
        for (dim2 = 0; dim2 < Mx; dim2++){
            /* Assignment of vertex index, which start from 0 */
            bot[0] = dim1*(Mx + 1)+dim2;
            bot[1] = dim1*(Mx + 1)+dim2 + 1;
            top[0] = bot[0] + Mx + 1;
            top[1] = bot[1] + Mx + 1;

            /* Assignment of EToV,
             * each element start from the left lower vertex */
            EToV[dim1 * Mx + dim2][0] = bot[0];
            EToV[dim1 * Mx + dim2][1] = bot[1];
            EToV[dim1 * Mx + dim2][2] = top[1];
            EToV[dim1 * Mx + dim2][3] = top[0];
        }
    }

    /* initialize */
    dg_grid* grid = dg_grid_create(cell, K, Nv, VX, VY, NULL, EToV);

    /* free memory */
    matrix_int_free(EToV);

    return grid;
}

/**
 * @brief Generate dg_grid structure from files.
 * @details
 * The mesh files include 'casename.node' and 'casename.ele',
 * while the 'casename' is defined as parameter.
 * */
dg_grid* dg_grid_read_file2d(dg_cell *cell, char *casename){
    char element_file[MAX_NAME_LENGTH];
    char vertex_file[MAX_NAME_LENGTH];
    strcpy(element_file, casename);
    strcat(element_file, ".ele");
    strcpy(vertex_file, casename);
    strcat(vertex_file, ".node");

    /* read the element file (EToV) */
    FILE *fp;
    if( (fp = fopen(element_file, "r")) == NULL ){
        fprintf(stderr, "%s (%d)\n"
                        "Unable to open element file %s.\n",
                __FUNCTION__,__LINE__,element_file);
    }
    int Nv, K, temp;
    fscanf(fp, "%d %d %d\n", &K, &Nv, &temp);
    // check element vertex
    if(dg_cell_Nv(cell) !=  Nv){
        fprintf(stderr, "%s (%d)\n"
                "The input element type is not correct!\n", __FILE__, __LINE__);
        exit(-1);
    }

    int **EToV = matrix_int_create(K, Nv);
    int **EToR = matrix_int_create(K, 2);
    register int k,n;
    for(k=0;k<K;k++){
        fscanf(fp, "%d", EToR[k]); //read index
        EToR[k][0] -= 1; // change index start from 0 (C style)
        for(n=0;n<Nv;n++){
            fscanf(fp, "%d", EToV[0]+k*Nv+n);
            EToV[k][n] -= 1; // change index start from 0 (C style)
        }
        fscanf(fp, "%d", EToR[k]+1); //read region ncid
    }
    fclose(fp);

    /* read vertex */
    if( (fp = fopen(vertex_file, "r")) == NULL ){
        fprintf(stderr, "mr_grid_read_file (%s): %d\n"
                        "Unable to open node file %s.\n",
                __FILE__,__LINE__,vertex_file);
    }
    int Nvert;
    fscanf(fp, "%d %d %d %d\n", &Nvert, &temp, &temp, &temp);
    double *vx = vector_double_create(Nvert);
    double *vy = vector_double_create(Nvert);
    for(n=0;n<Nvert;n++){
        fscanf(fp, "%d", &temp); //read vertex ncid
        fscanf(fp, "%lf", vx+n);
        fscanf(fp, "%lf", vy+n);
    }
    fclose(fp);

    /* initialize */
    dg_grid* grid = dg_grid_create(cell, K, Nvert, vx, vy, NULL, EToV);
    dg_grid_add_EToR(grid, K, EToR);

#if DEBUG
    char filename[20] = "mr_grid_read_file2d";
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    fp = create_log(filename, procid, nprocs);
    print_int_matrix2file(fp, "EToV", EToV, K, Nv);
    print_double_vector2file(fp, "vx", vx, Nvert);
    print_double_vector2file(fp, "vy", vy, Nvert);
    print_int_matrix2file(fp, "grid->EToV", grid->EToV, grid->K, Nv);
    print_double_vector2file(fp, "vx", grid->vx, grid->Nv);
    print_double_vector2file(fp, "vy", grid->vy, grid->Nv);
    fclose(fp);
#endif
    /* free memory */
    matrix_int_free(EToV);
    matrix_int_free(EToR);
    vector_double_free(vx);
    vector_double_free(vy);

    return grid;
}


void dg_grid_read_BSfile2d(dg_grid *grid, char *casename){
    char filename[MAX_NAME_LENGTH];
    strcpy(filename, casename);
    strcat(filename, ".edge");

    FILE *fp;
    dg_fopen(fp, filename, "Unable to open boundary condition file");
    int Nsurf, tmp, n;
    fscanf(fp, "%d %d\n", &Nsurf, &tmp);
#if DEBUG
    char testname[30] = "dg_grid_read_bcfile2d_test";
    FILE *fh = CreateLog(testname, mesh->grid->procid, mesh->grid->nprocs);
    fprintf(fh, "Nsurf=%d\n", Nsurf);
#endif
    int **SFToV = NULL;
    if(Nsurf){
        SFToV = matrix_int_create(Nsurf, 3);
        for(n=0;n<Nsurf;n++){
            fscanf(fp, "%d", &tmp);
            fscanf(fp, "%d %d %d", SFToV[n], SFToV[n]+1, SFToV[n]+2);
            SFToV[n][0] -= 1;
            SFToV[n][1] -= 1; // change to C type
        }
        dg_grid_add_BS2d(grid, Nsurf, SFToV);
#if DEBUG
        PrintIntMatrix2File(fh, "SFToV", SFToV, Nsurf, 3);
        PrintIntMatrix2File(fh, "ETBS", mesh->EToBS, mesh->grid->K, mesh->cell->Nfaces);
        fclose(fh);
#endif
        matrix_int_free(SFToV);
    }else{
        dg_grid_add_BS2d(grid, Nsurf, NULL);
    }
    fclose(fp);
    return;
}

void dg_grid_add_EToR(dg_grid *grid, int cell_num, int **RgionID){
    int k;
    int *EToR = grid->EToR;
    for(k=0;k<cell_num;k++){
        int cell_id = RgionID[k][0];
        EToR[ cell_id ] = RgionID[k][1];
    }
    return;
}