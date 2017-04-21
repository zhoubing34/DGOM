#include "dg_grid_BS.h"
#include "dg_grid_reader.h"

#define DEBUG 0
#if DEBUG
#include "Utility/unit_test.h"
#endif

/**
 * @brief Generation of uniform triangle grid.
 *
 * @details
 * Generate uniform triangle grid with specific coordinate range
 * and number of elements. The flag 'type' defines how the square
 * divides into two triangles.
 *
 * @param[in] cell standard cell object;
 * @param[in] Mx,My number of cells on x and y coordinate;
 * @param[in] xmin,xmax  range of x coordinate;
 * @param[in] ymin,ymax  range of y coordinate;
 * @param[in] type  flag for triangle dividing, 0 for '\' and 1 for '/';
 *
 * @return grid pointer to a new dg_grid object.
 */
dg_grid* dg_grid_create_uniform_tri(dg_cell *cell, int Mx, int My,
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
    double vx[Nv], vy[Nv];

    /* vertex coordinate */
    /* the order of vertex is from left to right, and from bottom to the top */
    int ind; /* vertex index */
    int row, col;
    for (row=0; row<My+1; row++){
        for (col=0; col<Mx+1; col++){
            ind = row*(Mx+1) + col;
            vx[ind] = (double)col/Mx * (xmax - xmin) + xmin;
            vy[ind] = (double)row/My * (ymax - ymin) + ymin;
        }
    }

    /* EToV connection */
    int top[2], bot[2];
    for (row=0; row<My; row++){
        for (col=0; col<Mx; col++){
            /* Assignment of vertex index, which start from 0 */
            bot[0] = row*(Mx + 1)+col;     /* bottom left  */
            bot[1] = row*(Mx + 1)+col + 1; /* bottom right */
            top[0] = bot[0] + (Mx + 1); /* top left  */
            top[1] = bot[1] + (Mx + 1); /* top right */

            /* Assignment of EToV */
            if ( type==1 ){ /* for '/' division */
                EToV[row * Mx * 2 + col][0] = bot[0];
                EToV[row * Mx * 2 + col][1] = top[1];
                EToV[row * Mx * 2 + col][2] = top[0];

                EToV[row * Mx * 2 + Mx + col][0] = bot[0];
                EToV[row * Mx * 2 + Mx + col][1] = bot[1];
                EToV[row * Mx * 2 + Mx + col][2] = top[1];
            }else if( type==0 ){ /* for '\' division */
                EToV[row * Mx * 2 + col][0] = bot[1];
                EToV[row * Mx * 2 + col][1] = top[1];
                EToV[row * Mx * 2 + col][2] = top[0];

                EToV[row * Mx * 2 + Mx + col][0] = bot[0];
                EToV[row * Mx * 2 + Mx + col][1] = bot[1];
                EToV[row * Mx * 2 + Mx + col][2] = top[0];
            }
        }
    }

    /* initialize */
    dg_grid* grid = dg_grid_create(cell, K, Nv, vx, vy, NULL, EToV);

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
dg_grid* dg_grid_create_uniform_quad(dg_cell *cell, int Mx, int My,
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
 * @brief create new dg_grid structure from input files.
 * @details
 * The mesh files include
 * - *casename*.node
 * - *casename*.ele
 * - *casename*.edge
 *
 * while the **casename** is an input parameter.
 * File format
 * @param cell pointer to dg_cell structure;
 * @param casename name of case;
 * */
dg_grid* dg_grid_create_from_file2d(dg_cell *cell, char *casename){

    /* read the element file */
    char element_file[MAX_NAME_LENGTH];
    strcpy(element_file, casename);
    strcat(element_file, ".ele");
    FILE *fp;
    if( (fp = fopen(element_file, "r")) == NULL ){
        fprintf(stderr, "%s (%d): Unable to open element file %s.\n", __FUNCTION__,__LINE__, element_file);
        exit(-1);
    }
    /* read element number */
    int nfld, K;
    fscanf(fp, "%d %d\n", &K, &nfld);
    /* check element vertex */
    if( (dg_cell_Nv(cell)+1) !=  nfld){
        fprintf(stderr, "%s (%d): The input element vertex number in %s is not correct!\n",
                __FUNCTION__, __LINE__, element_file);
        exit(-1);
    }
    fclose(fp);

    int **vlist = matrix_int_create(nfld, K);
    int **EToV = matrix_int_create(K, dg_cell_Nv(cell));
    switch (dg_cell_celltype(cell)){
        case TRIANGLE:
            if(read_int_file(element_file, K, nfld, vlist[0], vlist[1], vlist[2], vlist[3])){
                fprintf(stderr, "%s (%d): Format error in element file %s\n",
                        __FUNCTION__, __LINE__, element_file);
                exit(-1);
            }
            break;
        case QUADRIL:
            if(read_int_file(element_file, K, nfld, vlist[0], vlist[1], vlist[2], vlist[3], vlist[4])){
                fprintf(stderr, "%s (%d): Format error in element file %s\n",
                        __FUNCTION__, __LINE__, element_file);
                exit(-1);
            }
            break;
        default:
            fprintf(stderr, "%s (%d): Unknown cell type %d\n",
                    __FUNCTION__, __LINE__, dg_cell_celltype(cell));
            exit(-1);
    }
    /* change vertex index to C type */
    int n,k;
    for(k=0;k<K;k++){
        for(n=0;n<(nfld-1);n++){
            EToV[k][n] = vlist[n][k] - 1;
        }
    }
    ///todo: set EToR property

#if DEBUG
    printf("%s (%d): finish reading element file\n", __FUNCTION__, __LINE__);
#endif
    /* read node file */
    char vertex_file[MAX_NAME_LENGTH];
    strcpy(vertex_file, casename);
    strcat(vertex_file, ".node");
    if( (fp = fopen(vertex_file, "r")) == NULL ){
        fprintf(stderr, "%s (%d): Unable to open node file %s.\n",
                __FUNCTION__,__LINE__, vertex_file);
    }
    int Nvert;
    fscanf(fp, "%d\n", &Nvert);
    double *vx = vector_double_create(Nvert);
    double *vy = vector_double_create(Nvert);
    fclose(fp);
    /* read vertex */
    if( read_double_file(vertex_file, Nvert, 2, vx, vy) ){
        fprintf(stderr, "%s (%d): Format error in vertex file %s\n", __FUNCTION__, __LINE__, vertex_file);
    }

    /* create dg_grid */
    dg_grid* grid = dg_grid_create(cell, K, Nvert, vx, vy, NULL, EToV);
#if DEBUG
    printf("%s (%d): finish reading node file\n", __FUNCTION__, __LINE__);
#endif
    /* read boundary file */
    char edge_file[MAX_NAME_LENGTH];
    strcpy(edge_file, casename);
    strcat(edge_file, ".edge");
    if( (fp = fopen(edge_file, "r")) == NULL ){
        fprintf(stderr, "%s (%d): Unable to open edge file %s.\n",
                __FUNCTION__,__LINE__, edge_file);
    }
    int Nsurf;
    fscanf(fp, "%d\n", &Nsurf);
    int **slist = matrix_int_create(3, Nsurf);
    int **SFToV = matrix_int_create(Nsurf, 3);
    fclose(fp);
    if( read_int_file(edge_file, Nsurf, 3, slist[0], slist[1], slist[2]) ){
        fprintf(stderr, "%s (%d): Format error in edge file %s\n", __FUNCTION__, __LINE__, edge_file);
    }
    /* change vertex index to C type */
    for(k=0;k<Nsurf;k++){
        for(n=0;n<2;n++){
            SFToV[k][n] = slist[n][k] - 1;
        }
        SFToV[k][2] = slist[2][k];
    }
    if(Nsurf) { grid->set_EToBS(grid, Nsurf, SFToV); }
#if DEBUG
    printf("%s (%d): finish reading boundary file\n", __FUNCTION__, __LINE__);
#endif

#if DEBUG
    int procid, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &procid);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    fp = create_log(__FUNCTION__, procid, nprocs);
    print_int_matrix2file(fp, "EToV", EToV, K, nfld);
    print_double_vector2file(fp, "vx", vx, Nvert);
    print_double_vector2file(fp, "vy", vy, Nvert);
    print_int_matrix2file(fp, "slist", slist, 3, Nsurf);
    print_int_matrix2file(fp, "SFToV", SFToV, Nsurf, 3);
    print_int_matrix2file(fp, "grid->EToV", grid->EToV, grid->K, nfld);
    print_double_vector2file(fp, "grid->vx", grid->vx, grid->nfld);
    print_double_vector2file(fp, "grid->vy", grid->vy, grid->nfld);
    print_int_matrix2file(fp, "grid->EToBS", grid->EToBS, grid->K, nfld);
    fclose(fp);
#endif
    /* free memory */
    matrix_int_free(vlist);
    matrix_int_free(slist);
    matrix_int_free(EToV);
    matrix_int_free(SFToV);
    vector_double_free(vx);
    vector_double_free(vy);
    return grid;
}