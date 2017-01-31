#include <StandCell/sc_stdcell.h>
#include "mr_grid.h"

/**
 * @brief
 * Generation of uniform triangle grid.
 *
 * @details
 * Generate uniform triangle mesh with specific coordinate range and number of elements.
 * The flag 'type' defines how the square divides into two triangles.
 *
 * @param[in] shape standard cell object
 * @param[in] Mx number of elements on x coordinate
 * @param[in] My number of elements on x coordinate
 * @param[in] xmin  minimum value of x
 * @param[in] xmax  maximum value of x
 * @param[in] ymin  minimum value of y
 * @param[in] ymax  maximum value of y
 * @param[in] type  flag for triangle dividing, 0 for '\' and 1 for '/'
 *
 * @return grid geometry grid object
 */
geoGrid* mr_grid_createUniformGrid_tri
        (stdCell *shape, int Mx, int My,
         double xmin, double xmax,
         double ymin, double ymax, int type)
{
    /* check stand element */
    if(shape->type !=  TRIANGLE){
        printf("MultiRegions (mr_grid_createUniformGrid_tri): "
                       "the input element type %d is not triangle!\n", shape->type);
        exit(-1);
    }

    /* parameter */
    const int K = Mx*My*2;
    const int Nv = (Mx+1)*(My+1);

    /* allocation */
    int **EToV = matrix_int_create(K, 3);
    double VX[Nv], VY[Nv];

    /* vertex coordinate */
    int ind; /* vertex index */
    int dim1, dim2;
    for (dim1 = 0; dim1<My+1; dim1++){
        for (dim2 = 0; dim2<Mx+1; dim2++){
            ind = dim1*(Mx+1) + dim2; /* the order of vertex is from left to right, and from bottom to the top */
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
    geoGrid* grid = mr_grid_create(shape, K, Nv, VX, VY, NULL, EToV);

    /* free memory */
    matrix_int_free(EToV);
    return grid;
}

/**
 * @brief
 * Generation of uniform quadrilateral mesh.
 *
 * @details
 * Generate uniform quadrilateral mesh with specific coordinate range and number of elements.
 *
 * @param[in] shape standard cell object
 * @param[in] Mx      number of elements on x coordinate
 * @param[in] My      number of elements on x coordinate
 * @param[in] xmin    minimum value of x
 * @param[in] xmax    maximum value of x
 * @param[in] ymin    minimum value of y
 * @param[in] ymax    maximum value of y
 *
 * @return grid geometry grid object
 */
geoGrid* mr_grid_createUniformGrid_quad
        (stdCell *shape, int Mx, int My,
         double xmin, double xmax,
         double ymin, double ymax)
{
    /* check stand element */
    if(shape->type !=  QUADRIL){
        fprintf(stderr, "mr_grid_createUniformGrid_quad (%s): %d\n"
                "the input element type %d is not quadrilateral!\n",
                __FILE__, __LINE__, shape->type);
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
    geoGrid* grid = mr_grid_create(shape, K, Nv, VX, VY, NULL, EToV);

    /* free memory */
    matrix_int_free(EToV);

    return grid;
}


geoGrid* mr_grid_read_file2d(stdCell *shape, char *casename){
    char element_file[1024];
    char vertex_file[1024];
    strcpy(element_file, casename);
    strcat(element_file, ".ele");
    strcpy(vertex_file, casename);
    strcat(vertex_file, ".node");

    /* read the EToV */
    FILE *fp;
    if( (fp = fopen(element_file, "r")) == NULL ){
        fprintf(stderr, "mr_grid_read_file (%s): %d\n"
                        "Unable to open element file %s.\n",
                __FILE__,__LINE__,element_file);
    }
    int Nv, K, temp;
    fscanf(fp, "%d %d %d\n", &K, &Nv, &temp);
    // check element vertex
    if(shape->Nv !=  Nv){
        fprintf(stderr, "mr_grid_file2d (%s): %d\n"
                "the input element type is not correct!\n", __FILE__, __LINE__);
        exit(-1);
    }

    int **EToV = matrix_int_create(K, Nv);
    int k,n;
    for(k=0;k<K;k++){
        fscanf(fp, "%d", &temp); //read index
        for(n=0;n<Nv;n++){
            fscanf(fp, "%d", EToV[0]+k*Nv+n);
        }
        fscanf(fp, "%d", &temp); //read region id
    }
    fclose(fp);

    /* read vertex */
    if( (fp = fopen(vertex_file, "r")) == NULL ){
        fprintf(stderr, "mr_grid_read_file (%s): %d\n"
                        "Unable to open node file %s.\n",
                __FILE__,__LINE__,vertex_file);
    }
    fscanf(fp, "%d %d %d %d\n", &Nv, &temp, &temp, &temp);
    double *vx = vector_double_create(Nv);
    double *vy = vector_double_create(Nv);
    for(n=0;n<Nv;n++){
        fscanf(fp, "%d", &temp); //read vertex id
        fscanf(fp, "%lf", vx+n);
        fscanf(fp, "%lf", vy+n);
    }
    fclose(fp);

    /* initialize */
    geoGrid* grid = mr_grid_create(shape, K, Nv, vx, vy, NULL, EToV);

    /* free memory */
    matrix_int_free(EToV);
    vector_double_free(vx);
    vector_double_free(vy);

    return grid;
}