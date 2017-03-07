//
// Created by li12242 on 12/21/16.
//

#include "mr_mesh_bc.h"

#define DEBUG 0
#if DEBUG
#include "Utility/UTest.h"
#endif

/* count the number of uique elements in an array */
static int count_unique_integer(int len, int *list);
/* sort numbers from small to large */
static int mr_mesh_cmp(const void *a, const void *b);

/* find the boundary type id, in bcTypeList */
#define _surfTypeId(mesh, ind, typeid) do{\
typeid = ind;\
}while(0)

/**
 * @brief add boundary conditions from file.
 * @param mesh
 * @param casename
 */
void mr_mesh_read_bcfile2d(dg_mesh *mesh, char *casename){
    char filename[MAX_NAME_LENGTH];
    strcpy(filename, casename);
    strcat(filename, ".edge");

    FILE *fp;
    dg_fopen(fp, filename, "Unable to open boundary condition file");
    int Nsurf, tmp, n;
    fscanf(fp, "%d %d\n", &Nsurf, &tmp);
#if DEBUG
    char testname[30] = "mr_mesh_read_bcfile2d_test";
    FILE *fh = CreateLog(testname, mesh->grid->procid, mesh->grid->nprocs);
    fprintf(fh, "Nsurf=%d\n", Nsurf);
#endif
    int **SFToV = NULL;
    if(Nsurf>0){
        SFToV = matrix_int_create(Nsurf, 3);
        for(n=0;n<Nsurf;n++){
            fscanf(fp, "%d", &tmp);
            fscanf(fp, "%d %d %d", SFToV[n], SFToV[n]+1, SFToV[n]+2);
            SFToV[n][0] -= 1;
            SFToV[n][1] -= 1; // change to C type
        }
        mr_mesh_add_bc2d(mesh, Nsurf, SFToV);
        int Nobc = mesh->Nobc;
        mesh->obcfilename = (char**) calloc(Nobc, sizeof(char *));
        for(n=0;n<Nobc;n++){
            mesh->obcfilename[n] = (char*) calloc(MAX_NAME_LENGTH, sizeof(char));
            int ind = mesh->obcind[n];
            char indstr[20];
            sprintf(indstr, "%d", ind);
            strcpy(mesh->obcfilename[n], casename);
            strcat(mesh->obcfilename[n], ".obc");
            strcat(mesh->obcfilename[n], indstr);
            strcat(mesh->obcfilename[n], ".nc");
        }
#if DEBUG
        PrintIntMatrix2File(fh, "SFToV", SFToV, Nsurf, 3);
        PrintIntMatrix2File(fh, "ETBS", mesh->EToBS, mesh->grid->K, mesh->cell->Nfaces);
        fclose(fh);
#endif
        matrix_int_free(SFToV);
    }else{
        mr_mesh_add_bc2d(mesh, Nsurf, NULL);
    }
    fclose(fp);
    return;
}

/**
 * @brief add boundary conditions into the 2d mesh object
 *
 * @details
 * The `SFToV` matrix contains the vertex list and surface type indicator (integer) of each surface,
 * .e.g, [v1, v2, typeid].
 *
 * @param[in,out] mesh mesh object
 * @param[in] Nsurf number of surface
 * @param[in] SFToV surface to vertex list
 */
void mr_mesh_add_bc2d(dg_mesh *mesh, int Nsurf, int **SFToV){
    int k,f1,f2;

    const int K = mesh->grid->K;
    const int Nfaces = mesh->cell->Nfaces;
    const int Nv = mesh->cell->Nv; ///< number of vertex in each element
    const int Nvert = mesh->grid->Nv; ///< number of all the vertex

    int **EToV = mesh->grid->EToV;

    /* count obc number */
    int surfList[Nsurf], v[2];
    for(f1=0;f1<Nsurf;f1++){
        v[0] = SFToV[f1][0]; v[1] = SFToV[f1][1];
        qsort(v, 2, sizeof(int), mr_mesh_cmp);
        SFToV[f1][0] = v[0]; SFToV[f1][1] = v[1];

        surfList[f1] = SFToV[f1][2];
        if(SFToV[f1][2] == INNERLOC | SFToV[f1][2] == INNERBS ){
            printf("%s (%d): Error boundary type in SFToV[%d][2] = %d\n",
                   __FUNCTION__, __LINE__, f1, SFToV[f1][2]);
            printf("The boundary type indicator cannot be %d or %d:\n", INNERLOC, INNERBS);
            printf("   %d - local boundary surface[default]\n", INNERLOC);
            printf("   %d - parallel boundary surface\n", INNERBS);
            printf("please use other integers for boundary ID:\n");
            printf("   %d - slip wall\n", SLIPWALL);
            printf("   %d - non-slip wall\n", NSLIPWALL);
            printf("   other ids - different open boundaries\n");
        }
    }
    int Nbc = count_unique_integer(Nsurf, surfList);
    mesh->Nbc = Nbc;
    /* store the boundary type indicator (from smallest to largest) */
    mesh->bcind = vector_int_create(Nbc);
    mesh->bcind[0] = surfList[0];
    int sk = 1;
    for(f1=0;f1<(Nsurf-1);f1++){
        if(surfList[f1+1]-surfList[f1] != 0)
            mesh->bcind[sk++] = surfList[f1+1];
    }

    mesh->Nobc = 0;
    for(f1=0;f1<Nbc;f1++){
        if( mesh->bcind[f1]>NSLIPWALL ){ mesh->Nobc += 1; }
    }
    mesh->obcind = vector_int_create(mesh->Nobc);
    sk = 0;
    for(f1=0;f1<Nbc;f1++){
        if( mesh->bcind[f1]>NSLIPWALL ){ mesh->obcind[sk++] = mesh->bcind[f1]; }
    }

    /* now allocate the element to surface type matrix */
    mesh->EToBS = matrix_int_create(K, Nfaces);

    int t[2];
    for(k=0;k<K;k++){
        for(f2=0; f2<Nfaces; f2++){ // loop through all element surface
            /* inner surface (default) */
            mesh->EToBS[k][f2] = INNERLOC;

            /* the inner boundary surface */
            if(mesh->EToP[k][f2] != mesh->procid){
                mesh->EToBS[k][f2] = INNERBS;
                continue; // finish this loop
            }

            /* get the usr-defined boundary */
            int n1 = f2;
            int n2 = (f2+1)%Nv;
            t[0] = EToV[k][n1];
            t[1] = EToV[k][n2];
            qsort(t, 2, sizeof(int), mr_mesh_cmp);
            int t_temp = t[0]*Nvert + t[1];
#if 0
            if(mesh->procid == 0)
                printf("k=%d, f=%d, t1=%d, t2=%d, t_temp = %d, \n", k, f2, t[0], t[1], t_temp);
#endif
            for(f1=0; f1<Nsurf; f1++){
                v[0] = SFToV[f1][0];
                v[1] = SFToV[f1][1];
                int v_temp = v[0]*Nvert + v[1];
#if 0
                if(mesh->procid == 0)
                    printf("surfid=%d, v1=%d, v2=%d, v_temp = %d\n", f1, v[0], v[1], v_temp);
#endif
                if(t_temp == v_temp){ // compare the face
                    _surfTypeId(mesh, SFToV[f1][2], mesh->EToBS[k][f2]);
                    break; // jump out loop of SFToV
                }
            }
        }
    }
}

/**
 * @brief delete the boundary relative properties
 */
void mr_mesh_del_bc2d(dg_mesh *mesh){
    matrix_int_free(mesh->EToBS);
    vector_int_free(mesh->bcind);
    vector_int_free(mesh->obcind);
}

/* sort numbers from small to large */
static int mr_mesh_cmp(const void *a, const void *b){
    return (* (int *)a) - (* (int *)b);
}

/**
 * @brief count the number of uique elements in an array
 * @param len length of the array
 * @param list array of integer
 * @return n number of unique elements
 */
static int count_unique_integer(int len, int *list){

    if(len ==0) { return 0; }
    qsort(list, (size_t)len, sizeof(int), mr_mesh_cmp); // sort the list
    int n=1, i;
    for(i=0;i<(len-1);i++){
        if(list[i+1]-list[i] != 0) { n++; }
    }
    return n;
}