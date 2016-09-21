#include <stdio.h>
#include "LibUtilities/LibUtilities.h"

#define MAXLEN 1024
#define HEADLEN 20  // the length of title section in *.ele and *.edge files

/* private function */
int RemoveSuffix(char *fullname, char *filename);
int WriteVertexFile(int ndim, FILE *fp, char *casename, char **argv);
int WriteEleFile2d(int ndim, FILE *fp, char *casename, char **argv);

/**
 * @brief
 * Read mesh file of Gmsh and write to the result files.
 *
 * @details
 *
 * Usages:
 *      GmshReader 2 data.msh
 *
 * The output files include:
 *  1. `data.node`, vertex coordinate;
 *  2. `data.ele`, vertex index in each element;
 *  3. `data.edge` or `data.face`, boundary element.
 */
int main(int argc, char **argv){

    if(argc!=3){ // check input
        printf("Wrong number of input arguments\n"); exit(-1);
    }

    int ndim; // get dimensions
    str2int(argv[1], &ndim, "Wrong value of dimensions.");

    // open the mesh file and check exist
    FILE *fp = fopen(argv[2], "r");
    if(fp == NULL){
        printf("Can't open Gmsh file %s\n", argv[2]); exit(-1);
    }

    char casename[MAXLEN]; // get the filename
    if(RemoveSuffix(argv[2], casename) < 0)
        exit(-1);

    // read and write the vertex
    if(WriteVertexFile(ndim, fp, casename, argv)<0)
        exit(-1);

    // read and write the elements and boundaries
    if( WriteEleFile2d(ndim, fp, casename, argv)<0)
        exit(-1);

    // close file
    fclose(fp);

    return 0;
}

/* read and write the *.ele and *.edge files */
int WriteEleFile2d(int ndim, FILE *fp, char *casename, char **argv){
    // create output file
    char fullname[MAXLEN];
    strcpy(fullname, casename);
    strncat(fullname, ".ele", strlen(".ele"));

    FILE *elefp = fopen(fullname, "w");
    if(elefp == NULL){
        printf("Can't create ele file %s\n", fullname); return(-1);
    }

    strcpy(fullname, casename);
    strncat(fullname, ".edge", strlen(".edge"));

    FILE *edgefp = fopen(fullname, "w");
    if(edgefp == NULL){
        printf("Can't create edge file %s\n", fullname); return(-1);
    }

    // read title
    int i;
    char strbuff[MAXLEN];
    for(i=0;i<2;i++)
        fgets(strbuff, MAXLEN, fp);

    for(i=0;i<HEADLEN;i++){
        fprintf(edgefp, " ");
        fprintf(elefp, " ");
    }
    fprintf(edgefp, "\n");
    fprintf(elefp, "\n");

    if( !memcpy(strbuff, "$Elements", strlen("$Elements")) ){ // check file format
        printf("File format error, start of element section, %s", strbuff); return -1;
    }

    // read EToV
    int Ne;
    fgets(strbuff, MAXLEN, fp);
    if (sscanf(strbuff, "%d", &Ne) != 1){
        printf("Wrong No. of elements %s\n", strbuff); return -1;
    }
    int id,type,physid,triConunter,quadContour,edgeContour,temp;
    int node[4];
    for(triConunter=0,quadContour=0,edgeContour=0, i=0;i<Ne;i++){
        fgets(strbuff, MAXLEN, fp);
        sscanf(strbuff, "%d %d", &id, &type);
        if(type==1){ //edge
            edgeContour++;
            sscanf(strbuff, "%d %d %d %d %d %d %d",
                   &id, &type, &temp, &physid, &temp, node, node+1);
            fprintf(edgefp, "%d %d %d %d\n", edgeContour, node[0], node[1], physid);
        }else if(type==2){ //triangle
            triConunter++;
            sscanf(strbuff, "%d %d %d %d %d %d %d %d",
                   &id, &type, &temp, &physid, &temp, node, node+1, node+2);
            fprintf(elefp, "%d %d %d %d %d\n",
                    triConunter, node[0], node[1], node[2], physid);
        }else if(type==3){ //quadrilateral
            quadContour++;
            sscanf(strbuff, "%d %d %d %d %d %d %d %d %d",
                   &id, &type, &temp, &physid, &temp, node, node+1, node+2, node+3);
            fprintf(elefp, "%d %d %d %d %d %d\n",
                    quadContour, node[0], node[1], node[2], node[3], physid);
        }else{
            printf("Wrong element type %d in Elements %d", type, id); return -1;
        }
    }

    int eleContour=0;
    type = 0;
    if( (triConunter>0) & (quadContour>0) ){
        printf("Error: Mixed triangle/quad mesh encountered\n"); return -1;
    }else if(triConunter>0){
        eleContour = triConunter;
        type = 3;
    }else if(quadContour>0){
        eleContour = quadContour;
        type = 4;
    }

    printf("Nbc = %d\n", edgeContour);
    printf("Ne = %d\n", eleContour);

    // rewind file to write total number
    rewind(edgefp);
    sprintf(strbuff, "%d 1", edgeContour);
    if(strlen(strbuff) <= HEADLEN){ // check the length
        fprintf(edgefp, "%s", strbuff);
    }else{
        printf("The pre-assigned title length %d is not long enough!\n",HEADLEN);
        return -1;
    }

    fseek(edgefp, 0, SEEK_END);
    fprintf(edgefp, "# Produce by: %s %s %s", argv[0], argv[1], argv[2]);

    rewind(elefp);
    sprintf(strbuff, "%d %d 1", eleContour, type);
    if(strlen(strbuff) <= HEADLEN){ // check the length
        fprintf(elefp, "%s", strbuff);
    }else{
        printf("The pre-assigned title length %d is not long enough!\n",HEADLEN);
        return -1;
    }
    fseek(elefp, 0, SEEK_END);
    fprintf(elefp, "# Produce by: %s %s %s", argv[0], argv[1], argv[2]);

    fclose(edgefp);
    fclose(elefp);

    return 0;
}

/* read the node coordinate and write to .node file */
int WriteVertexFile(int ndim, FILE *fp, char *casename, char **argv){
    // create output file
    char fullname[MAXLEN];
    strcpy(fullname, casename);
    strncat(fullname, ".node", strlen(".node"));

    FILE *outfp = fopen(fullname, "w");
    if(outfp == NULL){
        printf("Can't create node file %s\n", fullname); return(-1);
    }

    // read the title of mesh file, tile 4th line
    char strbuff[MAXLEN];
    int i;
    for(i=0; i<4; i++)
        fgets(strbuff, MAXLEN, fp);

    if( !memcpy(strbuff, "$Nodes", strlen("$Nodes")) ){ // check file format
        printf("File format error, 4th line, %s", strbuff); return -1;
    }

    // # of vertex
    int Nv;
    fgets(strbuff, MAXLEN, fp);
    if (sscanf(strbuff, "%d", &Nv) != 1){
        printf("Wrong No. of vertex %s\n", strbuff); return -1;
    }
    printf("Nv = %d\n", Nv); // info
    fprintf(outfp, "%d %d 0 0\n", Nv, ndim); // write to node file

    // vertex
    if(ndim == 2) {
        float coor[2];
        int k;
        for (i = 0; i < Nv; i++) {
            fgets(strbuff, MAXLEN, fp);
            sscanf(strbuff, "%d %f %f",&k,coor,coor+1);
            fprintf(outfp, "%d %f %f\n", k, coor[0], coor[1]);
        }
    }else if(ndim == 3){
        float coor[3];
        int k;
        for (i=0; i<Nv; i++) {
            fgets(strbuff, MAXLEN, fp);
            sscanf(strbuff, "%d %f %f %f",&k,coor,coor+1,coor+2);
            fprintf(outfp, "%d %f %f %f\n", k, coor[0], coor[1], coor[2]);
        }
    }else{
        printf("Wrong value of dimensions %d\n", ndim); return -1;
    }

    // finish
    fprintf(outfp, "# Produce by: %s %s %s", argv[0], argv[1], argv[2]);
    fclose(outfp);
    return 0;
}


/* remove the suffix of filename */
int RemoveSuffix(char *fullname, char *filename){
    char *last_dot = strrchr(fullname, '.');
    if (last_dot != NULL && strrchr(fullname, '\\') < last_dot) {
        *last_dot = '\0';
        strcpy(filename, fullname);
    }else{
        printf("Wrong name of mesh file %s\n", fullname);
        return -1;
    }
    return 0;
}