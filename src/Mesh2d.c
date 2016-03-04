#include "mpi.h"
#include "fem.h"

/* vertices on each edge */
int vnum[3][2] = { { 0, 1}, {1, 2}, {2, 0}};

Mesh *ReadMesh2d(char *filename){

  int n;

  Mesh *mesh = (Mesh*) calloc(1, sizeof(Mesh));

  char buf[BUFSIZ];
  
  FILE *fp = fopen(filename, "r");

  /* assume modified Gambit neutral format */
  for(n=0;n<6;++n)
    fgets(buf, BUFSIZ, fp);

  fgets(buf, BUFSIZ, fp);
  sscanf(buf, "%d %d \n", &(mesh->Nv), &(mesh->K));
  mesh->Nverts = 3; /* assume triangles */
  mesh->Nfaces = 3; /* assume triangles */

  fgets(buf, BUFSIZ, fp);
  fgets(buf, BUFSIZ, fp);

  /* read vertex coordinates */
  double *VX = BuildVector(mesh->Nv);
  double *VY = BuildVector(mesh->Nv);
  for(n=0;n<mesh->Nv;++n){
    fgets(buf, BUFSIZ, fp);
    sscanf(buf, "%*d %lf %lf", VX+n, VY+n);
  }

  /* decide on parition */
  int procid, nprocs;

  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  mesh->procid = procid;
  mesh->nprocs = nprocs;

  /* assume this proc owns a block of elements */

  int Klocal, Kstart;
  int *Kprocs = (int*) calloc(nprocs, sizeof(int));
  int p;
  
  int **newEToV, *newKprocs;
  double **newVX, **newVY;
  
  Klocal = (int) ( (double)(mesh->K)/(double)nprocs );
  
  for(p=0;p<nprocs-1;++p){
    Kprocs[p] = Klocal;
  }
  Kprocs[p] = Klocal + mesh->K - nprocs*Klocal;
  
  
  Kstart= 0;
  for(p=0;p<procid;++p)
    Kstart += Kprocs[p];
  
  Klocal = Kprocs[procid];

  /* read element to vertex connectivity */
  fgets(buf, BUFSIZ, fp);
  fgets(buf, BUFSIZ, fp);
  mesh->EToV = BuildIntMatrix(Klocal, mesh->Nverts);
  mesh->GX = BuildMatrix(Klocal, mesh->Nverts);
  mesh->GY = BuildMatrix(Klocal, mesh->Nverts);

  int sk = 0;
  for(n=0;n<mesh->K;++n){
    fgets(buf, BUFSIZ, fp);
    if(n>=Kstart && n<Kstart+Klocal){
      sscanf(buf, "%*d %*d %*d %d %d %d", 
	     mesh->EToV[sk]+0, mesh->EToV[sk]+1, mesh->EToV[sk]+2);
      
      /* correct to 0-index */
      --(mesh->EToV[sk][0]);
      --(mesh->EToV[sk][1]);
      --(mesh->EToV[sk][2]);
      
      mesh->GX[sk][0] = VX[mesh->EToV[sk][0]];
      mesh->GX[sk][1] = VX[mesh->EToV[sk][1]];
      mesh->GX[sk][2] = VX[mesh->EToV[sk][2]];
      
      mesh->GY[sk][0] = VY[mesh->EToV[sk][0]];
      mesh->GY[sk][1] = VY[mesh->EToV[sk][1]];
      mesh->GY[sk][2] = VY[mesh->EToV[sk][2]];
      ++sk;
    }
  }
  fgets(buf, BUFSIZ, fp);
  fgets(buf, BUFSIZ, fp);

  mesh->K = Klocal;

  fclose(fp);

  return mesh;
  
}

void PrintMesh2d ( Mesh *mesh ){
  int n;
  printf("Mesh data: \n");
  printf("\n K = %d\n", mesh->K);
  printf("\n Nv = %d\n", mesh->Nv);
  printf("\n Nverts = %d\n", mesh->Nverts);
  printf("\n Node coordinates = \n");
  for(n=0;n<mesh->Nv;++n){
    printf("%d %g %g\n", n, mesh->GX[n], mesh->GY[n]);
  }
  printf("\n Element to vertex connectivity = \n");
  for(n=0;n<mesh->K;++n){
    printf("%d: %d %d %d \n", n, 
	   mesh->EToV[n][0], mesh->EToV[n][1], mesh->EToV[n][2]);
  }

}

void GeometricFactors2d(Mesh *mesh, int k,
		      double *drdx, double *dsdx, double *drdy, double *dsdy, 
		      double *J){

  double x1 = mesh->GX[k][0], y1 =  mesh->GY[k][0];
  double x2 = mesh->GX[k][1], y2 =  mesh->GY[k][1];
  double x3 = mesh->GX[k][2], y3 =  mesh->GY[k][2];
  
  /* compute geometric factors of the following afine map */
  /* x = 0.5*( -(r+s)*x1 + (1+r)*x2 + (1+s)*x3 ) */
  /* y = 0.5*( -(r+s)*y1 + (1+r)*y2 + (1+s)*y3 ) */

  double dxdr = (x2-x1)/2,  dxds = (x3-x1)/2;
  double dydr = (y2-y1)/2,  dyds = (y3-y1)/2;

  /* Jacobian of coordinate mapping */
  *J = -dxds*dydr + dxdr*dyds;
  
  if(*J<0)
    printf("warning: J = %lg\n", *J);
  
  /* inverted Jacobian matrix for coordinate mapping */
  *drdx =  dyds/(*J);
  *dsdx = -dydr/(*J);
  *drdy = -dxds/(*J);
  *dsdy =  dxdr/(*J);
}

void Normals2d(Mesh *mesh, int k, double *nx, double *ny, double *sJ){

  int f;

  double x1 = mesh->GX[k][0], y1 = mesh->GY[k][0];
  double x2 = mesh->GX[k][1], y2 = mesh->GY[k][1];
  double x3 = mesh->GX[k][2], y3 = mesh->GY[k][2];

  nx[0] =  (y2-y1);  ny[0] = -(x2-x1);
  nx[1] =  (y3-y2);  ny[1] = -(x3-x2);
  nx[2] =  (y1-y3);  ny[2] = -(x1-x3);

  for(f=0;f<3;++f){
    sJ[f] = sqrt(nx[f]*nx[f]+ny[f]*ny[f]);
    nx[f] /= sJ[f];
    ny[f] /= sJ[f];
    sJ[f] /= 2.;
  }
}
