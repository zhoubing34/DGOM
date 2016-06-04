#include "mpi.h"
#include "Maxwell2d/fem.h"

main(int argc, char **argv){

  Mesh *mesh;
  int procid, nprocs, maxNv;
  double minEz,maxEz, gminEz, gmaxEz;

  /* initialize MPI */
  MPI_Init(&argc, &argv);

  /* assign gpu */
  MPI_Comm_rank(MPI_COMM_WORLD, &procid);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if(procid==0){
    printf("--------------------------------\n");
    printf("            MIDG\n");
    printf("--------------------------------\n");
    printf("\n     2d TM mode Maxwells\n");
    printf("\n         N = %d \n", p_N);
    printf("--------------------------------\n");
  }

#ifdef CUDA  
  cudaSetDevice((procid+2)%4);
#endif

  /* (parallel) read part of fem mesh from file */
  mesh = ReadMesh2d(argv[1]);

  /* perform load balancing */
  LoadBalance2d(mesh);

  /* find element-element connections */
  FacePair2d(mesh, &maxNv);

  /* perform start up */
  StartUp2d(mesh);

  /* field storage (double) */
  double *Hx = (double*) calloc(mesh->K*p_Np, sizeof(double));
  double *Hy = (double*) calloc(mesh->K*p_Np, sizeof(double));
  double *Ez = (double*) calloc(mesh->K*p_Np, sizeof(double));

  /* initial conditions */
  int k,n, sk=0;
  for(k=0;k<mesh->K;++k){
    for(n=0;n<p_Np;++n) {
      Hx[sk] = 0;
      Hy[sk] = 0;
      Ez[sk] = sin(M_PI*mesh->x[k][n])*sin(M_PI*mesh->y[k][n]);
      ++sk;
    }
  }

  double dt;
#ifdef CUDA
  /* initialize GPU info */
  dt = InitGPU2d(mesh, p_Nfields);

  /* load data onto GPU */
  gpu_set_data2d(mesh->K, Hx, Hy, Ez);
#else
  /* initialize GPU info */
  dt = InitCPU2d(mesh, p_Nfields);

  /* load data onto CPU float storage */
  cpu_set_data2d(mesh, Hx, Hy, Ez);
#endif

  double gdt;
  MPI_Allreduce(&dt, &gdt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

  dt = .5*gdt/((p_N+1)*(p_N+1));

  if(procid==0)
    printf("dt=%f \n", dt);
  
  double FinalTime = .1;


  /* solve */
  void MaxwellsRun2d(Mesh *mesh, double FinalTime, double dt);
  MaxwellsRun2d(mesh, FinalTime, dt); 

#ifdef CUDA
  /* unload data from GPU */
  gpu_get_data2d(mesh->K, Hx, Hy, Ez);
#else
  cpu_get_data2d(mesh, Hx, Hy, Ez);
#endif

  /* find maximum & minimum values for Ez */
  minEz=Ez[0], maxEz=Ez[0];

  for(n=0;n<mesh->K*p_Np;++n) {
    minEz = (minEz>Ez[n])?Ez[n]:minEz;
    maxEz = (maxEz<Ez[n])?Ez[n]:maxEz;
  }

  MPI_Reduce(&minEz, &gminEz, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&maxEz, &gmaxEz, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(procid==0)
    printf("t=%f Ez in [ %f, %f ] \n", FinalTime, gminEz, gmaxEz );

  /* nicely stop MPI */
  MPI_Finalize();

  /* end game */
  exit(0);
}
