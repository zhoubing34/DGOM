/* -*- mode: C; c-basic-offset: 8; c-indent-level: 8; c-continued-statement-offset: 8; c-label-offset: -8; -*- */

#include <stdio.h>
#include <cuda.h>

texture<float4, 1, cudaReadModeElementType> t_LIFT;
texture<float4, 1, cudaReadModeElementType> t_DrDs;
texture<float, 1, cudaReadModeElementType> t_vgeo;
texture<float, 1, cudaReadModeElementType> t_Q;
texture<float, 1, cudaReadModeElementType> t_partQ;
texture<float, 1, cudaReadModeElementType> t_surfinfo;

static float *c_LIFT;
static float4 *c_DrDs;
static float2 *c_DrDs2;
static float *c_surfinfo;
static float *c_vgeo;
static float *c_Q; 
static float *c_partQ; 
static float *c_rhsQ; 
static float *c_resQ; 
static float *c_tmp;

extern "C"
{

#include "fem.h"

double InitGPU2d(Mesh *mesh, int Nfields){

  printf("Np = %d, BSIZE = %d\n", p_Np, BSIZE);

  /* Q  */
  int sz = mesh->K*(BSIZE)*Nfields*sizeof(float); 

  float *f_Q = (float*) calloc(mesh->K*BSIZE*Nfields, sizeof(float));
  cudaMalloc  ((void**) &c_Q, sz);
  cudaMalloc  ((void**) &c_rhsQ, sz);
  cudaMalloc  ((void**) &c_resQ, sz);
  cudaMalloc  ((void**) &c_tmp, sz);
  cudaMemcpy( c_Q,    f_Q, sz, cudaMemcpyHostToDevice);
  cudaMemcpy( c_rhsQ, f_Q, sz, cudaMemcpyHostToDevice);
  cudaMemcpy( c_resQ, f_Q, sz, cudaMemcpyHostToDevice);
  cudaMemcpy( c_tmp,  f_Q, sz, cudaMemcpyHostToDevice);
  
  cudaBindTexture(0,  t_Q, c_Q, sz); 
  
  sz = mesh->parNtotalout*sizeof(float);
  cudaMalloc((void**) &c_partQ, sz);
  cudaBindTexture(0,  t_partQ, c_partQ, sz); 

  /*  LIFT  */
  sz = p_Np*(p_Nfp)*(p_Nfaces+1)*sizeof(float);
#if 0
   float *f_LIFT = (float*) malloc(sz);
   int skL = 0;

   for(int m=0;m<p_Nfp*p_Nfaces;++m){
     for(int n=0;n<p_Np;++n){
       f_LIFT[skL++] = mesh->LIFT[n][m];
     }
   }
#else
   float *f_LIFT = (float*) malloc(sz);
   int skL = 0;
   for(int m=0;m<p_Nfp;++m){
     for(int n=0;n<p_Np;++n){
       for(int f=0;f<p_Nfaces;++f){
	 f_LIFT[skL++] = mesh->LIFT[0][p_Nfp*p_Nfaces*n+(f+p_Nfaces*m)];
       }
       ++skL;
     }
   }
#endif
   cudaMalloc  ((void**) &c_LIFT, sz);
   cudaMemcpy( c_LIFT, f_LIFT, sz, cudaMemcpyHostToDevice);

#if 1
   /* Bind the array to the texture */
   cudaBindTexture(0,  t_LIFT, c_LIFT, sz);

   /* DrDsDt */
   sz = BSIZE*BSIZE*4*sizeof(float);

   float* h_DrDs = (float*) calloc(BSIZE*BSIZE, sizeof(float4));
   int sk = 0;
   /* note transposed arrays to avoid "bank conflicts" */
   for(int n=0;n<p_Np;++n){ 
     for(int m=0;m<p_Np;++m){
       h_DrDs[4*(m+n*BSIZE)+0] = mesh->Dr[0][n+m*p_Np];
       h_DrDs[4*(m+n*BSIZE)+1] = mesh->Ds[0][n+m*p_Np];
#if (p_Np%2)==0
       h_DrDs[4*(m+n*BSIZE)+2] = mesh->Dr[0][n+1+m*p_Np];
       h_DrDs[4*(m+n*BSIZE)+3] = mesh->Ds[0][n+1+m*p_Np];
#endif

     }
   }
	   
   cudaMalloc  ((void**) &c_DrDs, sz);
   cudaMemcpy( c_DrDs, h_DrDs, sz, cudaMemcpyHostToDevice);

   /* Bind the array to the texture */
   cudaBindTexture(0,  t_DrDs, c_DrDs, sz);

   sz = BSIZE*BSIZE*2*sizeof(float);
   float* h_DrDs2 = (float*) calloc(BSIZE*BSIZE, sizeof(float2));
   sk = 0;
   /* note transposed arrays to avoid "bank conflicts" */
   for(int n=0;n<p_Np;++n){ 
     for(int m=0;m<p_Np;++m){
       h_DrDs2[2*(m+n*BSIZE)+0] = mesh->Dr[0][n+m*p_Np]; 
       h_DrDs2[2*(m+n*BSIZE)+1] = mesh->Ds[0][n+m*p_Np]; 
     }
   }
   cudaMalloc  ((void**) &c_DrDs2, sz);
   cudaMemcpy( c_DrDs2, h_DrDs2, sz, cudaMemcpyHostToDevice);
   
   free(h_DrDs);

   /* vgeo */
   double drdx, dsdx, drdy, dsdy, J;
   float *vgeo = (float*) calloc(4*mesh->K, sizeof(float));

   for(int k=0;k<mesh->K;++k){
     GeometricFactors2d(mesh, k, &drdx, &dsdx, &drdy, &dsdy, &J);
     vgeo[k*4+0] = drdx;
     vgeo[k*4+1] = drdy;
     vgeo[k*4+2] = dsdx;
     vgeo[k*4+3] = dsdy;
   }

   sz = mesh->K*4*sizeof(float);
   cudaMalloc  ((void**) &c_vgeo, sz);
   cudaMemcpy( c_vgeo, vgeo, sz, cudaMemcpyHostToDevice);
   cudaBindTexture(0,  t_vgeo, c_vgeo, sz);
   
   /* surfinfo (vmapM, vmapP, Fscale, Bscale, nx, ny, nz, 0) */
   sz = mesh->K*p_Nfp*p_Nfaces*6*sizeof(float); 
   float* h_surfinfo = (float*) malloc(sz); 
   
   /* local-local info */
   sk = 0;
   int skP = -1;
   double *nxk = BuildVector(mesh->Nfaces);
   double *nyk = BuildVector(mesh->Nfaces);
   double *sJk = BuildVector(mesh->Nfaces);

   double dt = 1e6;

   for(int k=0;k<mesh->K;++k){
     GeometricFactors2d(mesh, k, &drdx, &dsdx, &drdy, &dsdy, &J);     
     Normals2d(mesh, k, nxk, nyk, sJk);
     
     for(int f=0;f<mesh->Nfaces;++f){
       dt = min(dt, J/sJk[f]);
  
       for(int m=0;m<p_Nfp;++m){
	 int n = m + f*p_Nfp + p_Nfp*p_Nfaces*k;
	 int idM = mesh->vmapM[n];
	 int idP = mesh->vmapP[n];
	 int  nM = idM%p_Np; 
	 int  nP = idP%p_Np; 
	 int  kM = (idM-nM)/p_Np;
	 int  kP = (idP-nP)/p_Np;
	 idM = nM + Nfields*BSIZE*kM;
	 idP = nP + Nfields*BSIZE*kP;
	 
	 /* stub resolve some other way */
	 if(mesh->vmapP[n]<0){
	   idP = mesh->vmapP[n]; /* -ve numbers */
	 }
 
	 sk = 6*p_Nfp*p_Nfaces*k+m+f*p_Nfp;
	 h_surfinfo[sk + 0*p_Nfp*p_Nfaces] = idM;
	 h_surfinfo[sk + 1*p_Nfp*p_Nfaces] = idP;
	 h_surfinfo[sk + 2*p_Nfp*p_Nfaces] = sJk[f]/(2.*J);
	 h_surfinfo[sk + 3*p_Nfp*p_Nfaces] = (idM==idP)?-1.:1.; 
	 h_surfinfo[sk + 4*p_Nfp*p_Nfaces] = nxk[f];
	 h_surfinfo[sk + 5*p_Nfp*p_Nfaces] = nyk[f];
       }
     }
  }
   
   cudaMalloc  ((void**) &c_surfinfo, sz);
   cudaMemcpy( c_surfinfo, h_surfinfo, sz, cudaMemcpyHostToDevice);

   cudaBindTexture(0,  t_surfinfo, c_surfinfo, sz);

   free(h_surfinfo);

   sz = mesh->parNtotalout*sizeof(int);
   cudaMalloc((void**) &(mesh->c_parmapOUT), sz);
   cudaMemcpy(mesh->c_parmapOUT,  mesh->parmapOUT, sz, cudaMemcpyHostToDevice);

   return dt;
#endif
}



__global__ void MaxwellsGPU_VOL_Kernel2D(float *g_rhsQ, float2 *g_DrDs){

  /* fastest */
  __device__ __shared__ float s_Q[p_Nfields*BSIZE];
  __device__ __shared__ float s_facs[4];

  /* LOCKED IN to using Np threads per block */
  const int n = threadIdx.x;
  const int k = blockIdx.x;
  
  /* "coalesced"  */
  int m = n+k*p_Nfields*BSIZE;
  int id = n;
  s_Q[id] = tex1Dfetch(t_Q, m); m+=BSIZE; id+=BSIZE;
  s_Q[id] = tex1Dfetch(t_Q, m); m+=BSIZE; id+=BSIZE;
  s_Q[id] = tex1Dfetch(t_Q, m); 

#if 1
  if(p_Np<4 && n==0)
    for(m=0;m<4;++m)
      s_facs[m] = tex1Dfetch(t_vgeo, 4*k+m);
  else if((n<4) && (p_Np>=4))
    s_facs[n] = tex1Dfetch(t_vgeo, 4*k+n);
#else
  if(n==0)
    for(m=0;m<4;++m)
      s_facs[m] = tex1Dfetch(t_vgeo, 4*k+m);
#endif
  __syncthreads();

  float dHxdr=0,dHxds=0;
  float dHydr=0,dHyds=0;
  float dEzdr=0,dEzds=0;

  float Q;
  for(m=0;p_Np-m;){
    float4 D = tex1Dfetch(t_DrDs, n+m*BSIZE);

    id = m;
    Q = s_Q[id]; dHxdr += D.x*Q; dHxds += D.y*Q;  id += BSIZE;
    Q = s_Q[id]; dHydr += D.x*Q; dHyds += D.y*Q;  id += BSIZE;
    Q = s_Q[id]; dEzdr += D.x*Q; dEzds += D.y*Q;  
    ++m;

#if (p_Np%2) == 0
    id = m;
    Q = s_Q[id]; dHxdr += D.z*Q; dHxds += D.w*Q;  id += BSIZE;
    Q = s_Q[id]; dHydr += D.z*Q; dHyds += D.w*Q;  id += BSIZE;
    Q = s_Q[id]; dEzdr += D.z*Q; dEzds += D.w*Q;  
    ++m;
#endif
  }
  
  const float drdx= s_facs[0];
  const float drdy= s_facs[1];
  const float dsdx= s_facs[2];
  const float dsdy= s_facs[3];

  m = n+p_Nfields*BSIZE*k;
  if(n<BSIZE){
  g_rhsQ[m] = -(drdy*dEzdr+dsdy*dEzds); m += BSIZE;
  g_rhsQ[m] =  (drdx*dEzdr+dsdx*dEzds); m += BSIZE;
  g_rhsQ[m] =  (drdx*dHydr+dsdx*dHyds - drdy*dHxdr-dsdy*dHxds); 
  }
}

__global__ void MaxwellsGPU_SURF_Kernel2D(float *g_Q, float *g_rhsQ){

  __device__ __shared__ float s_fluxQ[p_Nfields*p_Nfp*p_Nfaces];

  /* LOCKED IN to using Np threads per block */
  const int n = threadIdx.x;
  const int k = blockIdx.x;
  int m;

  /* grab surface nodes and store flux in shared memory */
  if(n< (p_Nfp*p_Nfaces) ){
    /* coalesced reads (maybe) */
    m = 6*(k*p_Nfp*p_Nfaces)+n;
    const  int idM   = tex1Dfetch(t_surfinfo, m); m += p_Nfp*p_Nfaces;
           int idP   = tex1Dfetch(t_surfinfo, m); m += p_Nfp*p_Nfaces;
    const  float Fsc = tex1Dfetch(t_surfinfo, m); m += p_Nfp*p_Nfaces;
    const  float Bsc = tex1Dfetch(t_surfinfo, m); m += p_Nfp*p_Nfaces;
    const  float nx  = tex1Dfetch(t_surfinfo, m); m += p_Nfp*p_Nfaces;
    const  float ny  = tex1Dfetch(t_surfinfo, m); 

    /* check if idP<0  */
    float dHx=0, dHy=0, dEz=0;

    if(idP<0){
      idP = p_Nfields*(-1-idP);
      
      dHx = Fsc*(tex1Dfetch(t_partQ, idP+0) - tex1Dfetch(t_Q, idM+0*BSIZE));
      dHy = Fsc*(tex1Dfetch(t_partQ, idP+1) - tex1Dfetch(t_Q, idM+1*BSIZE));
      dEz = Fsc*(tex1Dfetch(t_partQ, idP+2) - tex1Dfetch(t_Q, idM+2*BSIZE));
    }
    else{
      dHx = Fsc*(    tex1Dfetch(t_Q, idP+0*BSIZE) - tex1Dfetch(t_Q, idM+0*BSIZE));
      dHy = Fsc*(    tex1Dfetch(t_Q, idP+1*BSIZE) - tex1Dfetch(t_Q, idM+1*BSIZE));
      dEz = Fsc*(Bsc*tex1Dfetch(t_Q, idP+2*BSIZE) - tex1Dfetch(t_Q, idM+2*BSIZE));
    }

    const float ndotdH = nx*dHx + ny*dHy;

    m = n;
    s_fluxQ[m] = -ny*dEz + dHx - ndotdH*nx; m += p_Nfp*p_Nfaces;
    s_fluxQ[m] =  nx*dEz + dHy - ndotdH*ny; m += p_Nfp*p_Nfaces;
    s_fluxQ[m] =  nx*dHy - ny*dHx + dEz;
  }

  /* make sure all element data points are cached */
  __syncthreads();

  if(n< (p_Np))
  {
    float rhsHx = 0, rhsHy = 0, rhsEz = 0;
    
    int sk = n;
    /* can manually unroll to 4 because there are 3 faces */
    for(m=0;p_Nfaces*p_Nfp-m;){
#if 0
      float4 L;
      L.x = tex1Dfetch(t_LIFT, n+    m*p_Np); 
      L.y = tex1Dfetch(t_LIFT, n+(m+1)*p_Np); 
      L.z = tex1Dfetch(t_LIFT, n+(m+2)*p_Np); 
#else
      float4 L = tex1Dfetch(t_LIFT, sk); sk+=p_Np;
#endif
      /* broadcast */
      int sk1 = m;
      rhsHx += L.x*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsHy += L.x*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsEz += L.x*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      ++m;

      /* broadcast */
      sk1 = m;
      rhsHx += L.y*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsHy += L.y*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsEz += L.y*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      ++m;

      /* broadcast */
      sk1 = m;
      rhsHx += L.z*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsHy += L.z*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      rhsEz += L.z*s_fluxQ[sk1]; sk1 += p_Nfp*p_Nfaces;
      ++m;

    }
    
    m = n+k*p_Nfields*BSIZE;
    g_rhsQ[m] += rhsHx; m += BSIZE;
    g_rhsQ[m] += rhsHy; m += BSIZE;
    g_rhsQ[m] += rhsEz; 

  }
}


__global__ void MaxwellsGPU_RK_Kernel2D(int Ntotal, float *g_resQ, float *g_rhsQ, float *g_Q, float fa, float fb, float fdt){
  
  int n = blockIdx.x * blockDim.x + threadIdx.x;
    
  if(n<Ntotal){
    float rhs = g_rhsQ[n];
    float res = g_resQ[n];
    res = fa*res + fdt*rhs;
    
    g_resQ[n] = res;
    g_Q[n]    += fb*res;
  }

} 


/* assumes data resides on device */
void MaxwellsKernel2d(Mesh *mesh, float frka, float frkb, float fdt){

  /* grab data from device and initiate sends */
  MaxwellsMPISend2d(mesh);

  int ThreadsPerBlock, BlocksPerGrid;	

  BlocksPerGrid   = mesh->K; 
  ThreadsPerBlock = p_Np; 
  
  /* evaluate volume derivatives */
  MaxwellsGPU_VOL_Kernel2D <<< BlocksPerGrid, ThreadsPerBlock >>>  
    (c_rhsQ, c_DrDs2);

  /* finalize sends and recvs, and transfer to device */
  MaxwellsMPIRecv2d(mesh, c_partQ);

  BlocksPerGrid = mesh->K;

  if( ( p_Nfp*p_Nfaces ) > (p_Np) )
    ThreadsPerBlock = p_Nfp*p_Nfaces;
  else
    ThreadsPerBlock = p_Np;

  /* evaluate surface contributions */
  MaxwellsGPU_SURF_Kernel2D <<< BlocksPerGrid, ThreadsPerBlock >>>
    (c_Q, c_rhsQ);

  int Ntotal = mesh->K*BSIZE*p_Nfields;
  
  ThreadsPerBlock = 256;
  BlocksPerGrid = (Ntotal+ThreadsPerBlock-1)/ThreadsPerBlock;

  /* update RK Step */
  MaxwellsGPU_RK_Kernel2D<<< BlocksPerGrid, ThreadsPerBlock >>> 
    (Ntotal, c_resQ, c_rhsQ, c_Q, frka, frkb, fdt);

}




void gpu_set_data2d(int K,
		  double *d_Hx, double *d_Hy, double *d_Ez){


  float *f_Q = (float*) calloc(K*p_Nfields*BSIZE,sizeof(float));
  
  /* also load into usual data matrices */
  
  for(int k=0;k<K;++k){
    for(int n=0;n<p_Np;++n)
       f_Q[n       +k*BSIZE*p_Nfields] = d_Hx[n+k*p_Np];
    for(int n=0;n<p_Np;++n)
      f_Q[n  +BSIZE+k*BSIZE*p_Nfields] = d_Hy[n+k*p_Np];
    for(int n=0;n<p_Np;++n)
      f_Q[n+2*BSIZE+k*BSIZE*p_Nfields] = d_Ez[n+k*p_Np];
  }
  
  cudaMemcpy(c_Q, f_Q, BSIZE*K*p_Nfields*sizeof(float), cudaMemcpyHostToDevice);
  
  free(f_Q);
}
  
void gpu_get_data2d(int K,
		  double *d_Hx, double *d_Hy, double *d_Ez){

  float *f_Q = (float*) calloc(K*p_Nfields*BSIZE,sizeof(float));
  
  cudaMemcpy(f_Q, c_Q, K*BSIZE*p_Nfields*sizeof(float), cudaMemcpyDeviceToHost);

  /* also load into usual data matrices */
  
  for(int k=0;k<K;++k){
    for(int n=0;n<p_Np;++n)
      d_Hx[n+k*p_Np] = f_Q[n        +k*BSIZE*p_Nfields];
    for(int n=0;n<p_Np;++n) 
      d_Hy[n+k*p_Np] = f_Q[n  +BSIZE+k*BSIZE*p_Nfields];
    for(int n=0;n<p_Np;++n)
      d_Ez[n+k*p_Np] = f_Q[n+2*BSIZE+k*BSIZE*p_Nfields];

  }

  free(f_Q);
}

__global__ void partial_get_kernel(int Ntotal, int *g_index, float *g_partQ){
  
  int n = blockIdx.x * blockDim.x + threadIdx.x;
    
  if(n<Ntotal)
    g_partQ[n] = tex1Dfetch(t_Q, g_index[n]);
  
} 

void get_partial_gpu_data2d(int Ntotal, int *g_index, float *h_partQ){

  int ThreadsPerBlock = 256;
  int BlocksPerGrid = (Ntotal+ThreadsPerBlock-1)/ThreadsPerBlock;

  partial_get_kernel <<< BlocksPerGrid, ThreadsPerBlock >>> (Ntotal, g_index, c_tmp);

  cudaMemcpy(h_partQ, c_tmp, Ntotal*sizeof(float), cudaMemcpyDeviceToHost);
}


}
