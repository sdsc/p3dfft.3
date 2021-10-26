
/*
This program exemplifies using P3DFFT++ library for taking a spectral derivative in a given dimension. 

This program initializes a 3D array with a 3D sine wave, then
performs 3D forward Fourier transform, takes derivative in spectral space (multiplies by wavenumbers), then backward transform,
and checks that the results are what we expect (sine in one dimension becomes cosine)

The program expects 'stdin' file in the working directory, with a single 
line of numbers : Nx,Ny,Nz,ndim,Nrep,idir. 
Here 
  Nx,Ny,Nz are 3D grid dimensions
  ndim is the number of dimensions of processor grid (2 most likely)
  Nrep is the # of repetitions of timing loop. 
  idir is the dimension of derivative (valid values are 1 through 3, and the logical dimension is specified, i.e. actual storage dimension may be different as specified by mem-order mapping).

Optionally, a file named 'dims' can also be provided to guide in the choice
of processor geometry in case of 2D decomposition. It should contain
two numbers in a line, with their product equal to the total number
of tasks. Otherwise processor grid geometry is chosen automatically.
For better performance, experiment with this setting, varying
iproc and jproc. In many cases, minimizing iproc gives best results.
Setting it to 1 corresponds to one-dimensional decomposition.
*/

#include "p3dfft.h"
#include "compiler_check.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))

void init_wave(double *,int[3],int *,int[3]);
void print_res(double *,int *,int *,int *);
void normalize(double *,size_t,int *);
double check_res(double*,int[3],int[3],int[3],int);
void  compute_deriv(double *,double *,int[3],int[3],int[3],int[3],int);

main(int argc,char **argv)
{
  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3],gdims2[3];
  int dmap1[] = {0,1,2};// Mapping data dimension X onto processor dimension X, and so on - 
  // Initialize final grid configuration
                     // Y onto Z and Z onto X
                     // this is a Z-pencil, since Px =1 - or at least one way to define it 
                     // (the other would be (2,1,0))
  int dmap2[] = {1,2,0}; // Mapping data dimension X onto processor dimension Y, 
  int Pgrid;
  int mem_order1[3];
  // Set up memory order for the final grid layout (for complex array in Fourier space). It is more convenient to have the storage order of the array reversed, this helps save on memory access bandwidth, and shouldn't affect the operations in the Fourier space very much, requiring basically a change in the loop order. However, note that as an alternative, it is possible to define the memory ordering the same as default (0,1,2). Note that the memory ordering is specified in C indices, i.e. starting from 0.
  int mem_order2[] = {1,2,0};
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  int sdims1[3],sdims2[3];
  size_t size1,size2;
  double *IN;
  Grid *Xpencil,*Zpencil;
  int glob_start1[3],glob_start2[3],glob2[3];
  double *OUT;
  int type_ids1[3];
  int type_ids2[3];
  Type3D type_rcc,type_ccr;
  double t=0.;
  double *FIN;
  double mydiff;
  double diff = 0.0;
  double gtavg=0.;
  double gtmin=INFINITY;
  double gtmax = 0.;
  int pdims[3],nx,ny,nz,n,ndim,idir;
  Plan3D trans_f,trans_b;
  FILE *fp;
  size_t workspace_host,workspace_dev;
  int nslices=1;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);


   if(myid == 0) {
     printf("P3DFFT++ C test program. Running on %d cores\n",nprocs);
     printf("GitVersion = %s\n", GIT_VERSION);
     printf("GitDate = %s\n", GIT_DATE);
     printf("Executable, %s, was compiled with %s (version %d) on %s at %s\n", __FILE__, COMPILER_DETECTED, COMPILER_V_DETECTED, __DATE__, __TIME__);
     if((fp=fopen("stdin", "r"))==NULL){
        printf("Cannot open file. Setting to default nx=ny=nz=128, ndim=2, n=1, idir=1.\n");
        nx=ny=nz=128; Nrep=1;ndim=2;idir=1;
     } else {
       fscanf(fp,"%d %d %d %d %d %d\n",&nx,&ny,&nz,&ndim,&idir,&Nrep,&nslices);
        fclose(fp);
     }
     printf("P3DFFT spectral derivative test, 3D wave input\n");
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n idir=%d\n%d slices\n",nx,ny,nz,ndim,Nrep,idir,nslices);
#else
     printf("Single precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n idir=%d%d slices\n",nx,ny,nz,ndim,Nrep,idir,nslices);
#endif
   }
   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Nrep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&idir,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nslices,1,MPI_INT,0,MPI_COMM_WORLD);
  // Establish 2D processor grid decomposition, either by reading from file 'dims' or by an MPI default

   if(ndim == 1) {
     pdims[0] = 1; pdims[1] = nprocs;
     if(myid == 0)
       printf("Using one-dimensional decomposition\n");
   }
   else if(ndim == 2) {
     if(myid == 0)
       printf("Using two-dimensional decomposition\n");
     fp = fopen("dims","r");
     if(fp != NULL) {
       if(myid == 0)
         printf("Reading proc. grid from file dims\n");
       fscanf(fp,"%d %d\n",pdims,pdims+1);
       fclose(fp);
       if(pdims[0]*pdims[1] != nprocs)
          pdims[1] = nprocs / pdims[0];
     }
     else {
       if(myid == 0)
          printf("Creating proc. grid with mpi_dims_create\n");
       pdims[0]=pdims[1]=0;
       MPI_Dims_create(nprocs,2,pdims);
       if(pdims[0] > pdims[1]) {
          pdims[0] = pdims[1];
          pdims[1] = nprocs/pdims[0];
       }
     }
   }

   if(myid == 0)
      printf("Using processor grid %d x %d\n",pdims[0],pdims[1]);

  // Set up work structures for P3DFFT

   p3dfft_setup(nslices); // Use 8 streams/slices

  //Set up 2 transform types for 3D transforms

  type_ids1[0] = P3DFFT_R2CFFT_D; // X dimension
  type_ids1[1] = P3DFFT_CFFT_FORWARD_D; // Y dimension
  type_ids1[2] = P3DFFT_CFFT_FORWARD_D; // Z dimension

  type_ids2[0] = P3DFFT_C2RFFT_D; // X dimension
  type_ids2[1] = P3DFFT_CFFT_BACKWARD_D; // Y dimension
  type_ids2[2] = P3DFFT_CFFT_BACKWARD_D; // Z dimension

  //Now initialize 3D transform types (forward and backward) 

  type_rcc = p3dfft_init_3Dtype(type_ids1); // Forward R2C type
  type_ccr = p3dfft_init_3Dtype(type_ids2); // Backward C2R type

  //Set up global dimensions of the grid

  gdims[0] = nx;
  gdims[1] = ny;
  gdims[2] = nz;
  for(i=0; i < 3;i++) {
     mem_order1[i] = i; // The simplest case of sequential ordering
  }

  p1 = pdims[0];
  p2 = pdims[1];

  // Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively

  pdims[2] = pdims[1];
  pdims[1] = pdims[0];
  pdims[0] = 1;

  Pgrid = p3dfft_init_proc_grid(pdims,MPI_COMM_WORLD);

  // Initialize the initial grid 
                     // this is an X pencil, since Px =1

  Xpencil = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap1,mem_order1);

  //Define the final processor grid. It can be the same as initial grid, however here it is different. 
  // The final grid may have different processor distribution and memory ordering, since for example many applications with convolution and those solving partial differential equations do not need the initial grid configuration in Fourier space. The flow of these applications is typically 1) transform from physical to Fourier space, 2) apply convolution or derivative calculation in Fourier space, and 3) inverse FFT to physical space. Since forward FFT's last step is 1D FFT in the third dimension, it is more efficient to leave this dimension local and stride-1, and since the first step of the inverse FFT is to start with the third dimension 1D FFT, this format naturally fits the algorithm and results in big savings of time due to elimination of several extra transposes.

  // Set up the final global grid dimensions (these will be different from the original dimensions in one dimension due to conjugate symmetry, since we are doing real-to-complex transform)

  for(i=0; i < 3;i++) 
    gdims2[i] = gdims[i];
  gdims2[0] = gdims2[0]/2+1;

  Zpencil = p3dfft_init_data_grid(gdims2,0,Pgrid,dmap2,mem_order2);

#ifdef CUDA
  //Set up the forward transform, based on the predefined 3D transform type and Xpencil and Zpencil. This is the planning stage, needed once as initialization.
  trans_f = p3dfft_plan_3Dtrans(Xpencil,Zpencil,type_rcc,&workspace_host,&workspace_dev,LocHost,LocHost);
  //Now set up the backward transform
  trans_b = p3dfft_plan_3Dtrans(Zpencil,Xpencil,type_ccr,&workspace_host,&workspace_dev,LocHost,LocHost);
#else
  //Set up the forward transform, based on the predefined 3D transform type and Xpencil and Zpencil. This is the planning stage, needed once as initialization

  trans_f = p3dfft_plan_3Dtrans(Xpencil,Zpencil,type_rcc,&workspace_host);
  //Now set up the backward transform
  trans_b = p3dfft_plan_3Dtrans(Zpencil,Xpencil,type_ccr,&workspace_host);
#endif

  // Find local dimensions in storage order, and also the starting position of the local array in the global array
  // Note: dimensions and global starts given by grid object are in physical coordinates, which need to be translated into storage coordinates:

  for(i=0;i<3;i++) {
    glob_start1[mem_order1[i]] = Xpencil->GlobStart[i];
    sdims1[mem_order1[i]] = Xpencil->Ldims[i];
  }

  size1 = MULT3(sdims1);//[0]*sdims1[1]*sdims1[2];

  //Now allocate initial and final arrays in physical space as real-valued 1D storage containing a contiguous 3D local array 
  IN=(double *) malloc(sizeof(double)*size1);
  FIN= (double *) malloc(sizeof(double) *size1);

  //Initialize the IN array with a sine wave in 3D, based on the starting positions of my local grid within the global grid
  init_wave(IN,gdims,sdims1,glob_start1);

  //Determine local array dimensions and allocate fourier space, complex-valued out array

  for(i=0;i<3;i++) {
    glob_start2[mem_order2[i]] = Zpencil->GlobStart[i];
    sdims2[mem_order2[i]] = Zpencil->Ldims[i];
    glob2[mem_order2[i]] = Zpencil->Gdims[i];
  }

  size2 = MULT3(sdims2);//[0]*sdims2[1]*sdims2[2];
  OUT=(double *) malloc(sizeof(double) *size2 *2);

  // Warm-up run, forward transform
  p3dfft_exec_3Dtrans_double(trans_f,IN,OUT,0);

  Nglob = MULT3(gdims);

  // timing loop

  for(i=0; i < Nrep;i++) {
    t -= MPI_Wtime();
    /* Forward R2C transform, combined with spectral derivative in idir dimension (C convention, starting with 0) */
    p3dfft_exec_3Dderiv_double(trans_f,IN,OUT,idir-1,0); // Forward real-to-complex 3D FFT
    t += MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);

    /*    compute_deriv(OUT,OUT,sdims2,glob_start2,glob2,mem_order2,idir);*/

    if(myid == 0)
      printf("Results after derivative: \n");
    print_res(OUT,gdims,sdims2,glob_start2);
    normalize(OUT,size2,gdims);
    t -= MPI_Wtime();
    p3dfft_exec_3Dtrans_double(trans_b,OUT,FIN,1); // Backward (inverse) complex-to-real 3D FFT
    t += MPI_Wtime();
  }

  /* Check that we recovered initial array */
  mydiff = check_res(FIN,gdims,sdims1,glob_start1,idir);
  diff = 0.;
  MPI_Reduce(&mydiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * Nglob *0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");
    printf("Max. diff. =%lg\n",diff);
  }

  MPI_Reduce(&t,&gtavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&t,&gtmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&t,&gtmax,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0)
    printf("Transform time (avg/min/max): %lf %lf %lf\n",gtavg/(nprocs*Nrep),gtmin/Nrep,gtmax/Nrep);

  MPI_Barrier(MPI_COMM_WORLD);

  printf("Freeing In, OUT, FIN\n");
  free(IN); free(OUT); free(FIN);

  // Clean up grid structures

  p3dfft_free_data_grid(Xpencil);
  p3dfft_free_data_grid(Zpencil);
  p3dfft_free_proc_grid(Pgrid);

  // Clean up all P3DFFT++ data

  printf("Cleaning up p3dfft\n");
  p3dfft_cleanup();

  MPI_Finalize();
}


void normalize(double *A,size_t size,int *gdims)
{
  size_t i;
  double f = 1.0/(((double) gdims[0])*((double) gdims[1])*((double) gdims[2]));
  
  for(i=0;i<size*2;i++)
    A[i] = A[i] * f;

}

/* Initialize a 3D array with 3d sine wave */
void init_wave(double *IN,int *gdims,int *ldims,int *gstart)
{
  double *sinx,*siny,*sinz,sinyz,*p;
  int x,y,z;
  double twopi = atan(1.0)*8.0;

  sinx = (double *) malloc(sizeof(double)*gdims[0]);
  siny = (double *) malloc(sizeof(double)*gdims[1]);
  sinz = (double *) malloc(sizeof(double)*gdims[2]);

   for(z=0;z < ldims[2];z++)
     sinz[z] = sin((z+gstart[2])*twopi/gdims[2]);
   for(y=0;y < ldims[1];y++)
     siny[y] = sin((y+gstart[1])*twopi/gdims[1]);
   for(x=0;x < ldims[0];x++)
     sinx[x] = sin((x+gstart[0])*twopi/gdims[0]);

   p = IN;
   for(z=0;z < ldims[2];z++)
     for(y=0;y < ldims[1];y++) {
       sinyz = siny[y]*sinz[z];
       for(x=0;x < ldims[0];x++)
          *p++ = sinx[x]*sinyz;
     }

   free(sinx); free(siny); free(sinz);
}

void print_res(double *A,int *gdims,int *ldims,int *gstart)
{
  int x,y,z;
  double *p;
  double Nglob;
  int imo[3],i,j;
  
  Nglob = gdims[0]*gdims[1];
  Nglob *= gdims[2];
  p = A;
  for(z=0;z < ldims[2];z++)
    for(y=0;y < ldims[1];y++)
      for(x=0;x < ldims[0];x++) {
	if(fabs(*p) > Nglob *1.25e-4 || fabs(*(p+1))  > Nglob *1.25e-4) 
    printf("(%d %d %d) %lg %lg\n",x+gstart[0],y+gstart[1],z+gstart[2],*p,*(p+1));
	p+=2;
      }
}

/* Assumes mem_order = (0,1,2), i.e. physical = storage */
double check_res(double *A,int gdims[3],int ldims[3],int glob_start[3],int idir)
{
  int x,y,z;
  double *p,mydiff;
  int imo[3],i,j,k;
  double *cosx,*cosy,*cosz,yz,diff;
  double *sinx,*siny,*sinz;

  double twopi = atan(1.0)*8.0;

  cosx = (double *) malloc(sizeof(double)*gdims[0]);
  cosy = (double *) malloc(sizeof(double)*gdims[1]);
  cosz = (double *) malloc(sizeof(double)*gdims[2]);
  sinx = (double *) malloc(sizeof(double)*gdims[0]);
  siny = (double *) malloc(sizeof(double)*gdims[1]);
  sinz = (double *) malloc(sizeof(double)*gdims[2]);

  for(z=0;z < ldims[2];z++) {
    cosz[z] = cos((z+glob_start[2])*twopi/gdims[2]);
    sinz[z] = sin((z+glob_start[2])*twopi/gdims[2]);
  }
  for(y=0;y < ldims[1];y++) {
    cosy[y] = cos((y+glob_start[1])*twopi/gdims[1]);
    siny[y] = sin((y+glob_start[1])*twopi/gdims[1]);
  }
  for(x=0;x < ldims[0];x++) {
    cosx[x] = cos((x+glob_start[0])*twopi/gdims[0]);
    sinx[x] = sin((x+glob_start[0])*twopi/gdims[0]);
  }

  mydiff = 0.;
  p = A;
  switch(idir) {
  case 1:

    for(z=0;z < ldims[2];z++)
      for(y=0;y < ldims[1];y++) {
	yz = sinz[z] * siny[y];
	for(x=0;x < ldims[0];x++) {
	  if((diff=fabs(*p - cosx[x] * yz)) > mydiff)
	    mydiff = diff;
	  p++;
	}
      }
    break;	

  case 2:

    for(z=0;z < ldims[2];z++)
      for(y=0;y < ldims[1];y++) {
	yz = sinz[z] * cosy[y];
	for(x=0;x < ldims[0];x++) {
	  if((diff=fabs(*p - sinx[x] * yz)) > mydiff)
	    mydiff = diff;
	  p++;
	}
      }
    break;
	
  case 3:

    for(z=0;z < ldims[2];z++)
      for(y=0;y < ldims[1];y++) {
	yz = cosz[z] * siny[y];
	for(x=0;x < ldims[0];x++) {
	  if((diff=fabs(*p - sinx[x] * yz)) > mydiff)
	    mydiff = diff;
	  p++;
	}
      }
	
  }
      
  return(mydiff);
}

void write_buf(double *buf,char *label,int sz[3],int mo[3], int taskid) {
  int i,j,k;
  FILE *fp;
  char str[80],filename[80];
  double *p= buf;

  strcpy(filename,label);
  sprintf(str,".%d",taskid);
  strcat(filename,str);
  fp=fopen(filename,"w");
  for(k=0;k<sz[mo[2]];k++)
    for(j=0;j<sz[mo[1]];j++)
      for(i=0;i<sz[mo[0]];i++) {
	if(abs(*p) > 1.e-7) {
	  fprintf(fp,"(%d %d %d) %lg\n",i,j,k,*p);
	}
	p++;
      }
  fclose(fp); 
}
