
/*
This program exemplifies using P3DFFT++ library for 3D real-to-complex FFT. 

This program initializes a 3D array with a 3D sine wave, then
performs 3D forward Fourier transform, then backward transform,
and checks that
the results are correct, namely the same as in the start except
for a normalization factor. It can be used both as a correctness
test and for timing the library functions.

The program expects 'stdin' file in the working directory, with
a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here 
  Nx,Ny,Nz are box (grid) dimensions.
  Ndim is the dimentionality of processor grid (1 or 2).
  Nrep is the number of repititions for the timing loop. 

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

void init_wave(double *,int[3],int *,int[3]);
void print_res(double *,int *,int *,int *);
void normalize(double *,long int,int *);
double check_res(double*,double *,int *);

main(int argc,char **argv)
{
  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3],gdims2[3];
  int pgrid1[3],pgrid2[3];
  int proc_order[3];
  int mem_order1[3];
  int mem_order2[] = {2,1,0};
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  int sdims1[3],ldims2[3];
  long int size1,size2;
  double *IN;
  Grid *grid1,*grid2;
  int glob_start1[3],glob_start2[3];
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
  int pdims[2],nx,ny,nz,n,ndim;
  Plan3D trans_f,trans_b;
  FILE *fp;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);


   if(myid == 0) {
     printf("P3DFFT++ C test program. Running on %d cores\n",nprocs);
     printf("GitVersion = %s\n", GIT_VERSION);
     printf("GitDate = %s\n", GIT_DATE);
     printf("Executable, %s, was compiled with %s (version %d) on %s at %s\n", __FILE__, COMPILER_DETECTED, COMPILER_V_DETECTED, __DATE__, __TIME__);
     if((fp=fopen("stdin", "r"))==NULL){
        printf("Cannot open file. Setting to default nx=ny=nz=128, ndim=2, n=1.\n");
        nx=ny=nz=128; Nrep=1;ndim=2;
     } else {
        fscanf(fp,"%d %d %d %d %d\n",&nx,&ny,&nz,&ndim,&Nrep);
        fclose(fp);
     }
     printf("P3DFFT test, 3D wave input, 3D real-to-complex FFT\n");
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n",nx,ny,nz,ndim,Nrep);
#else
     printf("Single precision\n (%d %d %d) grid\n %d proc. dimensions\n%d repetitions\n",nx,ny,nz,ndim,Nrep);
#endif
   }
   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Nrep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ndim,1,MPI_INT,0,MPI_COMM_WORLD);

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

  p3dfft_setup();

  //Set up 2 transform types for 3D transforms

  type_ids1[0] = P3DFFT_R2CFFT_D;
  type_ids1[1] = P3DFFT_CFFT_FORWARD_D;
  type_ids1[2] = P3DFFT_CFFT_FORWARD_D;

  type_ids2[0] = P3DFFT_C2RFFT_D;
  type_ids2[1] = P3DFFT_CFFT_BACKWARD_D;
  type_ids2[2] = P3DFFT_CFFT_BACKWARD_D;

  //Now initialize 3D transforms (forward and backward) with these types

  type_rcc = p3dfft_init_3Dtype(type_ids1);
  type_ccr = p3dfft_init_3Dtype(type_ids2);

  //Set up global dimensions of the grid

  gdims[0] = nx;
  gdims[1] = ny;
  gdims[2] = nz;
  for(i=0; i < 3;i++) {
    proc_order[i] = mem_order1[i] = i; // The simplest case of sequential ordering
  }

  p1 = pdims[0];
  p2 = pdims[1];

  // Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively

  pgrid1[0] = 1;
  pgrid1[1] = p1;
  pgrid1[2] = p2;
 
  //Define the final processor grid. It can be the same as initial grid, however here it is different. First, in real-to-complex and complex-to-real transforms the global grid dimensions change for example from (n0,n1,n2) to (n0/2+1,n1,n2), since most applications attempt to save memory by using the conjugate symmetry of the Fourier transform of real data. Secondly, the final grid may have different processor distribution and memory ordering, since for example many applications with convolution and those solving partial differential equations do not need the initial grid configuration in Fourier space. The flow of these applications is typically 1) transform from physical to Fourier space, 2) apply convolution or derivative calculation in Fourier space, and 3) inverse FFT to physical space. Since forward FFT's last step is 1D FFT in the third dimension, it is more efficient to leave this dimension local and stride-1, and since the first step of the inverse FFT is to start with the third dimension 1D FFT, this format naturally fits the algorithm and results in big savings of time due to elimination of several extra transposes.

  pgrid2[0] =p1;
  pgrid2[1] = p2;
  pgrid2[2] = 1;

  // Set up the final global grid dimensions (these will be different from the original dimensions in one dimension since we are doing real-to-complex transform)

  for(i=0; i < 3;i++) 
    gdims2[i] = gdims[i];
  gdims2[0] = gdims2[0]/2+1;

  //Initialize initial and final grids, based on the above information

  grid1 = p3dfft_init_grid(gdims,-1,pgrid1,proc_order,mem_order1,MPI_COMM_WORLD); 

  grid2 = p3dfft_init_grid(gdims2,0,pgrid2,proc_order,mem_order2,MPI_COMM_WORLD); 

  //Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization.
  trans_f = p3dfft_plan_3Dtrans(grid1,grid2,type_rcc);

  //Now set up the backward transform

  trans_b = p3dfft_plan_3Dtrans(grid2,grid1,type_ccr);

  // Find local dimensions in storage order, and also the starting position of the local array in the global array

  for(i=0;i<3;i++) {
    glob_start1[mem_order1[i]] = grid1->glob_start[i];
    sdims1[mem_order1[i]] = grid1->ldims[i];
  }

  size1 = sdims1[0]*sdims1[1]*sdims1[2];

  //Now allocate initial and final arrays in physical space as real-valued 1D storage containing a contiguous 3D local array 
  IN=(double *) malloc(sizeof(double)*size1);
  FIN= (double *) malloc(sizeof(double) *size1);

  //Initialize the IN array with a sine wave in 3D, based on the starting positions of my local grid within the global grid


  init_wave(IN,gdims,sdims1,glob_start1);

  //Determine local array dimensions and allocate fourier space, complex-valued out array

  for(i=0;i<3;i++) {
    glob_start2[mem_order2[i]] = grid2->glob_start[i];
    ldims2[mem_order2[i]] = grid2->ldims[i];
  }

  size2 = ldims2[0]*ldims2[1]*ldims2[2];
  OUT=(double *) malloc(sizeof(double) *size2 *2);

  // Warm-up run, forward transform
  p3dfft_exec_3Dtrans_double(trans_f,IN,OUT,0);

  Nglob = gdims[0]*gdims[1];
  Nglob *= gdims[2];

  // timing loop

  for(i=0; i < Nrep;i++) {
    t -= MPI_Wtime();
    p3dfft_exec_3Dtrans_double(trans_f,IN,OUT,0); // Forward real-to-complex 3D FFT
    t += MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0)
      printf("Results of forward transform: \n");
    print_res(OUT,gdims,ldims2,glob_start2);
    normalize(OUT,size2,gdims);
    t -= MPI_Wtime();
    p3dfft_exec_3Dtrans_double(trans_b,OUT,FIN,1); // Backward (inverse) complex-to-real 3D FFT
    t += MPI_Wtime();
  }

  mydiff = check_res(IN,FIN,sdims1);
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
  free(IN); free(OUT); free(FIN);

  // Clean up grid structures

  p3dfft_free_grid(grid1);
  p3dfft_free_grid(grid2);

  // Clean up all P3DFFT++ data

  p3dfft_cleanup();

  MPI_Finalize();
}

void normalize(double *A,long int size,int *gdims)
{
  long int i;
  double f = 1.0/(((double) gdims[0])*((double) gdims[1])*((double) gdims[2]));
  
  for(i=0;i<size*2;i++)
    A[i] = A[i] * f;

}

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

double check_res(double *A,double *B,int *ldims)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  int imo[3],i,j;
  p1 = A;
  p2 = B;
  
  mydiff = 0.;
  for(z=0;z < ldims[2];z++)
    for(y=0;y < ldims[1];y++)
      for(x=0;x < ldims[0];x++) {
	if(fabs(*p1 - *p2) > mydiff)
	  mydiff = fabs(*p1-*p2);
	p1++;
	p2++;
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
