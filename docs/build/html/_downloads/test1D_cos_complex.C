/*
This program exemplifies the use of 1D transforms in P3DFFT++, using cosine 1D transform (DCT-1), for complex-valued arrays. 1D transforms are performed on 3D arrays, in the dimension specified as an argument. This could be an isolated 1D transform or a stage in a multidimensional transform. This function can do local transposition, i.e. arbitrary input and output memory ordering. However it does not do an inter-processor transpose (see test_transMPI for that). 

This program initializes a 3D array with a 1D cosine wave, then
performs cosine transform twice and checks that
the results are correct, namely the same as in the start except
for a normalization factor. It can be used both as a correctness
test and for timing the library functions.

The program expects 'trans.in' file in the working directory, with
a single line of numbers : Nx,Ny,Nz,dim,Nrep,mem-order-in(1)-(3),mem-order-out(1)-(3). 
Here 
  Nx,Ny,Nz are 3D grid dimensions (Note: the dimension of cosine transform must be odd).
  dim is the dimension of 1D transform (valid values are 0 through 2, and the logical dimension is specified, i.e. actual storage dimension may be different as specified by mem-order mapping)
  Nrep is the number of repititions. 
  mem-order-in are 3 values for the memory order of the input grid, valid values of each is 0 - 2, not repeating.
  mem-order-out is the memory order of the output grid. 

Optionally a file named 'dims' can also be provided to guide in the choice
of processor geometry in case of 2D decomposition. It should contain
two numbers in a line, with their product equal to the total number
of tasks. Otherwise processor grid geometry is chosen automatically.
For better performance, experiment with this setting, varying
iproc and jproc. In many cases, minimizing iproc gives best results.
Setting it to 1 corresponds to one-dimensional decomposition.

If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu
*/

#include "p3dfft.h"
#include "compiler_check.h"
#include <math.h>
#include <stdio.h>

using namespace p3dfft;

void init_wave1D(double *,int[3],int *,int,int);
void print_res(double *,int *,int *,int *, int);
void normalize(double *,long int,int *,int);
double check_res(double*,double *,int *);

main(int argc,char **argv)
{

using namespace p3dfft;

  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3],gdims2[3];
  int pgrid1[3],pgrid2[3];
  int proc_order[3];
  int mem_order1[3];
  int mem_order2[3];
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  int *ldims,*ldims2;
  long int size1,size2;
  double *IN;
  int glob_start1[3],glob_start2[3];
  double *OUT;
  int type_ids1;
  int type_ids2;
  Type3D type_rcc,type_ccr;
  double t=0.;
  double *FIN;
  double mydiff;
  double diff = 0.0;
  double gtavg=0.;
  double gtmin=INFINITY;
  double gtmax = 0.;
  int pdims[2],nx,ny,nz,n,dim,cnt;
  FILE *fp;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  // Read input parameters

   if(myid == 0) {
     printf("P3DFFT++ test1. Running on %d cores\n",nprocs);
     printf("GitVersion = %s\n", GIT_VERSION);
     printf("GitDate = %s\n", GIT_DATE);
     printf("Executable, %s, was compiled with %s (version %d) on %s at %s\n", __FILE__, COMPILER_DETECTED, COMPILER_V_DETECTED, __DATE__, __TIME__);
     if((fp=fopen("trans.in", "r"))==NULL){
        printf("Cannot open input file. Setting to default nx=ny=nz=128, dim=0, n=1.\n");
        nx=ny=nz=128; Nrep=1;dim=0;
     } else {
        fscanf(fp,"%d %d %d %d %d\n",&nx,&ny,&nz,&dim,&Nrep);
        fscanf(fp,"%d %d %d\n",mem_order1,mem_order1+1,mem_order1+2);
        fscanf(fp,"%d %d %d\n",mem_order2,mem_order2+1,mem_order2+2);
        fclose(fp);
     }
     printf("C++ cosine transform, complex values, 1D wave input, 1D FFT\n");
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n dimension of transform: %d\n%d repetitions\n",nx,ny,nz,dim,n);
#else
     printf("Single precision\n (%d %d %d) grid\n dimension of transform %d\n%d repetitions\n",nx,ny,nz,dim,n);
#endif
   }

   // Broadcast input parameters

   MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&Nrep,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&dim,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mem_order1,3,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&mem_order2,3,MPI_INT,0,MPI_COMM_WORLD);

  //! Establish 2D processor grid decomposition, either by readin from file 'dims' or by an MPI default

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

   if(myid == 0)
      printf("Using processor grid %d x %d\n",pdims[0],pdims[1]);

  // Set up work structures for P3DFFT

  setup();

  //Set up transform types for 1D cosine transform

  type_ids1 = type_ids2 = DCT1_COMPLEX_D;

  //Set up global dimensions of the grid

  gdims[0] = nx;
  gdims[1] = ny;
  gdims[2] = nz;

  cnt = 1;
  for(i=0; i < 3;i++) 
    proc_order[i] = i;

  p1 = pdims[0];
  p2 = pdims[1];

  // Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively

  cnt=0;
  for(i=0;i<3;i++)
    if(i == dim)
      pgrid1[i] = 1;
    else
      pgrid1[i] = pdims[cnt++];

  // Set up the final global grid dimensions (these will be different from the original dimensions in one dimension since we are doing real-to-complex transform)

  for(i=0; i < 3;i++) 
    gdims2[i] = gdims[i];

  //Initialize initial and final grids, based on the above information
  // There is no conjugate symmetry (-1) since the arrays are real-valued
  grid grid1(gdims,-1,pgrid1,proc_order,mem_order1,MPI_COMM_WORLD); 

  grid grid2(gdims2,-1,pgrid1,proc_order,mem_order2,MPI_COMM_WORLD); 

  //Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization.
  transplan<complex_double,complex_double> trans_f(grid1,grid2,type_ids1,dim);

  //Now set up the backward transform

  transplan<complex_double,complex_double> trans_b(grid2,grid1,type_ids2,dim);

  //Determine local array dimensions. 

  ldims = grid1.ldims;
  size1 = ldims[0]*ldims[1]*ldims[2];

  //Now allocate initial and final arrays in physical space (real-valued)
  IN=(double *) malloc(sizeof(double)*size1*2);
  FIN= (double *) malloc(sizeof(double) *size1*2);

  int sdims1[3],sdims2[3];
  for(i=0;i<3;i++) {
    sdims1[mem_order1[i]] = grid1.ldims[i];
    sdims2[mem_order2[i]] = grid2.ldims[i];
    glob_start1[mem_order1[i]] = grid1.glob_start[i];
    glob_start2[mem_order2[i]] = grid2.glob_start[i];
  }

  //Initialize the IN array with a sine wave in 3D

  int ld = mem_order1[dim];  // Storage mapping of dimension of transform
  init_wave1D(IN,gdims,sdims1,dim,ld);

  //Determine local array dimensions and allocate fourier space, complex-valued out array

  ldims2 = grid2.ldims;
  size2 = ldims2[0]*ldims2[1]*ldims2[2];
  OUT=(double *) malloc(sizeof(double) *size2*2);


  // Execution of forward transform

  trans_f.exec((char *) IN,(char *) OUT,false);

  Nglob = gdims[0]*gdims[1]*gdims[2];

  if(myid == 0)
    printf("Results of forward transform: \n");
  print_res(OUT,gdims,sdims2,glob_start2,dim);
  normalize(OUT,sdims2[0]*sdims2[1]*sdims2[2],gdims,dim);

  // Execution of backward transform
  trans_b.exec((char *) OUT,(char *) FIN,true);

  mydiff = check_res(IN,FIN,sdims1);
  //  printf("%d: my diff =%lf\n",myid,mydiff);
  diff = 0.;
  MPI_Reduce(&mydiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * Nglob *0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");
    printf("Max. diff. =%lg\n",diff);
  }

  free(IN); free(OUT); free(FIN);

  // Clean up all P3DFFT++ data


  p3dfft_cleanup();

  MPI_Finalize();
}

void normalize(double *A,long int size,int *gdims,int dim)
{
  long int i;
  double f = 0.5/((double) (gdims[dim]-1));
  
  for(i=0;i<size*2;i++)
    A[i] = A[i] * f;

}

// Initialize 3D array with a 1D cosine wave in the specified dimension
void init_wave1D(double *IN,int *gdims,int *sdims,int dim, int ld)
{
  double *cos_coords,*p;
  int i,x,y,z;
  double pi = atan(1.0)*4.0;

  cos_coords = (double *) malloc(sizeof(double)*gdims[dim]);

  for(x=0;x < gdims[dim];x++)
    cos_coords[x] = cos(x*pi/(gdims[dim]-1));

   p = IN;
   switch(ld) {
   case 0:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++) {
	   *p++ = cos_coords[x];
	   *p++ = 0.0;
	 }
       
     break;

   case 1:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++) {
	   *p++ = cos_coords[y]; 
	   *p++ = 0.0;
	 }
       
     break;
     
   case 2:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++) {
	   *p++ = cos_coords[z];
	   *p++ = 0.0;
	 }
       
     break;
   default:
     break;
   }

   free(cos_coords); 
}

void print_res(double *A,int *gdims,int *sdims,int *gstart, int dim)
{
  int x,y,z;
  double *p;
  double N;

  N = gdims[dim];
  p = A;
  for(z=0;z < sdims[2];z++)
    for(y=0;y < sdims[1];y++)
      for(x=0;x < sdims[0];x++) {
	if(fabs(*p) + fabs(*(p+1))> N *1.25e-4) 
	  printf("(%d %d %d) %lg %lg\n",x+gstart[0],y+gstart[1],z+gstart[2],*p,*(p+1));
	p+=2;
      }
}

double check_res(double *A,double *B,int *sdims)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  p1 = A;
  p2 = B;

  mydiff = 0.0;
  for(z=0;z < sdims[2];z++)
    for(y=0;y < sdims[1];y++)
      for(x=0;x < sdims[0]*2;x++) {
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
