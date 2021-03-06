/*
This program exemplifies the use of 1D transforms in P3DFFT++, using sine 1D transform (DST-1), for real-valued arrays. 1D transforms are performed on 3D arrays, in the dimension specified as an argument. This could be an isolated 1D transform or a stage in a multidimensional transform. This function can do local transposition, i.e. arbitrary input and output memory ordering. However it does not do an inter-processor transpose (see test_transMPI for that). 

This program initializes a 3D array with a 1D sine wave, then
performs sine transform twice and checks that
the results are correct, namely the same as in the start except
for a normalization factor. It can be used both as a correctness
test and for timing the library functions.

The program expects 'trans.in' file in the working directory, with
a single line of numbers : Nx,Ny,Nz,dim,Nrep,mem-order-in(1)-(3),mem-order-out(1)-(3). 
Here 
  Nx,Ny,Nz are 3D grid dimensions (Note: the dimension of sine transform must be odd).
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
void  check_res_forward(double *OUT,int sdims[3],int dim,int glob_start[3], int myid);

int main(int argc,char **argv)
{

using namespace p3dfft;

  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3];
  int dmap[3];
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
  int pdims[3],nx,ny,nz,n,dim,cnt;
  FILE *fp;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  {
  // Read input parameters

   if(myid == 0) {
     printf("P3DFFT++ test. Running on %d cores\n",nprocs);
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
     printf("C++ test, sine transform DST-1, 1D wave input, 1D FFT\n");
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
       fscanf(fp,"%d %d\n",pdims+1,pdims+2);
       fclose(fp);
       if(pdims[1]*pdims[2] != nprocs)
          pdims[2] = nprocs / pdims[1];
     }
     else {
       if(myid == 0)
          printf("Creating proc. grid with mpi_dims_create\n");
       pdims[1]=pdims[2]=0;
       MPI_Dims_create(nprocs,2,pdims+1);
       if(pdims[1] > pdims[2]) {
          pdims[1] = pdims[2];
          pdims[2] = nprocs/pdims[1];
       }
     }


   if(myid == 0)
      printf("Using processor grid %d x %d\n",pdims[0],pdims[1]);

  // Set up work structures for P3DFFT

  setup();

  //Set up transform types for 1D sine transform

  type_ids1 = type_ids2 = DST1_REAL_D;

  //Set up global dimensions of the grid

  gdims[0] = nx;
  gdims[1] = ny;
  gdims[2] = nz;

  pdims[0] = 1;

  // Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively

  ProcGrid pgrid(pdims,MPI_COMM_WORLD);

  //Initialize initial and final grids, based on the above information
  // For initial grid, intended for real-valued array, there is no conjugate symmetry, i.e. -1
  cnt=1;
  for(i=0;i<3;i++)
    if(i == dim)
      dmap[i] = 0;
    else
      dmap[i] = cnt++;


  // Initialize the initial grid 

  DataGrid grid1(gdims,-1,&pgrid,dmap,mem_order1);

  // For final grid, intended for complex-valued array, there will be conjugate symmetry in the dimension of the transform (dim) since it is a R2C transform
  DataGrid grid2(gdims,-1,&pgrid,dmap,mem_order2);

  //Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization.
  transplan<double,double> trans_f(grid1,grid2,type_ids1,dim);

  //Now set up the backward transform

  transplan<double,double> trans_b(grid2,grid1,type_ids2,dim);

  //Determine local array dimensions. 

  ldims = grid1.Ldims;
  size1 = ldims[0]*ldims[1]*ldims[2];

  //Now allocate initial and final arrays in physical space (real-valued)
  IN=(double *) malloc(sizeof(double)*size1);
  FIN= (double *) malloc(sizeof(double) *size1);

  int sdims1[3],sdims2[3];
  for(i=0;i<3;i++) {
    sdims1[mem_order1[i]] = grid1.Ldims[i];
    sdims2[mem_order2[i]] = grid2.Ldims[i];
    glob_start1[mem_order1[i]] = grid1.GlobStart[i];
    glob_start2[mem_order2[i]] = grid2.GlobStart[i];
  }

  //Initialize the IN array with a sine wave in 3D

  int ld = mem_order1[dim];  // Storage mapping of dimension of transform
  init_wave1D(IN,gdims,sdims1,dim,ld);
  int ar_dim2 = mem_order2[dim];

  //Determine local array dimensions and allocate fourier space, complex-valued out array

  ldims2 = grid2.Ldims;
  size2 = ldims2[0]*ldims2[1]*ldims2[2];
  OUT=(double *) malloc(sizeof(double) *size2);


  // Execution of forward transform

  trans_f.exec((char *) IN,(char *) OUT,false);

  Nglob = gdims[0]*gdims[1]*gdims[2];

  //  if(myid == 0)
  //  printf("Results of forward transform: \n");
  //print_res(OUT,gdims,sdims2,glob_start2,dim);
  normalize(OUT,sdims2[0]*sdims2[1]*sdims2[2],gdims,dim);
  check_res_forward(OUT,sdims2,ar_dim2,glob_start2,myid);

  // Execution of backward transform
  trans_b.exec((char *) OUT,(char *) FIN,true);

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

  free(IN); free(OUT); free(FIN);

  // Clean up all P3DFFT++ data

  p3dfft_cleanup();

  }

  MPI_Finalize();
}

void  check_res_forward(double *OUT,int sdims[3],int dim,int glob_start[3], int myid)
{
  int it[3];
  double ans,d,diff,cdiff=0;
  double *p=OUT;

  for(it[2]=0;it[2] < sdims[2];it[2]++)
    for(it[1]=0;it[1] < sdims[1];it[1]++)
      for(it[0]=0;it[0] < sdims[0];it[0]++) {
	if(it[dim] + glob_start[dim] == 0) {
	  ans = 0.5;
	}
	else ans = 0.0;
	d = fabs(*p++ - ans);
	if(cdiff < d)
	  cdiff = d;
      }

	   
  diff = 0.;
  MPI_Reduce(&cdiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * sdims[dim] *0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");
    printf("Max. diff. =%lg\n",diff);
  }
	

}

void normalize(double *A,long int size,int *gdims,int dim)
{
  long int i;
  double f = 0.5/((double) (gdims[dim]+1));
  
  for(i=0;i<size;i++)
    A[i] = A[i] * f;

}

// Initialize 3D array with a 1D sine wave in the specified dimension
void init_wave1D(double *IN,int *gdims,int *sdims, int dim, int ld)
{
  double *sin_coords,*p;
  int i,x,y,z;
  double pi = atan(1.0)*4.0;

  sin_coords = (double *) malloc(sizeof(double)*gdims[dim]);

  for(x=0;x < gdims[dim];x++)
    sin_coords[x] = sin((x+1)*pi/(gdims[dim]+1));

   p = IN;
   switch(ld) {
   case 0:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++)
	   *p++ = sin_coords[x];
       
     break;

   case 1:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++)
	   *p++ = sin_coords[y];
       
     break;
     
   case 2:

     for(z=0;z < sdims[2];z++)
       for(y=0;y < sdims[1];y++) 
	 for(x=0;x < sdims[0];x++)
	   *p++ = sin_coords[z];
       
     break;
   default:
     break;
   }

   free(sin_coords); 
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
	if(fabs(*p) > N *1.25e-4) 
    printf("(%d %d %d) %lg\n",x+gstart[0],y+gstart[1],z+gstart[2],*p);
	p++;
      }
}

double check_res(double *A,double *B,int *sdims)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  p1 = A;
  p2 = B;

  mydiff = 0.;
  for(z=0;z < sdims[2];z++)
    for(y=0;y < sdims[1];y++)
      for(x=0;x < sdims[0];x++) {
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
