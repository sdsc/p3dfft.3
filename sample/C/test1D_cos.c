/*
This program exemplifies the use of 1D transforms in P3DFFT++, using cosine 1D transform (DCT-1), for real-valued arrays. 1D transforms are performed on 3D arrays, in the dimension specified as an argument. This could be an isolated 1D transform or a stage in a multidimensional transform. This function can do local transposition, i.e. arbitrary input and output memory ordering. However it does not do an inter-processor transpose (see test_transMPI for that). 


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

void init_wave1D(double *IN,int *mydims,int *gstart, int ar_dim);
void print_res(double *A,int *mydims,int *gstart, int N     );
void normalize(double *,long int,double);
double check_res(double *A,double *B,int *mydims);
void  check_res_forward(double *OUT,int sdims[3],int dim,int glob_start[3], int myid);

int main(int argc,char **argv)
{
  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3];
  int dmap[3],Pgrid;
  int mem_order1[3];
  int mem_order2[3];
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  int *ldims1,*ldims2;
  long int size1,size2;
  double *IN;
  Grid *grid1,*grid2;
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
  int pdims[3],nx,ny,nz,n,dim,cnt,ar_dim,ar_dim2,sdims1[3],sdims2[3];
  Plan3D trans_f,trans_b;
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
     printf("P3DFFT test, 1D wave input, 1D cosine transform\n");
#ifndef SINGLE_PREC
     printf("Double precision\n (%d %d %d) grid\n dimension of transform: %d\n%d repetitions\n",nx,ny,nz,dim,Nrep);
#else
     printf("Single precision\n (%d %d %d) grid\n dimension of transform %d\n%d repetitions\n",nx,ny,nz,dim,Nrep);
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

  p3dfft_setup();

  //Set up 2 transform types for 3D transforms

  type_ids1 = type_ids2 = P3DFFT_DCT1_REAL_D;

  //Set up global dimensions of the grid

  gdims[0] = nx;
  gdims[1] = ny;
  gdims[2] = nz;

  p1 = pdims[0];
  p2 = pdims[1];

  // Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively
  
  pdims[0] = 1;
  pdims[1] = p1;
  pdims[2] = p2;

  cnt=1;
  for(i=0;i<3;i++)
    if(i == dim)
      dmap[i] = 0;
    else
      dmap[i] = cnt++;

  // Set up the final global grid dimensions (these will be different from the original dimensions in one dimension since we are doing real-to-complex transform)

  Pgrid = p3dfft_init_proc_grid(pdims,MPI_COMM_WORLD);

  for(i=0; i < 3;i++) {
    if(i == dim) 
      ar_dim = mem_order1[i]; // Find storage dimension corresponding to dim
  }

  //Initialize initial and final grids, based on the above information
  // No conjugate symmetry (-1)
  grid1 = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap,mem_order1);
  grid2 = p3dfft_init_data_grid(gdims,-1,Pgrid,dmap,mem_order2);

  //Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization. Last argument is for out-of-place transform (input and output spaces different).
  trans_f = p3dfft_plan_1Dtrans(grid1,grid2,type_ids1,dim);

  //Now set up the backward transform

  trans_b = p3dfft_plan_1Dtrans(grid2,grid1,type_ids2,dim);

  //Determine local array dimensions. 

  ldims1 = grid1->Ldims;
  size1 = ldims1[0]*ldims1[1]*ldims1[2];

  // Find local dimensions in storage order, and also the starting position of the local array in the global array
  // Note: dimensions and global starts given by grid object are in physical coordinates, which need to be translated into storage coordinates:
  for(i=0;i<3;i++) {
    sdims1[mem_order1[i]] = ldims1[i];
    glob_start1[mem_order1[i]] = grid1->GlobStart[i];
  }
  //Now allocate initial and final arrays in physical space (real-valued)
  IN=(double *) malloc(sizeof(double)*size1);
  FIN= (double *) malloc(sizeof(double) *size1);

  //Initialize the IN array with a sine wave in 3D

  init_wave1D(IN,sdims1,glob_start1,ar_dim);

  //Determine local array dimensions and allocate fourier space, complex-valued out array

  ldims2 = grid2->Ldims;
  size2 = ldims2[0]*ldims2[1]*ldims2[2];
  OUT=(double *) malloc(sizeof(double) *size2);

  for(i=0;i < 3;i++) {
    sdims2[mem_order2[i]] = ldims2[i];
    glob_start2[mem_order2[i]] = grid2->GlobStart[i];
  }

  ar_dim2 = mem_order2[dim];
  //  double Nglob = ldims2[0]*ldims2[1]*ldims2[2];
  for(i=0;i<Nrep;i++) {

  // Execute forward transform
    p3dfft_exec_1Dtrans_double(trans_f,IN,OUT,0);

    //     if(myid == 0)
    // printf("Results of forward transform: \n");
    //print_res(OUT,sdims2,glob_start2,sdims1[ar_dim]);
  normalize(OUT,(long int) ldims2[0]*ldims2[1]*ldims2[2],0.5/((double) sdims1[ar_dim]-1));
  check_res_forward(OUT,sdims2,ar_dim2,glob_start2,myid);

  // Execute backward transform
  p3dfft_exec_1Dtrans_double(trans_b,OUT,FIN,1);
  }

  mydiff = check_res(IN,FIN,sdims1);

  diff = 0.;
  MPI_Reduce(&mydiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * sdims1[ar_dim] *0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");
    printf("Max. diff. =%lg\n",diff);
  }

  free(IN); free(OUT); free(FIN);

  // Clean up grid structures 

  p3dfft_free_data_grid(grid1);
  p3dfft_free_data_grid(grid2);
  p3dfft_free_proc_grid(Pgrid);

  // Clean up all P3DFFT++ data

  p3dfft_cleanup();

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
	if(it[dim] + glob_start[dim] == 1) {
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
      printf("Forward results are incorrect\n");
    else
      printf("Forward results are correct\n");
    printf("Max. diff. =%lg\n",diff);
  }
	

}

void normalize(double *A,long int size,double f)
{
  long int i;
  
  for(i=0;i<size;i++)
    A[i] = A[i] * f;

}

// Initialize a 3D array with a cosine wave in specified dimension
void init_wave1D(double *IN,int *mydims,int *gstart, int ar_dim)
{
  double *cos_coords,*p;
  int x,y,z;
  double pi = atan(1.0)*4.0;

  cos_coords = (double *) malloc(sizeof(double)*mydims[ar_dim]);

  for(x=0;x < mydims[ar_dim];x++)
    cos_coords[x] = cos(x*pi/(mydims[ar_dim]-1));

   p = IN;
   switch(ar_dim) {
   case 0:

     for(z=0;z < mydims[2];z++)
       for(y=0;y < mydims[1];y++) 
	 for(x=0;x < mydims[0];x++)
	   *p++ = cos_coords[x];
       
     break;

   case 1:

     for(z=0;z < mydims[2];z++)
       for(y=0;y < mydims[1];y++) 
	 for(x=0;x < mydims[0];x++)
	   *p++ = cos_coords[y];
       
     break;
     
   case 2:

     for(z=0;z < mydims[2];z++)
       for(y=0;y < mydims[1];y++) 
	 for(x=0;x < mydims[0];x++)
	   *p++ = cos_coords[z];
       
     break;
   default:
     break;
   }

   free(cos_coords); 
}

void print_res(double *A,int *mydims,int *gstart, int N)
{
  int x,y,z;
  double *p;
  int imo[3],i,j;
  
  p = A;
  for(z=0;z < mydims[2];z++)
    for(y=0;y < mydims[1];y++)
      for(x=0;x < mydims[0];x++) {
	if(fabs(*p) > N *1.25e-4 ) 
	  printf("(%d %d %d) %lg\n",x+gstart[0],y+gstart[1],z+gstart[2],*p);
	p++;
      }
}

double check_res(double *A,double *B,int *mydims)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  int imo[3],i,j;
  p1 = A;
  p2 = B;
  
  mydiff = 0.;
  for(z=0;z < mydims[2];z++)
    for(y=0;y < mydims[1];y++)
      for(x=0;x < mydims[0];x++) {
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
