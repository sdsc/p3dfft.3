#include "p3dfft.h"
//#include "defs.h"
#include <math.h>
#include <stdio.h>
//#include "templ.C"
//#include "exec.C"

  void init_wave(double *,int[3],int *,int[3]);
void print_res(double *,int *,int *,int *, int *);
  void normalize(double *,long int,int *);
double check_res(double*,double *,int *, int *);

main(int argc,char **argv)
{
  int N=64;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3],gdims2[3];
  int pgrid1[3],pgrid2[3];
  int proc_order[3];
  int mem_order[3];
  int mem_order2[] = {2,1,0};
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  //  void inv_mo(int[3],int[3]);
  //void write_buf(double *,char *,int[3],int[3],int);
  int *ldims,*ldims2;
  long int size1,size2;
  double *IN;
  Grid *grid1,*grid2;
  int *glob_start,*glob_start2;
  double *OUT;
  // Set up transform types for 3D transform
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
  int pdims[2];
  Plan3D trans_f,trans_b;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  printf("P3DFFT++ test1. Running on %d cores\n",nprocs);

  p3dfft_setup();

  type_ids1[0] = P3DFFT_R2CFFT_D;
  type_ids1[1] = P3DFFT_CFFT_FORWARD_D;
  type_ids1[2] = P3DFFT_CFFT_FORWARD_D;
  type_ids2[0] = P3DFFT_C2RFFT_D;
  type_ids2[1] = P3DFFT_CFFT_BACKWARD_D;
  type_ids2[2] = P3DFFT_CFFT_BACKWARD_D;

  for(i=0; i < 3;i++) {
    gdims[i] = N;
    proc_order[i] = mem_order[i] = i;
  }
  /*
  pdims[0]=pdims[1]=0;
  MPI_Dims_create(nprocs,2,pdims);
  //  p1 = floor(sqrt(nprocs));
  //p2 = nprocs / p1;
  p1 = pdims[0];
  p2 = pdims[1];
  */
  pgrid1[0] = 1;
  pgrid1[1] = 1;
  pgrid1[2] = nprocs;


  /*  grid1 = (Grid *) malloc(sizeof(Grid));
  for(i=0; i < 3;i++) {
    grid1->gdims[i] = gdims[i];
    grid1->pgrid[i] = pgrid1[i];
    grid1->proc_order[i] = proc_order[i];
    grid1->mem_order[i] = mem_order[i];
  }
  grid1->mpi_comm_glob = MPI_COMM_WORLD;
  */
 
  pgrid2[0] = 1;
  pgrid2[1] = nprocs;
  pgrid2[2] = 1;
  for(i=0; i < 3;i++) 
    gdims2[i] = gdims[i];
  gdims2[0] = gdims2[0]/2+1;

  /*
  grid2 = (Grid *) malloc(sizeof(Grid));
  for(i=0; i < 3;i++) {
    grid2->gdims[i] = gdims2[i];
    grid2->pgrid[i] = pgrid2[i];
    grid2->proc_order[i] = proc_order[i];
    grid2->mem_order[i] = mem_order2[i];
  }
  grid2->mpi_comm_glob = MPI_COMM_WORLD;
  */

  printf("Initiating grid1\n");
  
  grid1 = p3dfft_init_grid(gdims,pgrid1,proc_order,mem_order,MPI_COMM_WORLD); 

  printf("Initiating grid2\n");
  grid2 = p3dfft_init_grid(gdims2,pgrid2,proc_order,mem_order2,MPI_COMM_WORLD); 

  printf("Initializing type RCC\n");
  type_rcc = p3dfft_init_3Dtype(type_ids1);
  printf("Initializing type CCR\n");
  type_ccr = p3dfft_init_3Dtype(type_ids2);

  printf("Plan rcc\n");
  // Set up 3D transforms, including stages and plans, for forward trans.
  trans_f = p3dfft_plan_3Dtrans(grid1,grid2,type_rcc,1);

  // Set up 3D transforms, including stages and plans, for backward trans.
  printf("Plan ccr\n");
  trans_b = p3dfft_plan_3Dtrans(grid2,grid1,type_ccr,1);

  ldims = grid1->ldims;
  size1 = ldims[0]*ldims[1]*ldims[2];
  printf("Allocating IN\n");
  IN=(double *) malloc(sizeof(double)*size1);

  printf("Initiating wave\n");
  glob_start = grid1->glob_start;
  //  printf("%d: grid1 sglobal starts: %d %d %d\n",myid,glob_start[0],glob_start[1],glob_start[2]);
  init_wave(IN,gdims,ldims,glob_start);
  //inv_mo(mem_order,imo1);
  //write_buf(IN,"Init",ldims,imo1,myid);

  ldims2 = grid2->ldims;
  glob_start2 = grid2->glob_start;
  size2 = ldims2[0]*ldims2[1]*ldims2[2];
  printf("allocating OUT, size=%d\n",size2);
  OUT=(double *) malloc(sizeof(double) *size2 *2);

  // Warm-up run, forward transform
  p3dfft_exec_3Dtrans_double(trans_f,IN,OUT,0);

  FIN= (double *) malloc(sizeof(double) *size1);
  Nglob = gdims[0]*gdims[1]*gdims[2];

  for(i=0; i < Nrep;i++) {
    printf("Executing rcc\n");
    t -= MPI_Wtime();
    p3dfft_exec_3Dtrans_double(trans_f,IN,OUT,0);
    t += MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0)
      printf("Results of forward transform: \n");
    print_res(OUT,gdims,ldims2,glob_start2,mem_order2);
    normalize(OUT,size2,gdims);
    printf("Executing ccr\n");
    t -= MPI_Wtime();
    p3dfft_exec_3Dtrans_double(trans_b,OUT,FIN,0);
    t += MPI_Wtime();
  }

  mydiff = check_res(IN,FIN,ldims,mem_order);
  MPI_Reduce(&mydiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * Nglob *0.25)
      printf("Results are incorrect\n");
    else
      printf("Results are correct\n");
    printf("Max. diff. =%d\n",diff);
  }

  MPI_Reduce(&gtavg,&t,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&gtmin,&t,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&gtmax,&t,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid = 0)
    printf("Transform time (avg/min/max): %lf %lf %lf\n",gtavg/nprocs,gtmin,gtmax);

  free(IN); free(OUT); free(FIN);
  p3dfft_free_grid(grid1);
  p3dfft_free_grid(grid2);
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

void print_res(double *A,int *gdims,int *ldims,int *gstart, int *mo)
{
  int x,y,z;
  double *p;
  double Nglob;
  int imo[3],i,j;
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(mo[i] == j)
	imo[j] = i;

  Nglob = gdims[0]*gdims[1]*gdims[2];
  p = A;
  for(z=0;z < ldims[imo[2]];z++)
    for(y=0;y < ldims[imo[1]];y++)
      for(x=0;x < ldims[imo[0]];x++) {
	if(abs(*p) > Nglob *1.25e-4 || abs(*(p+1))  > Nglob *1.25e-4) 
    printf("(%d %d %d) %lg %lg\n",x+gstart[imo[0]],y+gstart[imo[1]],z+gstart[imo[2]],*p,*(p+1));
	p+=2;
      }
}

double check_res(double *A,double *B,int *ldims, int *mo)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  int imo[3],i,j;
  p1 = A;
  p2 = B;
  
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      if(mo[i] == j)
	imo[j] = i;

  mydiff = 0.;
  for(z=0;z < ldims[imo[2]];z++)
    for(y=0;y < ldims[imo[1]];y++)
      for(x=0;x < ldims[imo[0]];x++) {
	if(abs(*p1 - *p2) > mydiff)
	  mydiff = abs(*p1-*p2);
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
