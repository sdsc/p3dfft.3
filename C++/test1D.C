#include "p3dfft.h"
#include <math.h>
//#include "templ.C"
//#include "exec.C"

using namespace p3dfft;

  void init_wave(double *,int[3],int *,int[3]);
void print_res(complex_double *,int *,int *,int *, int *);
  void normalize(complex_double *,long int,int *);
double check_res(double*,double *,int *, int *);
void write_buf(double *buf,char *label,int sz[3],int mo[3], int taskid);

main(int argc,char **argv)
{
  using namespace p3dfft;

  int N=128;
  int Nrep = 1;
  int myid,nprocs;
  int gdims[3],gdims2[3];
  int pgrid1[3],pgrid2[3];
  int proc_order[3];
  int mem_order[3];
  int i,j,k,x,y,z,p1,p2;
  double Nglob;
  int imo1[3];
  //  void inv_mo(int[3],int[3]);
  //void write_buf(double *,char *,int[3],int[3],int);

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  cout << "P3DFFT++ test1. Running on " << nprocs << "cores" << endl;

  setup();

  cout << "Passed p3dfft_setup" << endl;
  cout << "types1D[0]: dt1=" << types1D[0]->dt1 << endl;

  for(i=0; i < 3;i++) {
    gdims[i] = N;
    proc_order[i] = mem_order[i] = i;
  }
  p1 = 1;
  p2 = nprocs;
  pgrid1[0] = 1;
  pgrid1[1] = p1;
  pgrid1[2] = p2;

  cout << "Initiating grid1" << endl;
  grid grid1(gdims,pgrid1,proc_order,mem_order,MPI_COMM_WORLD);  
  pgrid2[0] =p1;
  pgrid2[1] = p2;
  pgrid2[2] = 1;
  for(i=0; i < 3;i++) 
    gdims2[i] = gdims[i];
  gdims2[0] = gdims2[0]/2+1;
  int mem_order2[] = {2,1,0};
  cout << "Initiating grid2" << endl;
  grid grid2(gdims2,pgrid2,proc_order,mem_order2,MPI_COMM_WORLD);  

  int *ldims = &grid1.ldims[0];
  int size1 = ldims[0]*ldims[1]*ldims[2];
  cout << "Allocating IN" << endl;
  double *IN=new double[size1];
  cout << "Initiating wave" << endl;
  printf("%d: grid1 sglobal starts: %d %d %d\n",myid,grid1.glob_start[0],grid1.glob_start[1],grid1.glob_start[2]);
  init_wave(IN,gdims,ldims,grid1.glob_start);
  inv_mo(mem_order,imo1);
  write_buf(IN,"Init",ldims,imo1,myid);

  ldims = &grid2.ldims[0];
  long int size2 = ldims[0]*ldims[1]*ldims[2];
  cout << "allocating OUT" << endl;
  complex_double *OUT=new complex_double[size2];
  
  // Set up transform types for 3D transform
  int type_ids1[3] = {R2CFFT_D,CFFT_FORWARD_D,CFFT_FORWARD_D};
  int type_ids2[3] = {C2RFFT_D,CFFT_BACKWARD_D,CFFT_BACKWARD_D};

  trans_type3D type_rcc(type_ids1); 
  trans_type3D type_ccr(type_ids2); 
  // Set up 3D transforms, including stages and plans, for forward trans.
  transform3D<double,complex_double> trans_f(grid1,grid2,&type_rcc,false);
  // Set up 3D transforms, including stages and plans, for backward trans.
  transform3D<complex_double,double> trans_b(grid2,grid1,&type_ccr,false);

  //  trans_f.exec(IN,OUT,0);

  double t=0.;
  double *FIN=new double[size1];
  Nglob = gdims[0]*gdims[1]*gdims[2];

  for(i=0; i < Nrep;i++) {
    t -= MPI_Wtime();
    trans_f.exec(IN,OUT,0);
    t += MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid == 0)
      cout << "Results of forward transform: "<< endl;
    print_res(OUT,gdims,grid2.ldims,grid2.glob_start,mem_order2);
    normalize(OUT,size2,gdims);
    t -= MPI_Wtime();
    trans_b.exec(OUT,FIN,0);
    t += MPI_Wtime();
  }

  double mydiff = check_res(IN,FIN,grid1.ldims,mem_order);
  double diff = 0.0;
  MPI_Reduce(&mydiff,&diff,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid == 0) {
    if(diff > 1.0e-14 * Nglob *0.25)
      cout << "Results are incorrect" << endl;
    else
      cout << "Results are correct" << endl;
    cout << "Max. diff. =" << diff << endl;
  }

  double gtavg=0.;
  double gtmin=INFINITY;
  double gtmax = 0.;
  MPI_Reduce(&gtavg,&t,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&gtmin,&t,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&gtmax,&t,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  if(myid = 0)
    cout << "Transform time (avg/min/max): %lf %lf %lf" << gtavg/nprocs << gtmin << gtmax << endl;

  delete [] IN,OUT,FIN;
  cleanup();
  MPI_Finalize();


/*
  fopen("input","r");
  fscanf("%d %d %d %d %d\n",
  */
}

void normalize(complex_double *A,long int size,int *gdims)
{
  long int i;
  double f = 1.0/(((double) gdims[0])*((double) gdims[1])*((double) gdims[2]));
  
  for(i=0;i<size;i++)
    A[i] = A[i] * f;

}

void init_wave(double *IN,int *gdims,int *ldims,int *gstart)
{
  double *sinx,*siny,*sinz,sinyz,*p;
  int x,y,z;
  double twopi = atan(1.0)*8.0;

  sinx = new double[gdims[0]];
  siny = new double[gdims[1]];
  sinz = new double[gdims[2]];

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

   delete [] sinx,siny,sinz;
}

void print_res(complex_double *A,int *gdims,int *ldims,int *gstart, int *mo)
{
  int x,y,z;
  complex_double *p;
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
	if(abs(*p) > Nglob *1.25e-4) 
	  printf("(%d %d %d) %lg %lg\n",x+gstart[imo[0]],y+gstart[imo[1]],z+gstart[imo[2]],p->real(),p->imag());
	p++;
      }
}

double check_res(double *A,double *B,int *ldims, int *mo)
{
  int x,y,z;
  double *p1,*p2,mydiff;
  p1 = A;
  p2 = B;

  int imo[3],i,j;
  
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
