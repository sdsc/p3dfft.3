#include "p3dfft.h"
#include <string.h>

int P3DFFT_EMPTY_TYPE,P3DFFT_R2CFFT_S,P3DFFT_R2CFFT_D,P3DFFT_C2RFFT_S,P3DFFT_C2RFFT_D,P3DFFT_CFFT_FORWARD_S,P3DFFT_CFFT_FORWARD_D,P3DFFT_CFFT_BACKWARD_S,P3DFFT_CFFT_BACKWARD_D,P3DFFT_COSTRAN_REAL_S,P3DFFT_COSTRAN_REAL_D,P3DFFT_SINTRAN_REAL_S,P3DFFT_SINTRAN_REAL_D,P3DFFT_COSTRAN_COMPLEX_S,P3DFFT_COSTRAN_COMPLEX_D,P3DFFT_SINTRAN_COMPLEX_S,P3DFFT_SINTRAN_COMPLEX_D,P3DFFT_CHEB_REAL_S,P3DFFT_CHEB_REAL_D,P3DFFT_CHEB_COMPLEX_S,P3DFFT_CHEB_COMPLEX_D;

namespace p3dfft {

using namespace std;


vector<Plan *> Plans;
vector<gen_trans_type *> types1D;
vector<gen_transform3D *> stored_trans3D;
vector<stage *> stored_trans1D;
vector<trans_type3D> types3D;
vector<grid> stored_grids;

  //  extern "C" {
int EMPTY_TYPE,R2CFFT_S,R2CFFT_D,C2RFFT_S,C2RFFT_D,CFFT_FORWARD_S,CFFT_FORWARD_D,CFFT_BACKWARD_S,CFFT_BACKWARD_D,COSTRAN_REAL_S,COSTRAN_REAL_D,SINTRAN_REAL_S,SINTRAN_REAL_D,COSTRAN_COMPLEX_S,COSTRAN_COMPLEX_D,SINTRAN_COMPLEX_S,SINTRAN_COMPLEX_D,CHEB_REAL_S,CHEB_REAL_D,CHEB_COMPLEX_S,CHEB_COMPLEX_D;

  // }

int padd;

void setup()
{
  const char *name;
  int isign;
  gen_trans_type *p;
  /*
  float *pf;
  double *pd;
  complex *pc;
  complex_double *pcd;
  
  get_type(pf,type_float);
  get_type(pd,type_double);
  get_type(pc,type_complex);
  get_type(pcd,type_complex_double);
  */

  int types_count=0;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type" << endl;
#endif
  name = "Empty Type";
  p = new gen_trans_type(name,0);
  delete p;
  p = new gen_trans_type(name,0);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE = types_count;
  types_count++;

  /*

  cout << "Empty Type Float" << endl;
  name = "Empty Type Float";
  p = new trans_type1D<float,float>(name,NULL);
  types1D.push_back(p);
  EMPTY_TYPE_FLOAT_ID = types_count;
  types_count++;

  cout << "Empty Type DOUBLE" << endl;
  name = "Empty Type Double";
  p = new trans_type1D<double,double>(name,NULL);
  types1D.push_back(p);
  EMPTY_TYPE_DOUBLE_ID = types_count;
  types_count++;

  cout << "Empty Type Complex" << endl;
  name = "Empty Type Complex";
  p = new trans_type1D<mycomplex,mycomplex>(name,NULL);
  types1D.push_back(p);
  EMPTY_TYPE_COMPLEX_ID = types_count;
  types_count++;

  cout << "Empty Type Double Complex" << endl;
  name = "Empty Type Double Complex";
  p = new trans_type1D<complex_double,complex_double>(name,NULL);
  types1D.push_back(p);
  EMPTY_TYPE_DOUBLE_COMPLEX_ID = types_count;
  types_count++;
  */

#ifdef DEBUG
  cout << "p3dft_setup: adding R2C single type" << endl;
#endif

  name = "Real-to-complex Fourier Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<float,mycomplex>(name,(long (*)(...) ) plan_r2c_s);
#endif
  types1D.push_back(p);
  R2CFFT_S =  types_count;
  types_count++;
  

#ifdef DEBUG
  cout << "p3dft_setup: adding R2C double type" << endl;
#endif

  name ="Real-to-complex Fourier Transform double precision";
#ifdef FFTW
  p = new trans_type1D<double,complex_double>(name,(long (*)(...) ) plan_r2c_d);
#endif
  types1D.push_back(p);
  R2CFFT_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2R single type" << endl;
#endif
  name = "Complex-to-real Fourier Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,float>(name,(long (*)(...) ) plan_c2r_s);
#endif
  types1D.push_back(p);
  C2RFFT_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2R double type" << endl;
#endif

  name = "Complex-to-real Fourier Transform, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,double>(name,(long (*)(...) ) plan_c2r_d);
#endif
  types1D.push_back(p);
  C2RFFT_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C forward single type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, forward transform, singple precision";
#ifdef FFTW
  isign = FFTW_FORWARD;
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_c2c_s,(void (*)(...) )exec_c2c_s,isign);
#endif
  types1D.push_back(p);
  CFFT_FORWARD_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C forward double type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, forward transform, double precision";
#ifdef FFTW
  isign = FFTW_FORWARD;
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_c2c_d,(void (*)(...)) exec_c2c_d,isign);
#endif

  types1D.push_back(p);
  CFFT_FORWARD_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C backward single type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, backward transform, single precision";
#ifdef FFTW
  isign = FFTW_BACKWARD;
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_c2c_s,(void (*)(...)) exec_c2c_s,isign);
#endif
  types1D.push_back(p);
  CFFT_BACKWARD_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C double backward type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, backward transform, double precision";
#ifdef FFTW
  isign = FFTW_BACKWARD;
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_c2c_d,(void (*)(...)) exec_c2c_d,isign);
#endif
  types1D.push_back(p);
  CFFT_BACKWARD_D = types_count;
  types_count++;

  /*
  name = "Real-valued Chebyshev Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_cos_s,(void (*)(...)) scheb_r);
#endif
  types1D.push_back(p);
  */

  /*
  name = "Complex-valued Chebyshev Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myscos_plan,
       mydcos_plan,scheb_c,dcheb_c);
#endif
  types1D.push_back(p);
  */

#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R single type" << endl;
#endif

  name = "Real-valued Cosine Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_cos_s);
#endif
  types1D.push_back(p);
  COSTRAN_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R double type" << endl;
#endif

  name = "Real-valued Cosine Transform, double precision";
#ifdef FFTW
p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_cos_d);
#endif
  types1D.push_back(p);
  COSTRAN_REAL_D = types_count;
  types_count++;

  /*
  name = "Complex-valued Cosine Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myscos_plan,
		  mydcos_plan,fftwf_exec_r2r,fftw_exec_r2r);
#endif
  types1D.push_back(p);
  */

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R single type" << endl;
#endif

  name = "Real-valued Sine Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_sin_s);
#endif
  types1D.push_back(p);
  SINTRAN_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R double type" << endl;
#endif

  name = "Real-valued Sine Transform, double precision";
#ifdef FFTW
p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_sin_d);
#endif
  types1D.push_back(p);
  SINTRAN_REAL_D = types_count;
  types_count++;

  /*
  name = "Complex-valued Sine Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myssin_plan,
		  mydsin_plan,fftwf_exec_r2r,fftw_exec_r2r);
#endif
  types1D.push_back(p);
  */



  /*
  name = "Real-valued Chebyshev Transform, double precision";
#ifdef FFTW
				    p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_cos_d,(void (*)(...) )dcheb_r);
#endif
  types1D.push_back(p);
  */
  /*
  name = "Complex-valued Chebyshev Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myscos_plan,
       mydcos_plan,scheb_c,dcheb_c);
#endif
  types1D.push_back(p);
  */


  /*
  name = "Complex-valued Cosine Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myscos_plan,
		  mydcos_plan,fftwf_exec_r2r,fftw_exec_r2r);
#endif
  types1D.push_back(p);
  */

  /*
  name = "Complex-valued Sine Transform";
  dt1 = 2;dt2 = 2;
#ifdef FFTW
  p = new trans_type1D(name,dt1,dt2,myssin_plan,
		  mydsin_plan,fftwf_exec_r2r,fftw_exec_r2r);
#endif
  types1D.push_back(p);
  */

}


void cleanup()
{
  
  //  vector<grid>::iterator it1=stored_grids.begin();
  stored_grids.erase(stored_grids.begin(),stored_grids.end());

  vector<gen_trans_type *>::iterator it1=types1D.begin();
  while(it1 != types1D.end()) {
    delete *it1;
    it1 = types1D.erase(it1);
  }

  //  for(vector<trans_type3D>::iterator it=types3D.begin();it != types3D.end();it++) {
    types3D.erase(types3D.begin(),types3D.end());

  //  for(vector<Plan *>::iterator it=Plans.begin();it != Plans.end();it++) {
  vector<Plan *>::iterator it=Plans.begin();
  while(it != Plans.end()) {
    delete *it;
    it = Plans.erase(it);
  }

  vector<stage *>::iterator it2=stored_trans1D.begin();
  while(it2 != stored_trans1D.end()) {
    delete *it2;
    it2 = stored_trans1D.erase(it2);
  }

  vector<gen_transform3D *>::iterator it3=stored_trans3D.begin();
  while(it3 != stored_trans3D.end()) {
    delete *it3;
    it3 = stored_trans3D.erase(it3);
  }
    

  //   types1D.erase(types1D.begin(),types1D.end());
  //types3D.erase(types3D.begin(),types3D.end());
   //Plans.erase(Plans.begin(),Plans.end());
  // stored_grids.erase(stored_grids.begin(),stored_grids.end());
   //stored_trans3D.erase(stored_trans3D.begin(),stored_trans3D.end());
   //stored_trans1D.erase(stored_trans1D.begin(),stored_trans1D.end());
}

  
/*
long plan_r2c_s(const int *n,int howmany,float *in,const int *inembed,int istride,int idist,mycomplex *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftwf_plan_many_dft_r2c(1,n,howmany,in,inembed,istride,idist,(fftwf_complex *) out,onembed,ostride,odist,fft_flag));
#endif
}

long plan_r2c_d(const int *n,int howmany,double *in,const int *inembed,int istride,int idist,complex_double *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftw_plan_many_dft_r2c(1,n,howmany,in,inembed,istride,idist,(fftw_complex *) out,onembed,ostride,odist,fft_flag));
#endif
}
long plan_c2r_s(const int *n,int howmany,mycomplex *in,const int *inembed,int istride,int idist,float *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftwf_plan_many_dft_c2r(1,n,howmany,(fftwf_complex *) in,inembed,istride,idist,out,onembed,ostride,odist,fft_flag));
#endif
}

long plan_c2r_d(const int *n,int howmany,complex_double *in,const int *inembed,int istride,int idist,double *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftw_plan_many_dft_c2r(1,n,howmany,(fftw_complex *)in,inembed,istride,idist,out,onembed,ostride,odist,fft_flag));
#endif
}
long plan_c2c_s(const int *n,int howmany,mycomplex *in,const int *inembed,int istride,int idist,mycomplex *out,const int *onembed,int ostride,int odist,int isign,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftwf_plan_many_dft(1,n,howmany,(fftwf_complex *) in,inembed,istride,idist,(fftwf_complex *) out,onembed,ostride,odist,isign,fft_flag));
#endif
}

long plan_c2c_d(const int *n,int howmany,complex_double *in,const int *inembed,int istride,int idist,complex_double *out,const int *onembed,int ostride,int odist,int isign,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((long) fftw_plan_many_dft(1,n,howmany,(fftw_complex *) in,inembed,istride,idist,(fftw_complex *) out,onembed,ostride,odist,isign,fft_flag));
#endif
}
*/

void exec_r2c_s(long plan,float *in,mycomplex *out)
{
  fftwf_execute_dft_r2c((fftwf_plan) plan,in,(fftwf_complex *) out);
}
void exec_r2c_d(long plan,double *in,complex_double *out)
{
  fftw_execute_dft_r2c((fftw_plan) plan,in,(fftw_complex *) out);
}
void exec_c2r_s(long plan,mycomplex *in,float *out)
{
  fftwf_execute_dft_c2r((fftwf_plan) plan,(fftwf_complex *) in,out);
}
void exec_c2r_d(long plan,complex_double *in,double *out)
{
  fftw_execute_dft_c2r((fftw_plan) plan,(fftw_complex *) in, out);
}
void exec_c2c_s(long plan,mycomplex *in,mycomplex *out)
{
  fftwf_execute_dft((fftwf_plan) plan,(fftwf_complex *) in,(fftwf_complex *) out);
}
void exec_c2c_d(long plan,complex_double *in,complex_double *out)
{
  fftw_execute_dft((fftw_plan) plan,(fftw_complex *) in,(fftw_complex *) out);
}


long plan_cos_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT00;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


long plan_cos_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT00;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_sin_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT00;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_sin_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT00;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}

/*
grid newgrid(const grid &gr,const trans_type1D &type,int d;)
{
  int *gdims,int *pgrid,int *proc_order;int *mem_order;
  MPI_Comm mpicomm;

  gdims = gr.gdims;
  pgrid = gr.pgrid;
  proc_order = gr.proc_order;
  mem_order = gr.mem_order;

  if(type.dt1 < type.dt2) // Real-to-complex
    gdims[d] = gdims[2]/2+1;
  else if(type.dt2 < type.dt1) // Complex-to-real
    gdims[d] = (gdims[d]-1)*2;

  grid tmpgrid = new grid(gdims,pgrid,proc_order,mem_order,mpicomm);
} 
*/


grid::grid(int gdims_[3],int pgrid_[3],int proc_order_[3],int mem_order_[3],
	   MPI_Comm mpicomm_)
	   //,int prec_,int dt_)
{
  int i,j;
  int myid[2];
  MPI_Comm mpi_comm_tmp;

  MPI_Comm_dup(mpicomm_,&mpi_comm_glob);
  nd=0;
  P[0]=P[1]=P[2]=1;
  D[0]=D[1]=D[2]=L[0]=L[1]=L[2]=-1;
  for(i=0;i <3; i++) {
    pgrid[i] = pgrid_[i];
    gdims[i] = gdims_[i];
    proc_order[i] = proc_order_[i];
    mem_order[i] = mem_order_[i];
    if(pgrid[i] > 1) {
      P[proc_order[nd]] = pgrid[i];
      D[proc_order[nd]] = i;
      nd++;
    }
  }
  if(nd == 0) {
    nd = 1;
  }

  //  prec = prec_; dt = dt_;
  pm = max(P[0],P[1]);
  pm = max(pm,P[2]);

  MPI_Comm_rank(mpicomm_,&taskid);
  MPI_Comm_size(mpicomm_,&numtasks);

  // Set up communicators
  if(nd >1) {
    int periodic[nd];
    int reorder=0;

    for(i=0;i < nd;i++)
      periodic[i]=1;
    MPI_Cart_create(mpicomm_,nd,P,periodic,reorder,&mpi_comm_cart);
    MPI_Cart_coords(mpi_comm_cart,taskid,nd,grid_id_cart);



    int remain_dims[3];
    if(nd == 2) {
      remain_dims[0] = 1;
      remain_dims[1] = 0;
      // Create COLUMN sub comm.
      MPI_Cart_sub(mpi_comm_cart,remain_dims,&mpi_comm_tmp);
      mpicomm[0]=mpi_comm_tmp;

      remain_dims[0] = 0;
      remain_dims[1] = 1;
      // Create ROW sub comm.
      MPI_Cart_sub(mpi_comm_cart,remain_dims,mpicomm+1);

      MPI_Comm_rank(mpicomm[0],myid);
      MPI_Comm_rank(mpicomm[1],myid+1);
#ifdef DEBUG
      printf("%d: myid=%d %d\n",taskid,myid[0],myid[1]);
#endif
      //int tmp;
      //MPI_Comm_size(mpicomm[0],&tmp);
      //      printf("grid: Size of mpicomm[0]=%d\n",tmp);
      
    }
    else if(nd == 3) {
      remain_dims[0] = 1;
      remain_dims[1] = 0;
      remain_dims[2] = 0;
      MPI_Cart_sub(mpi_comm_cart,remain_dims,mpicomm+2);
      remain_dims[1] = 1;
      remain_dims[0] = 0;
      remain_dims[2] = 0;
      MPI_Cart_sub(mpi_comm_cart,remain_dims,mpicomm+1);
      remain_dims[2] = 1;
      remain_dims[1] = 0;
      remain_dims[0] = 0;
      MPI_Cart_sub(mpi_comm_cart,remain_dims,mpicomm);

    }
    //    MPI_Comm_free(&mpi_comm_cart);
  }
  else { //nd=1
    grid_id_cart[0] = taskid;
    *mpicomm = mpicomm_;
    //    InitPencil1D();
  }
  j=0;
  for(i=0;i<3;i++) 
    if(pgrid[i] == 1)
      grid_id[i] = 0;
    else
      grid_id[i] = grid_id_cart[proc_order[j++]];

  /*
  st = new int[pm][3];
  sz = new int[pm][3];
  en = new int[pm][3];

  st = new int**[nd];
  sz = new int**[nd];
  en = new int**[nd];
  */
  for(i=0;i<nd;i++) {
    st[i] = new int*[P[i]];
    sz[i] = new int*[P[i]];
    en[i] = new int*[P[i]];
    for(j=0; j < P[i]; j++) {
      st[i][j] = new int[3];
      sz[i][j] = new int[3];
      en[i][j] = new int[3];
    }
  }

  InitPencil();
  is_set = true;
}

grid::grid(const grid &rhs)
{
  
  if(rhs.is_set) {
    is_set = true;
    mpi_comm_glob = rhs.mpi_comm_glob;
    //    MPI_Comm_dup(rhs.mpi_comm_cart,&mpi_comm_cart);
    //    prec = rhs.prec;
    int i,j,m,l;
    pm = rhs.pm;
    nd = rhs.nd;

    /*    st = new int[pm][3];
    sz = new int[pm][3];
    en = new int[pm][3];
    */


    for(i=0;i<3; i++) {
      gdims[i] = rhs.gdims[i];
      ldims[i] = rhs.ldims[i];
      pgrid[i] = rhs.pgrid[i];
      proc_order[i] = rhs.proc_order[i];
      mem_order[i] = rhs.mem_order[i];
      L[i] = rhs.L[i];
      D[i] = rhs.D[i];
      grid_id[i] = rhs.grid_id[i];
    }
    for(i=0;i<nd; i++) {
      mpicomm[i] = rhs.mpicomm[i];
      P[i] = rhs.P[i];
    }
    taskid = rhs.taskid;
    /*    
    st = new int**[nd];
    sz = new int**[nd];
    en = new int**[nd];
    */
    for(i=0;i<nd;i++) {
      st[i] = new int*[P[i]];
      sz[i] = new int*[P[i]];
      en[i] = new int*[P[i]];
      for(j=0; j < P[i]; j++) {
	st[i][j] = new int[3];
	sz[i][j] = new int[3];
	en[i][j] = new int[3];
	for(int k=0;k<3; k++) {
	  st[i][j][k] = rhs.st[i][j][k];
	  sz[i][j][k] = rhs.sz[i][j][k];
	  en[i][j][k] = rhs.en[i][j][k];
	}
      }
    }
    
  }
  else
    is_set = false;
}

grid::~grid() 
{
  if(is_set) {
    for(int i=0;i < nd;i++) {
      for(int j=0;j<P[i];j++) {
	delete [] st[i][j];
	delete [] sz[i][j];
	delete [] en[i][j];
      }
      delete [] st[i];
      delete [] sz[i];
      delete [] en[i];
    }
    /*
    delete [] st;
    delete [] sz;
    delete [] en;
    */
    //    MPI_Comm_free(&mpi_comm_cart);
    //for(int i=0;i<nd;i++)
    // MPI_Comm_free(&mpicomm[i]);
  }
}

void grid::InitPencil()
{
  int i,j,k,pm,l,size,nl,nu,data,proc,li,myproc,pdim,ploc,p;

  // First determine local dimension(s)
  ploc=0;
  for(k=0; k < 3; k++) {
    for(i=0;i<nd;i++)
      if(D[i] == k)
	break;
    if(i >= nd) {
      L[ploc++] = k;
      i = -1;
    }
  }
      //      i = proc_order[pdim++];

  // Distribute the dimension each communicator spans
  for(i=0;i<nd;i++) {
    j = D[i];
    data = gdims[j]; //[li];
    proc = pgrid[j];
    size = data/proc;
    nu = data%proc;
    nl = proc-nu;
    
    st[i][0][j] = 0;
    sz[i][0][j] = size;
    en[i][0][j] = size;
    for(p=1; p < nl; p++) {
      st[i][p][j] = st[i][p-1][j] +size;
      sz[i][p][j] = size;
      en[i][p][j] = en[i][p-1][j] +size;
    }
    size++;
    for(;p < proc; p++) {
      st[i][p][j] = st[i][p-1][j] +sz[i][p-1][j];
      sz[i][p][j] = size;
      en[i][p][j] = en[i][p-1][j] +sz[i][p-1][j];
    }
    //st[i][p][j] = st[i][p-1][j] + size;
    //sz[i][p][j] = size;
    en[i][p-1][j] = data;
    sz[i][p-1][j] = data - st[i][p-1][j];

    // Assign local size for this spanning dimension
    myproc = grid_id[j];
    ldims[j] = sz[i][myproc][j];
    glob_start[j] = st[i][myproc][j];
  }
  //      en[proc-1][i] = data;
      //sz[proc-1][i] = data -st[proc-1][i];
      //      ldims[li] = sz[myid[li]][li];

  // Assign sizes for local dimensions
  for(i=0;i<ploc;i++) {
    k = L[i];
    ldims[k] = gdims[k];
    glob_start[k] = 0;
    for(j=0;j<nd;j++) 
      for(p=0;p<P[j];p++) {
	sz[j][p][k] = ldims[k];
	en[j][p][k] = ldims[k];
	st[j][p][k] = 0;
      }
  }

  // Next, figure out local sizes for all remaining dimensions (non-local and non-spanning) 
  for(i=0;i<nd;i++) {
    k = D[i];
    for(j=0;j<nd;j++)
      if(j != i) { // non-spanning
	l = D[j];
	myproc = grid_id[l];
	for(p=0;p<P[i];p++) {
	  sz[i][p][l] = sz[j][myproc][l];
	  st[i][p][l] = st[j][myproc][l];
	  en[i][p][l] = en[j][myproc][l];
	}
      }	  
  }

}

  /*
void grid::InitPencil1D() 
{
  int i,j,pm,l,size,nl,nu,data,proc,li,nloc,d;
  
  d = D[0];

  // Local dims
  nloc = 0;
  for(i=0; i < 3; i++) 
    if(i != d) {
	L[nloc] = i;
	nloc++;
	sz[0][i] = en[0][i] = ldims[i] = gdims[i];
	st[0][i] = 0;
      }
  if(nloc != 2)
    printf("ERror in InitPencil1D: expected 2 local dimensions, got%d\n",nloc);

  // split dimension
  data = gdims[d];
  proc = P[0];
  size = data/proc;
  nu = data%proc;
  nl = proc-nu;

  st[0][d] = 0;
  sz[0][d] = en[0][d] = size;
  for(j=1; j < nl; j++) {
    st[j][d] = st[j-1][d] +size;
    sz[j][d] = size;
    en[j][d] = en[j-1][d] +size;
      }
  size++;
  for(;j < proc; j++) {
    st[j][d] = st[j-1][d] + size;
    sz[j][d] = size;
    en[j][d] = en[j-1][d] +size;
  }
  en[proc-1][d] = data;
  sz[proc-1][d] = data -st[proc-1][d];
  ldims[d] = sz[taskid][d];

}
  */
  
/*
variable::variable(const variable &rhs)
{
  ext = rhs.ext;
  prec = rhs.prec;
  datatype=rhs.datatype;
  decomp = rhs.decomp;
  data = rhs.data;
  alloc=false;
}
variable::variable(const grid &rhs,int precision,int dt)
{
  prec = precision; datatype=dt;
  ext = prec*dt;
  decomp = &rhs;
  int *ldims=decomp-> ldims;
  int n = ldims[0]*ldims[1]*ldims[2]*ext;
  data = new char[n]; 
  alloc=true;
}

variable::variable~()
{
  if(alloc)
    delete data;

}
*/

// Function for defining a new 1D trans_type
/*
template <class Type1,class Type2> trans_type1D<Type1,Type2>::trans_type1D(const char name[],const Plan &plan) {

  gen_trans_type(name,plan.isign,plan.fft_r2r_kind,plan.fft_flags); 

  int prec2;
  if(typeid(Type1) == type_float) {
    prec = 4;
    dt1 = 1;
  }
  else if(typeid(Type1) == type_double) {
    prec = 8;
    dt1 = 1;
  }
  else if(typeid(Type1) == type_complex) {
    prec = 4;
    dt1 = 2;
  }
  else if(typeid(Type1) == type_complex_double) {
    prec = 8;
    dt1 = 2;
  }

  if(typeid(Type2) == type_float) {
    prec2 = 4;
    dt2 = 1;
  }
  else if(typeid(Type2) == type_double) {
    prec2 = 8;
    dt2 = 1;
  }
  else if(typeid(Type2) == type_complex) {
    prec2 = 4;
    dt2 = 2;
  }
  else if(typeid(Type2) == type_complex_double) {
    prec2 = 8;
    dt2 = 2;
  }
  if(prec != prec2)
    cout << "Error in trans_type1D: precisions don't match!" << endl;

  exec = plan.exec;
  doplan = plan.doplan;
  }
*/

trans_type3D::trans_type3D(int types_IDs[3]) //,const char name_[]) 
{
  //  name = new char[sizeof(name_)+1];
  //strcpy(name,name_);
  
  //  types = new trans_type1D[3];
  for(int i=0; i < 3; i++) {
    if(types1D.begin() + types_IDs[i] > types1D.end()) {
      cout << "Error in init trans_plan: type %d is undefined\n" << types_IDs[i] << endl;
      return;
    }
    //    trans_type1D *p=new trans_type1D(*types1D[types_IDs[i]]);
    types[i] = types_IDs[i];
  }
  dt1 = types1D[types[0]]->dt1;
  dt2 = types1D[types[2]]->dt2;
  prec = types1D[types[0]]->prec;
  is_set = true;
}

trans_type3D::trans_type3D(gen_trans_type *types_[3]) //st char name_[]) 
{
  //name = new char[sizeof(name_)+1];
  //rcpy(name,name_);
  for(int i=0; i < 3; i++) {
    int j=0;
    for(vector<gen_trans_type *>::iterator it=types1D.begin(); it < types1D.end(); it++,j++)
      if(types_[i] == *it)
	break;
    if(j >= types1D.size())
      types1D.push_back(types_[i]);
    types[i] = j;
  }  
  dt1 = types1D[types[0]]->dt1;
  dt2 = types1D[types[2]]->dt2;
  prec = types1D[types[0]]->prec;
  is_set = true;

}

trans_type3D::trans_type3D(const trans_type3D &tr) {
  //  name = new char[sizeof(tr.name)+1];
  //strcpy(name,tr.name);
  for(int i=0;i<3;i++) 
    types[i] = tr.types[i];

  dt1 = tr.dt1;
  dt2 = tr.dt2;
  is_set = true;
  prec = tr.prec;
} 


trans_type3D::~trans_type3D()
{
  //delete [] name;
}


template <class Type1,class Type2>  trans_type1D<Type1,Type2>::trans_type1D(const char *name_,long (*doplan_)(...),void (*exec_)(...),int isign) : gen_trans_type(name_,isign) {
  int prec2;


#ifdef FFTW
  if(typeid(Type1) ==typeid(float)) {
    dt1 = 1;
    prec = 4;
    if(typeid(Type2) ==type_float) {
	dt2 = 1;
	prec2 = 4;
	//	doplan = fftwf_plan_many_r2r;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...)) fftwf_execute_r2r;
    }
    else if(typeid(Type2) ==type_complex) {
	dt2 = 2;
	prec2 = 4;
	//doplan = fftwf_plan_many_r2c;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...)) fftwf_execute_dft_r2c;
      }
  }
  else if(typeid(Type1) ==type_double) {
      dt1 = 1;
      prec = 8;
      if(typeid(Type2) ==type_double) {
	dt2 = 1;
	prec2 = 8;
	//	doplan = fftw_plan_many_r2r;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftw_execute_r2r;
      }
      else if(typeid(Type2) ==type_complex_double) {
	dt2 = 2;
	prec2 = 8;
	// doplan = fftw_plan_many_r2c;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftw_execute_dft_r2c;
      }
  }
  else if(typeid(Type1) ==type_complex) {
      dt1 = 2;
      prec = 4;
      if(typeid(Type2) ==type_float) {
	dt2 = 1;
	prec2 = 4;
	//doplan = fftwf_plan_many_c2r;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftwf_execute_dft_c2r;
      }
      else if(typeid(Type2) ==type_complex) {
	dt2 = 2;
	prec2 = 4;
	//doplan = fftwf_plan_many;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftwf_execute_dft;
      }
  }
  else if(typeid(Type1) ==type_complex_double) {

      dt1 = 2;
      prec = 8;
      if(typeid(Type2) ==type_double) {
	dt2 = 1;
	prec2 = 8;
	//doplan = fftw_plan_many_c2r;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftw_execute_dft_c2r;
      }
      else if(typeid(Type2) ==type_complex_double) {
	dt2 = 2;
	prec2 = 8;
	//doplan = fftw_plan_many;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))fftw_execute_dft;
      }
    }
#endif    

  doplan = doplan_;

    if(!doplan)
      cout << "Error in trans_type1D: no suitable doplan" << endl;    
    if(prec != prec2)
      cout << "Error in trans_type1D: precisions don't match!" << endl;
  }

}

