//#ifndef P3DFFT
//#define P3DFFT
#include <iostream>
#include "mpi.h"
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <typeinfo>
#include <complex>
#ifdef FFTW
#include "fftw3.h"
const int DEF_FFT_FLAGS=FFTW_MEASURE;
#endif
using namespace std;

static int mytype,mpireal,mpicomplex;
const int MAX_NS=50;
static int ls;

static const int REAL=1;
static const int COMPLEX=2;

static const int TRANS_ONLY=1;
static const int MPI_ONLY=2;
static const int TRANSMPI=3;

extern int EMPTY_TYPE,R2CFFT_S,R2CFFT_D,C2RFFT_S,C2RFFT_D,CFFT_FORWARD_S,CFFT_FORWARD_D,CFFT_BACKWARD_S,CFFT_BACKWARD_D,COSTRAN_REAL_S,COSTRAN_REAL_D,SINTRAN_REAL_S,SINTRAN_REAL_D,COSTRAN_COMPLEX_S,COSTRAN_COMPLEX_D,SINTRAN_COMPLEX_S,SINTRAN_COMPLEX_D,CHEB_REAL_S,CHEB_REAL_D,CHEB_COMPLEX_S,CHEB_COMPLEX_D;
/*
const int R2CFFT_S=0;
const int R2CFFT_D=1;
const int C2RFFT_S=2;
const int C2RFFT_D=3;
const int CFFT_FORWARD_S=4;
const int CFFT_FORWARD_D=5;
const int CFFT_BACKWARD_S=6;
const int CFFT_BACKWARD_D=7;
const int COSTRAN_REAL_S=8;
const int COSTRAN_REAL_D=9;
const int COSTRAN_COMPLEX_S=10;
const int COSTRAN_COMPLEX_D=11;
const int SINTRAN_REAL_S=12;
const int SINTRAN_REAL_D=13;
const int SINTRAN_COMPLEX_S=14;
const int SINTRAN_COMPLEX_D=15;
//const int CHEB_REAL_S=8;
//const int CHEB_REAL_D=9;
//const int CHEB_COMPLEX_S=10;
//const int CHEB_COMPLEX_D=11;
*/


const int CACHEPAD=32768;
const int CACHE_BL=32768;
#ifdef FFTW
//typedef fftwf_complex complex;
//typedef complex<float> mycomplex;
//typedef fftw_complex complex_double;
typedef fftw_plan lib_plan_double_type;
typedef fftwf_plan lib_plan_type;
#endif
//#else
//#include <complex>
typedef complex<float> mycomplex;
typedef complex<double> complex_double;

#define type_float typeid(float)
#define type_double typeid(double)
#define type_complex typeid(mycomplex)
#define type_complex_double typeid(complex_double)

using namespace std;

void p3dfft_setup();
void p3dfft_clean();
void rel_change(int *,int *,int *);

void scheb_r(long,float *,float *);
void dcheb_r(long,double *,double *);
void scheb_c(long,mycomplex *,mycomplex *);
void dcheb_c(long,complex_double *,complex_double *);
void exec_r2c_s(long,float *,mycomplex *);
void exec_r2c_d(long,double *,complex_double *);
void exec_c2r_s(long,mycomplex *,float *);
void exec_c2r_d(long,complex_double *,double *);
void exec_c2c_s(long,mycomplex *,mycomplex *);
void exec_c2c_d(long,complex_double *,complex_double *);

class gen_trans_type {
 public :
  char *name;
  int isign;
  bool is_set,is_empty;
  int dt1,dt2;   //Datatype before and after
  int prec;
  inline gen_trans_type(const char *name_,int isign_=0)
  {
    name = new char[strlen(name_)+1];
    strcpy(name,name_);
    is_set = true;
    isign = isign_;
  }
  ~gen_trans_type() { delete [] name;}
  bool operator==(const gen_trans_type &) const;
  inline bool operator==(const gen_trans_type &rhs) {
    if(rhs.isign == isign && rhs.is_set == is_set && rhs.dt1 == dt1 &&
       rhs.dt2 == dt2)
      return(true);
    else
      return(false);
  }

};

#ifdef FFTW
#define plan_r2c_s fftwf_plan_many_dft_r2c
#define plan_r2c_d fftw_plan_many_dft_r2c
#define plan_c2r_s fftwf_plan_many_dft_c2r
#define plan_c2r_d fftw_plan_many_dft_c2r
#define plan_c2c_s fftwf_plan_many_dft
#define plan_c2c_d fftw_plan_many_dft
long plan_cos_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_cos_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_sin_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_sin_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);

#endif

/*
long plan_r2c_s(...);
long plan_r2c_d(...);
long plan_c2r_s(...);
long plan_c2r_d(...);
long plan_c2c_s(...);
long plan_c2c_d(...);
*/



template <class Type1,class Type2>  class trans_type1D  : public gen_trans_type{

  int ID;
  public :

typedef long (*doplan_type)(const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,...);

 long (*doplan)(...);

 // void (*exec)(long,Type1 *,Type2 *);
 void (*exec)(...);

  trans_type1D(const char *name,  long (*doplan_)(...),void (*exec)(...)=NULL,int isign=0);
  inline int getID() {return(ID);}

  trans_type1D(const trans_type1D &rhs) 
    {
      ID = 0;
      dt1 = rhs.dt1;
      dt2 = rhs.dt2;
      name = new char[sizeof(rhs.name)];
      strcpy(name,rhs.name);
      splan = rhs.splan;
      dplan = rhs.dplan;
      sexec = rhs.sexec;
      dexec = rhs.dexec;
      isign = rhs.isign;
      is_set = rhs.is_set;
    }
  
  ~trans_type1D()
    {
      delete [] name;
    }

};



class grid;
template <class Type1,class Type2> class Plantype;



class stage {
  //  friend class transform3D;
 protected :
  int prec;
  int mo1[3],mo2[3];
 public :
  bool inplace;
  int dt1,dt2;
  int dims1[3],dims2[3];
  bool is_set;
  //  bool is_trans = false;
  //bool is_mpi = false;
  stage *next;
  int kind;
  stage() {next = NULL;};
  void myexec(char *in,char *out);
};

template <class Type> class MPIplan  : public stage {
  //    friend class transform3D;
    //    friend class stage;
 protected:
  int numtasks,taskid,commid;
  grid *grid1,*grid2;
  MPI_Comm mpicomm;
  int *SndCnts,*RcvCnts,*SndStrt,*RcvStrt;
  int d1,d2,du; // Dimension ranks: from local, to local, and unchanged

  void pack_sendbuf(Type *sendbuf,Type *in);
  void unpack_recvbuf(Type *out,Type * recvbuf);

 public :

  MPIplan(const grid &gr1,const grid &gr2,MPI_Comm comm,int d1,int d2);
  MPIplan() {};
  ~MPIplan();
  void exec(char *in,char *out);
};

template <class Type1,class Type2>   class transplan : public stage { 

 protected:

  Plantype<Type1,Type2> *plan;
#ifdef ESSL
  double *work1,double *work2;
#endif
  grid *grid1,*grid2;
  int N,m,istride,idist,ostride,odist,isign;
  int *inembed,*onembed;
  unsigned fft_flag;

  int trans_dim;  // Rank of dimension of the transform
  //  trans_type1D<Type1,Type2> 



  public :

  long lib_plan;
  trans_type1D<Type1,Type2> *trans_type;
  transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d, bool inplace_);
  transplan(const grid &gr1,const grid &gr2,int type_ID,int d, bool inplace_); 
  transplan() {};
  ~transplan() {delete grid1,grid2;};
  void reorder_in(Type1 *in,int mo1[3],int mo2[3],int dims1[3]);
  void reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int dims1[3]);
  void reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1);
  long find_plan(trans_type1D<Type1,Type2> *type);
  void exec(char *in,char *out);
};

template <class Type1,class Type2>   class trans_MPIplan : public transplan<Type1,Type2>, public MPIplan<Type2> { 

 private : 

  void pack_sendbuf_trans(Type2 *sendbuf,char *in);
  //  void unpack_recv(Type2 *out,Type2 * recvbuf);

  public :

  trans_MPIplan(const grid &gr1,const grid &intergrid,const grid &gr2,MPI_Comm mpicomm,int d1,int d2,const gen_trans_type *type,int trans_dim_,bool inplace_);
  ~trans_MPIplan();
  void exec(char *in,char *out);
  };

class grid {
 private :

  void InitPencil2D(int myid[2]);
  void InitPencil1D();
  int pm; // Maximum of the three processor dimensions


 public :
  //  int datatype;  //Datatype of each element (1 - real, 2 - complex)
  //  int prec;  //Precision (4/8)
  int taskid;
  int nd;  //number of dimensions the volume is split over
  int gdims[3];  //Global dimensions
  int mem_order[3];  //Memory ordering inside the data volume
  int ldims[3];  //Local dimensions on THIS processor
  int pgrid[3];  //Processor grid
  int proc_order[3];   //Ordering of tasks in processor grid, e.g. (1,2,3) : first dimension - adjacent tasks,then second, then third dimension
  int P[3];  //Processor grid size (in inverse order of split dimensions, i.e. rows first, then columns etc
  int D[3];  //Ranks of Dimensions of physical grid split over rows and columns correspondingly
  int L[3];  //Rank of Local dimension (p=1)
  int grid_id[3];  //Position of this pencil/cube in the processor grid
  int glob_start[3]; // Starting coords of this cube in the global grid
  MPI_Comm mpi_comm_glob; // Global MPi communicator we are starting from
  MPI_Comm mpicomm[3]; //MPI communicators for each dimension 
  int (*st)[3],(*sz)[3],(*en)[3];  // Lowest, size and uppermost location in 3D, for each processor in subcommunicator
  bool is_set;
  grid(int gdims_[3],int pgrid_[3],int proc_order_[3],int mem_order[3],
	    MPI_Comm mpicomm_);
  grid(const grid &rhs);
  grid() {};
  ~grid();
  void set_mo(int mo[3]) {for(int i=0;i<3;i++) mem_order[i] = mo[i];};

template<class Type>  friend  class MPIplan;
 template<class Type1,class Type2>  friend  class trans_plan;
 template<class Type1,class Type2>  friend  class trans_MPIplan;
};  


class Plan {
  /*  
  lib_plan libplan;
  lib_plan_double libplan_double;
  */
  //  friend class trans_type1D;
 public :
  long libplan;
  Plan() {};
  };

template <class Type1,class Type2> class Plantype : public Plan
{
  int dt1,dt2;
  int prec;
  int N,m;
  bool inplace;
  int istride,idist,ostride,odist;
  int *inembed,*onembed;
  int isign,fft_flag;
  int typeID1,typeID2;

 public:
  long (*doplan)(...);
		 //int rank,const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,...);
  //  void (*exec)(long,Type1 *in,Type2 *out);
  void (*exec)(...);

  //int rank,const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,

  inline Plantype(long (*doplan_)(...),void (*exec_)(...),int N_,int m_,bool inplace_,int istride_,int idist_,int ostride_,int odist_,int *inembed_=NULL,int *onembed_=NULL,int isign_=0,unsigned fft_flag_=DEF_FFT_FLAGS) 
  { doplan = doplan_;
    exec = exec_; 
    N = N_;m=m_;inplace=inplace_;istride = istride_;istride = istride_;idist = idist_;
    ostride = ostride_;odist = odist_;isign = isign_;fft_flag = fft_flag_;
    inembed = inembed_;
    onembed = onembed_;
    int prec2;
    if(typeid(Type1) == typeid(float)) {
      dt1 = 1; prec=4;
    }
    else if(typeid(Type1) == typeid(double)) {
      dt1 = 1; prec = 8;
    }
    else if(typeid(Type1) == typeid(mycomplex)) {
      dt1 = 2; prec = 4;
    }
    else if(typeid(Type1) == typeid(complex_double)) {
      dt1 = 2; prec = 8;
    }
    if(typeid(Type2) == typeid(float)) {
      dt2 = 1; prec2=4;
    }
    else if(typeid(Type2) == typeid(double)) {
      dt2 = 1; prec2 = 8;
    }
    else if(typeid(Type2) == typeid(mycomplex)) {
      dt2 = 2; prec2 = 4;
    }
    else if(typeid(Type2) == typeid(complex_double)) {
      dt2 = 2; prec2 = 8;
    }
    if(prec != prec2) 
      cout << "Error in Plantype: precisions don't match" << endl;
  }
  inline Plantype(const Plantype &rhs);
  inline ~Plantype();
  template <class Type1,class Type2>  friend long transplan<class Type1,class Type2>::find_plan(trans_type1D<Type1,Type2> *type);
 //  template <class Type1,class Type2> friend  long trans_MPIplan<Type1,Type2>::find_plan(trans_type1D<Type1,Type2> *type);

};


class trans_type3D {
  char *name;
  int dt1,dt2;
  int prec;

  //  friend class transform3D;

  bool is_set;
  //  trans_type1D **types;
 protected :
  int types[3];
 public :
  trans_type3D(gen_trans_type *types_[3],const char name_[]);
  trans_type3D(int types_IDs[3],const char name_[]);
  ~trans_type3D();
 template <class Type1,class Type2>    friend class transform3D;
 friend int print_type3D(const trans_type3D *);
};

/*
class variable { // Do we need this?
  int ext,datatype,prec;
  void *data;
  grid *decomp;<Type1,Type2>
  bool alloc;
  variable(const variable &rhs);
  variable(const grid &gr);
  variable~();
};
*/

  stage *init_transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d, bool inplace,int prec);
  stage *init_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2, int dt,int prec);
  stage *init_trans_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2, const gen_trans_type *type,int d, bool inplace,int prec);

class gen_transform3D
{
};

template<class Type1,class Type2> class transform3D : public gen_transform3D
{
  stage *Stages;
//  MPIplan MPI_plans[7];
//  trans_plan trans_plans[7];
//  int nstages=0;
  bool inplace,is_set;
  grid grid1,grid2;
  //  int trans_type[7];
  int prec; // Precision
  int dt1,dt2;   //Datatype before and after
  friend class stage;

 public:

  void exec(Type1 *in,Type2 *out,int OW);

  transform3D(const grid &grid1_, const grid &grid2_,const trans_type3D *type,bool inplace_); 
  ~transform3D();
};

//extern int ntrans;
//extern int npl;
//static Ntrans_types=0;
//const int MAX_TYPES=50;

//trans_type1D types1D[MAX_TYPES];

extern vector<Plan*> Plans;
extern vector <gen_trans_type*> types1D;
extern vector <gen_transform3D*> stored_trans3D;

//transform3D stored_trans3D[MAX_NS];

extern int padd;
const int gblock=1;

template <class Type> void reorder_out(Type *in,Type *out,int mo1[3],int mo2[3],int dims1[3]);
template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int dims1[3]);


//#endif
