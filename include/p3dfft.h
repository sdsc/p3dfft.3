#ifndef P3DFFT_3_H
#define P3DFFT_3_H
/*
Title: P3DFFT++ library

Authors: Dmitry Pekurovsky

Copyright (c) 2006-2019 

The Regents of the University of California.

All Rights Reserved.                        

 

    Permission to use, copy, modify and  distribute  any part

    of this software for  educational,  research  and  non-profit

    purposes, by individuals or non-profit organizations,

    without fee,  and  without a written  agreement is

    hereby granted,  provided  that the  above  copyright notice,

    this paragraph  and the following  three  paragraphs appear in

    all copies.       

 

    For-profit organizations desiring to use this software and others

    wishing to incorporate this  software into commercial

    products or use it for  commercial  purposes should contact the:    

          Office of Innovation & Commercialization 

          University of California San Diego

          9500 Gilman Drive,  La Jolla,  California, 92093-0910        

          Phone: (858) 534-5815

          E-mail: innovation@ucsd.edu

 

    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE

    TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR    

    CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT

    OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF

    CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH

    DAMAGE.

 

    THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND

    THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE        

    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 

    THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND    

    EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR

    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES

    OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR

    THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,        

    TRADEMARK OR OTHER RIGHTS.
*/

#ifdef HAVE_CONFIG_H
#include "p3dfft3config.h"
#endif

#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>


#ifdef FFTW
#include "fftw3.h"
#if defined FFTW_FLAG_MEASURE
const int DEF_FFT_FLAGS=FFTW_MEASURE;
#elif defined FFTW_FLAG_ESTIMATE
const int DEF_FFT_FLAGS=FFTW_ESTIMATE;
#elif defined FFTW_FLAG_PATIENT
const int DEF_FFT_FLAGS=FFTW_PATIENT;
#endif
#endif


extern int P3DFFT_EMPTY_TYPE_SINGLE,P3DFFT_EMPTY_TYPE_DOUBLE,P3DFFT_EMPTY_TYPE_SINGLE_COMPLEX,P3DFFT_EMPTY_TYPE_DOUBLE_COMPLEX;
extern int P3DFFT_R2CFFT_S,P3DFFT_R2CFFT_D,P3DFFT_C2RFFT_S,P3DFFT_C2RFFT_D,P3DFFT_CFFT_FORWARD_S,P3DFFT_CFFT_FORWARD_D,P3DFFT_CFFT_BACKWARD_S,P3DFFT_CFFT_BACKWARD_D;

extern int P3DFFT_DCT1_REAL_S,P3DFFT_DCT1_REAL_D,P3DFFT_DST1_REAL_S,P3DFFT_DST1_REAL_D,P3DFFT_DCT1_COMPLEX_S,P3DFFT_DCT1_COMPLEX_D,P3DFFT_DST1_COMPLEX_S,P3DFFT_DST1_COMPLEX_D;
extern int P3DFFT_DCT2_REAL_S,P3DFFT_DCT2_REAL_D,P3DFFT_DST2_REAL_S,P3DFFT_DST2_REAL_D,P3DFFT_DCT2_COMPLEX_S,P3DFFT_DCT2_COMPLEX_D,P3DFFT_DST2_COMPLEX_S,P3DFFT_DST2_COMPLEX_D;
extern int P3DFFT_DCT3_REAL_S,P3DFFT_DCT3_REAL_D,P3DFFT_DST3_REAL_S,P3DFFT_DST3_REAL_D,P3DFFT_DCT3_COMPLEX_S,P3DFFT_DCT3_COMPLEX_D,P3DFFT_DST3_COMPLEX_S,P3DFFT_DST3_COMPLEX_D;
extern int P3DFFT_DCT4_REAL_S,P3DFFT_DCT4_REAL_D,P3DFFT_DST4_REAL_S,P3DFFT_DST4_REAL_D,P3DFFT_DCT4_COMPLEX_S,P3DFFT_DCT4_COMPLEX_D,P3DFFT_DST4_COMPLEX_S,P3DFFT_DST4_COMPLEX_D;
extern int P3DFFT_DCT1_REAL_S,P3DFFT_DCT1_REAL_D,P3DFFT_DST1_REAL_S,P3DFFT_DST1_REAL_D,P3DFFT_DCT1_COMPLEX_S,P3DFFT_DCT1_COMPLEX_D,P3DFFT_DST1_COMPLEX_S,P3DFFT_DST1_COMPLEX_D;
extern int P3DFFT_DCT2_REAL_S,P3DFFT_DCT2_REAL_D,P3DFFT_DST2_REAL_S,P3DFFT_DST2_REAL_D,P3DFFT_DCT2_COMPLEX_S,P3DFFT_DCT2_COMPLEX_D,P3DFFT_DST2_COMPLEX_S,P3DFFT_DST2_COMPLEX_D;
extern int P3DFFT_DCT3_REAL_S,P3DFFT_DCT3_REAL_D,P3DFFT_DST3_REAL_S,P3DFFT_DST3_REAL_D,P3DFFT_DCT3_COMPLEX_S,P3DFFT_DCT3_COMPLEX_D,P3DFFT_DST3_COMPLEX_S,P3DFFT_DST3_COMPLEX_D;
extern int P3DFFT_DCT4_REAL_S,P3DFFT_DCT4_REAL_D,P3DFFT_DST4_REAL_S,P3DFFT_DST4_REAL_D,P3DFFT_DCT4_COMPLEX_S,P3DFFT_DCT4_COMPLEX_D,P3DFFT_DST4_COMPLEX_S,P3DFFT_DST4_COMPLEX_D;
extern int P3DFFT_DCT1_REAL_S,P3DFFT_DCT1_REAL_D,P3DFFT_DST1_REAL_S,P3DFFT_DST1_REAL_D,P3DFFT_DCT1_COMPLEX_S,P3DFFT_DCT1_COMPLEX_D,P3DFFT_DST1_COMPLEX_S,P3DFFT_DST1_COMPLEX_D;
extern int P3DFFT_DCT2_REAL_S,P3DFFT_DCT2_REAL_D,P3DFFT_DST2_REAL_S,P3DFFT_DST2_REAL_D,P3DFFT_DCT2_COMPLEX_S,P3DFFT_DCT2_COMPLEX_D,P3DFFT_DST2_COMPLEX_S,P3DFFT_DST2_COMPLEX_D;
extern int P3DFFT_DCT3_REAL_S,P3DFFT_DCT3_REAL_D,P3DFFT_DST3_REAL_S,P3DFFT_DST3_REAL_D,P3DFFT_DCT3_COMPLEX_S,P3DFFT_DCT3_COMPLEX_D,P3DFFT_DST3_COMPLEX_S,P3DFFT_DST3_COMPLEX_D;
extern int P3DFFT_DCT4_REAL_S,P3DFFT_DCT4_REAL_D,P3DFFT_DST4_REAL_S,P3DFFT_DST4_REAL_D,P3DFFT_DCT4_COMPLEX_S,P3DFFT_DCT4_COMPLEX_D,P3DFFT_DST4_COMPLEX_S,P3DFFT_DST4_COMPLEX_D;
extern int P3DFFT_DCT1_REAL_S,P3DFFT_DCT1_REAL_D,P3DFFT_DST1_REAL_S,P3DFFT_DST1_REAL_D,P3DFFT_DCT1_COMPLEX_S,P3DFFT_DCT1_COMPLEX_D,P3DFFT_DST1_COMPLEX_S,P3DFFT_DST1_COMPLEX_D;
extern int P3DFFT_DCT2_REAL_S,P3DFFT_DCT2_REAL_D,P3DFFT_DST2_REAL_S,P3DFFT_DST2_REAL_D,P3DFFT_DCT2_COMPLEX_S,P3DFFT_DCT2_COMPLEX_D,P3DFFT_DST2_COMPLEX_S,P3DFFT_DST2_COMPLEX_D;
extern int P3DFFT_DCT3_REAL_S,P3DFFT_DCT3_REAL_D,P3DFFT_DST3_REAL_S,P3DFFT_DST3_REAL_D,P3DFFT_DCT3_COMPLEX_S,P3DFFT_DCT3_COMPLEX_D,P3DFFT_DST3_COMPLEX_S,P3DFFT_DST3_COMPLEX_D;
extern int P3DFFT_DCT4_REAL_S,P3DFFT_DCT4_REAL_D,P3DFFT_DST4_REAL_S,P3DFFT_DST4_REAL_D,P3DFFT_DCT4_COMPLEX_S,P3DFFT_DCT4_COMPLEX_D,P3DFFT_DST4_COMPLEX_S,P3DFFT_DST4_COMPLEX_D;
  //P3DFFT_CHEB_REAL_S,P3DFFT_CHEB_REAL_D,P3DFFT_CHEB_COMPLEX_S,P3DFFT_CHEB_COMPLEX_D;

//#define Grid int
#define Type3D int
#define Type1D int
#define Plan3D int
#define Plan1D int

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <typeinfo>
#include <complex>


namespace p3dfft {

using namespace std;

int arcmp(int *A,int *B,int N);

static int ls;

static const int REAL=1;
static const int COMPLEX=2;

static const int TRANS_ONLY=1;
static const int MPI_ONLY=2;
static const int TRANSMPI=3;

extern int EMPTY_TYPE_SINGLE,EMPTY_TYPE_DOUBLE,EMPTY_TYPE_SINGLE_COMPLEX,EMPTY_TYPE_DOUBLE_COMPLEX;
 extern int R2CFFT_S,R2CFFT_D,C2RFFT_S,C2RFFT_D,CFFT_FORWARD_S,CFFT_FORWARD_D,CFFT_BACKWARD_S,CFFT_BACKWARD_D;
 extern int DCT1_REAL_S,DCT1_REAL_D,DST1_REAL_S,DST1_REAL_D,DCT1_COMPLEX_S,DCT1_COMPLEX_D,DST1_COMPLEX_S,DST1_COMPLEX_D;
 extern int DCT2_REAL_S,DCT2_REAL_D,DST2_REAL_S,DST2_REAL_D,DCT2_COMPLEX_S,DCT2_COMPLEX_D,DST2_COMPLEX_S,DST2_COMPLEX_D;
 extern int DCT3_REAL_S,DCT3_REAL_D,DST3_REAL_S,DST3_REAL_D,DCT3_COMPLEX_S,DCT3_COMPLEX_D,DST3_COMPLEX_S,DST3_COMPLEX_D;
 extern int DCT4_REAL_S,DCT4_REAL_D,DST4_REAL_S,DST4_REAL_D,DCT4_COMPLEX_S,DCT4_COMPLEX_D,DST4_COMPLEX_S,DST4_COMPLEX_D;

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

#ifdef MKL_BLAS
#define MKL_Complex8 mycomplex
#define MKL_Complex16 complex_double
#include <mkl.h>
#endif

#define type_float typeid(float)
#define type_double typeid(double)
#define type_complex typeid(mycomplex)
#define type_complex_double typeid(complex_double)

void setup();
void cleanup();

const int CACHEPAD=32768;
const int CACHE_BL=32768;
 const int VECBLOCK=128;

void rel_change(int *,int *,int *);
void inv_mo(int mo[3],int imo[3]);

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
void exec_r2r_s(long,float *,float *);
void exec_r2r_d(long,double *,double *);
void exec_r2r_complex_s(long,float *,float *);
void exec_r2r_complex_d(long,double *,double *);

 template <class Type> void blas_trans(size_t rows,size_t cols,const double alpha,const Type *A,size_t lda,Type *B,size_t ldb);

class gen_trans_type {
 public :
  char *name;
  int isign;
  bool is_set,is_empty=false;
  int dt1,dt2;   //Datatype before and after
  int prec;
  inline gen_trans_type(const char *name_,int isign_=0)
  {
    name = new char[strlen(name_)+1];
    strcpy(name,name_);
    is_set = true;
    isign = isign_;
  }
  inline gen_trans_type(const char *name_,int dt1_, int dt2_,int prec_,int isign_=0)
  {
    name = new char[strlen(name_)+1];
    strcpy(name,name_);
    is_set = true;
    isign = isign_;
    dt1 = dt1_;
    dt2 = dt2_;
    prec = prec_;
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
long plan_dct1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


long plan_dct2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


long plan_dct3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


long plan_dct4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dct4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
long plan_dst4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);

#endif

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
      ID = -1;
      dt1 = rhs.dt1;
      dt2 = rhs.dt2;
      prec = rhs.prec;
      name = new char[sizeof(rhs.name)];
      strcpy(name,rhs.name);
      //    splan = rhs.splan;
      //splan = rhs.dplan;
      //sexec = rhs.sexec;
      //dexec = rhs.dexec;
      isign = rhs.isign;
      is_set = rhs.is_set;
      is_empty = rhs.is_empty;
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
  //  int mo1[3],mo2[3];
 public :
  int stage_prec;
  //  bool inplace;
  int dt1,dt2;
  int dims1[3],dims2[3];
  //  bool is_set;
  //  bool is_trans = false;
  //bool is_mpi = false;
  stage *next;
  int kind;
  stage() {next = NULL;};
  void myexec(char *in,char *out,bool OW);
  void myexec_deriv(char *in,char *out, bool OW);
};


template <class Type> class MPIplan : public stage {
  //    frien class transform3D;
    //    friend class stage;
 protected:
  int prec;
  int numtasks,taskid;
  grid *grid1,*grid2;
  MPI_Comm mpicomm;
  int *SndCnts,*RcvCnts,*SndStrt,*RcvStrt;
  int d1,d2,du; // Dimension ranks: from local, to local, and unchanged
  int comm_id; //which communicator is this? 0=row, 1=column etc
  int mo1[3],mo2[3];

  void pack_sendbuf(Type *sendbuf,Type *in);
  void unpack_recvbuf(Type *out,Type * recvbuf);
  bool is_set;

 public :

  MPIplan(const grid &gr1,const grid &gr2,MPI_Comm comm,int d1,int d2, int prec_);
  MPIplan() {};
  ~MPIplan();
  void exec(char *in,char *out);
  template <class Type1,class Type2> friend class trans_MPIplan;
  };

 template <class Type1,class Type2>   class trans_MPIplan;

template <class Type1,class Type2>   class transplan : public stage {

 protected:

  int prec;
#ifdef ESSL
  double *work1,double *work2;
#endif
  grid *grid1,*grid2;
  int N,m,istride,idist,ostride,odist,isign;
  int *inembed,*onembed;
  unsigned fft_flag;
  void compute_deriv_loc(Type2 *in,Type2 *out,int dims[3]);

  public :

  int trans_dim;  // Rank of dimension of the transform
  int mo1[3],mo2[3];
  bool is_set;
  trans_type1D<Type1,Type2> *trans_type;
  Plantype<Type1,Type2> *plan;
  transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d); //, bool inplace_);
  transplan(const grid &gr1,const grid &gr2,int type_ID,int d); //, bool inplace_); 
  transplan() {};
  ~transplan() {delete grid1,grid2;};
  //  void reorder_in(Type1 *in,int mo1[3],int mo2[3],int dims1[3]);
  void reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int dims1[3]);
  void reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW);
  void reorder_deriv(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW);
  void find_plan(trans_type1D<Type1,Type2> *type);
  void exec(char *in,char *out,bool OW=false);
  void exec_deriv(char *in,char *out, bool OW=false);
  int find_m(int *mo1,int *mo2,int *dims1,int *dims2,int trans_dim);

  //template <class Type1,class Type2>  
  friend class trans_MPIplan<Type1,Type2>;

  };

template <class Type1,class Type2>   class trans_MPIplan : public stage {
  

 private : 

  void pack_sendbuf_trans(Type2 *sendbuf,char *in,bool OW);
  void pack_sendbuf_deriv(Type2 *sendbuf,char *in, bool OW);
  //  void unpack_recv(Type2 *out,Type2 * recvbuf);
  

  public :

  bool is_set;
  transplan<Type1,Type2>* trplan;
  MPIplan<Type2>* mpiplan;

  trans_MPIplan(const grid &gr1,const grid &intergrid,const grid &gr2,MPI_Comm mpicomm,int d1,int d2,const gen_trans_type *type,int trans_dim_); //,bool inplace_);
  ~trans_MPIplan();
  void exec(char *in,char *out, bool OW);
  void exec_deriv(char *in,char *out, bool OW);

  template <class TypeIn1,class TypeOut1> friend class transplan;
  template <class Type> friend class MPIplan;

  };

template <class Type>  void write_buf(Type *buf,char *label,int sz[3],int mo[3]);

class grid {
 private :

  //  void InitPencil2D(int myid[2]);
  void InitPencil();
  int pm; // Maximum of the three processor dimensions


 public :
  //  int datatype;  //Datatype of each element (1 - real, 2 - complex)
  //  int prec;  //Precision (4/8)
  int taskid,numtasks;
  int nd;  //number of dimensions the volume is split over
  int gdims[3];  //Global dimensions
  int dim_conj_sym; // Dimension of conjugate symmetry, where we store N/2+1 of the data after R2C transform due to conjugate symmety; =-1 if none
  int mem_order[3];  //Memory ordering inside the data volume
  int ldims[3];  //Local dimensions on THIS processor
  int pgrid[3];  //Processor grid
  int proc_order[3];   //Ordering of tasks in processor grid, e.g. (1,2,3) : first dimension - adjacent tasks,then second, then third dimension
  int P[3];  //Processor grid size (in inverse order of split dimensions, i.e. rows first, then columns etc
  int D[3];  //Ranks of Dimensions of physical grid split over rows and columns correspondingly
  int L[3];  //Rank of Local dimension (p=1)
  int grid_id[3];  //Position of this pencil/cube in the processor grid
  int grid_id_cart[3];
  int glob_start[3]; // Starting coords of this cube in the global grid
  MPI_Comm mpi_comm_glob; // Global MPi communicator we are starting from
  MPI_Comm mpi_comm_cart;
  MPI_Comm mpicomm[3]; //MPI communicators for each dimension 
  //  int (*st)[3],(*sz)[3],(*en)[3];  // Lowest, size and uppermost location in 3D, for each processor in subcommunicator
  int **st[3],**sz[3],**en[3];  // Lowest, size and uppermost location in 3D, for each processor in subcommunicator
  bool is_set;
  grid(int gdims_[3],int dim_conj_sym_,int pgrid_[3],int proc_order_[3],int mem_order[3],
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
  long libplan_in,libplan_out,libplan_inout;
  Plan() {};
  };

template <class Type1,class Type2> class Plantype : public Plan
{
  int dt1,dt2;
  int prec;
  int N,m;
  //  bool inplace;
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

  inline Plantype(long (*doplan_)(...),void (*exec_)(...),int N_,int m_,int istride_,int idist_,int ostride_,int odist_,int *inembed_=NULL,int *onembed_=NULL,int isign_=0,unsigned fft_flag_=DEF_FFT_FLAGS) 
  { doplan = doplan_;
    exec = exec_; 
    N = N_;m=m_;istride = istride_;istride = istride_;idist = idist_;
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
  //template <class _TypeIn,class _TypeOut> 
  friend void transplan<Type1,Type2>::find_plan(trans_type1D<Type1,Type2> *type); 
 //  template <class Type1,class Type2> friend  long trans_MPIplan<Type1,Type2>::find_plan(trans_type1D<Type1,Type2> *type);

};


class trans_type3D {
 public :
  char *name;
  int dt1,dt2;
  int prec;

  //  friend class transform3D;

  bool is_set;
  //  trans_type1D **types;
  int types[3];
  trans_type3D(gen_trans_type *types_[3]); //,const char name_[]="");
  trans_type3D(int types_IDs[3]); // ,const char name_[]="");
  trans_type3D(const trans_type3D &rhs);
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

  stage *init_transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d,int prec);
  stage *init_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2, int dt,int prec);
  stage *init_trans_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2, const gen_trans_type *type,int d, int prec);

class gen_transform3D
{
 public:
  int prec; // Precision
  int dt1,dt2;   //Datatype before and after
  bool OW;
};

template<class Type1,class Type2> class transform3D : public gen_transform3D
{
  stage *Stages;
//  MPIplan MPI_plans[7];
//  trans_plan trans_plans[7];
//  int nstages=0;
  bool is_set;
  grid *grid1,*grid2;
  //  int trans_type[7];
  friend class stage;

 public:

  void exec(Type1 *in,Type2 *out, bool OW=false);
  void exec_deriv(Type1 *in,Type2 *out,int idir, bool OW=false);

  transform3D(const grid &grid1_, const grid &grid2_,const trans_type3D *type);
  ~transform3D();
};

//extern int ntrans;
//extern int npl;
//static Ntrans_types=0;
//const int MAX_TYPES=50;

//trans_type1D types1D[MAX_TYPES];


//transform3D stored_trans3D[MAX_NS];

extern int padd;
const int gblock=1;

template <class Type> void reorder_out(Type *in,Type *out,int mo1[3],int mo2[3],int dims1[3]);
template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int dims1[3]);

//template <class Type> void compute_deriv_loc(Type *in,Type *out,int dims[3],bool r2c); 
template <class Type> void compute_deriv(Type *in,Type *out,grid *gr,int idir); 
// template <class Type> void compute_deriv(Type *in,Type *out,grid *gr,int idir);


extern vector<Plan *> Plans;
extern vector<gen_trans_type *> types1D;
extern vector<gen_transform3D *> stored_trans3D;
extern vector<stage *> stored_trans1D;
extern vector<trans_type3D> types3D;
extern vector<grid> stored_grids;

template<class Type> gen_trans_type *empty_type();
template<class Type> gen_trans_type *empty_type()
{
  gen_trans_type *t;
  if(typeid(Type) == type_float)
    t = types1D[EMPTY_TYPE_SINGLE];
  else   if(typeid(Type) == type_double)
    t = types1D[EMPTY_TYPE_DOUBLE];
  if(typeid(Type) == type_complex)
    t = types1D[EMPTY_TYPE_SINGLE_COMPLEX];
  if(typeid(Type) == type_complex_double)
    t = types1D[EMPTY_TYPE_DOUBLE_COMPLEX];
  return t;
}

 //CHEB_REAL_S,CHEB_REAL_D,CHEB_COMPLEX_S,CHEB_COMPLEX_D;


template class transform3D<float,float>;
template class transform3D<double,double>;
template class transform3D<mycomplex,float>;
template class transform3D<complex_double,double>;
template class transform3D<float,mycomplex>;
template class transform3D<double,complex_double>;
template class transform3D<mycomplex,mycomplex>;
template class transform3D<complex_double,complex_double>;

/*
template class trans_type1D<float,float>;
template class trans_type1D<double,double>;
template class trans_type1D<mycomplex,float>;
template class trans_type1D<complex_double,double>;
template class trans_type1D<float,mycomplex>;
template class trans_type1D<double,complex_double>;
template class trans_type1D<mycomplex,mycomplex>;
template class trans_type1D<complex_double,complex_double>;
*/
template class transplan<float,float>;
template class transplan<double,double>;
template class transplan<mycomplex,float>;
template class transplan<complex_double,double>;
template class transplan<float,mycomplex>;
template class transplan<double,complex_double>;
template class transplan<mycomplex,mycomplex>;
template class transplan<complex_double,complex_double>;

#ifdef TIMERS
  class timer {
  protected:
    double reorder_deriv;
    double reorder_trans;
    double reorder_out;
    double reorder_in;
    double trans_exec;
    double trans_deriv;
    double packsend;
    double packsend_trans;
    double packsend_deriv;
    double unpackrecv;
    double alltoall;
  public:
    void init();
    void print(MPI_Comm);

    template<class Type1,class Type2> friend  class transform3D;
    template<class Type1,class Type2> friend  class transplan;
    //    friend void transplan<Type1,Type2>::exec(char *,char *);
    template<class Type1,class Type2> friend  class trans_MPIplan;
    template<class Type1> friend class  MPIplan;
    //    friend void transform3D.exec;
    //    friend void transplan::exec(char *,char *);
    //friend void trans_MPIplan.exec;
    //friend void MPIplan.exec;
    //friend class stage;

  };

  extern  timer timers;

  //Timers
#endif

  // Namespace
}


extern "C" {
#endif

#include "Cwrap.h"
#include "Fwrap.h"

  //#ifdef __cplusplus
  //}
  //#endif

  //#ifdef __cplusplus
  //extern "C" {
  //#endif

#ifdef __cplusplus
}
#endif


#endif
