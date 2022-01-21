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

#elif defined ESSL
#include "essl.h"
const int DEF_FFT_FLAGS=0;
extern "C" int enotrm(int &,int &);
extern "C" typedef int (*FN) (int &,int &);
#elif defined CUDA

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cufftXt.h>
#include <device_launch_parameters.h>
#include <curand_mtgp32_kernel.h>
#ifdef CUTENSOR
#include <cutensor.h>
#elif defined CUTT
#include <cutt.h>
#endif

const int DEF_FFT_FLAGS=0;

#endif

#define LocDevice 1
#define LocHost 2
#define NONE 0
#define SLICE 1
#define FULL 2

size_t inline  MULT3(int *X);
size_t inline  MULT3(int *X)
{size_t sz=1; for(int i=0;i<3;i++) sz *= X[i]; 
  return(sz);
} 

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

#ifdef CUDA
#include <helper_cuda.h>
#include <helper_functions.h>
#endif

#ifdef MKL_BLAS
#define MKL_Complex8 mycomplex
#define MKL_Complex16 complex_double
#include <mkl.h>
#endif

namespace p3dfft {

using namespace std;

static const int REAL=1;
static const int COMPLEX=2;

static const int TRANS_ONLY=1;
static const int MPI_ONLY=2;
static const int TRANSMPI=3;
const int CACHEPAD=32768;
  const int CACHE_BL=32768;
const int VECBLOCK=128;
//static int ls;

extern int EMPTY_TYPE_SINGLE,EMPTY_TYPE_DOUBLE,EMPTY_TYPE_SINGLE_COMPLEX,EMPTY_TYPE_DOUBLE_COMPLEX;
extern int R2CFFT_S,R2CFFT_D,C2RFFT_S,C2RFFT_D,CFFT_FORWARD_S,CFFT_FORWARD_D,CFFT_BACKWARD_S,CFFT_BACKWARD_D;
extern int DCT1_REAL_S,DCT1_REAL_D,DST1_REAL_S,DST1_REAL_D,DCT1_COMPLEX_S,DCT1_COMPLEX_D,DST1_COMPLEX_S,DST1_COMPLEX_D;
extern int DCT2_REAL_S,DCT2_REAL_D,DST2_REAL_S,DST2_REAL_D,DCT2_COMPLEX_S,DCT2_COMPLEX_D,DST2_COMPLEX_S,DST2_COMPLEX_D;
extern int DCT3_REAL_S,DCT3_REAL_D,DST3_REAL_S,DST3_REAL_D,DCT3_COMPLEX_S,DCT3_COMPLEX_D,DST3_COMPLEX_S,DST3_COMPLEX_D;
extern int DCT4_REAL_S,DCT4_REAL_D,DST4_REAL_S,DST4_REAL_D,DCT4_COMPLEX_S,DCT4_COMPLEX_D,DST4_COMPLEX_S,DST4_COMPLEX_D;

 extern int nslices;

#ifndef CUDA
 typedef int event_t;
#endif
//#else
//#include <complex>
typedef complex<float> mycomplex;
typedef complex<double> complex_double;


#define type_float typeid(float)
#define type_double typeid(double)
#define type_complex typeid(mycomplex)
#define type_complex_double typeid(complex_double)

void setup(int nslices=1);
void cleanup();

int arcmp(int *A,int *B,int N);
int excl(int,int),dist(int);
void divide_work(size_t *offset,size_t *mysize,int dims[3],int nslices,int split_dim);
void divide_dims(int **offset,int **mysize,int dims[3],int nslices,int split_dim);
  void compl_mo(int mo[3]);
bool cmpmo(int mo[3],int rhs);
void rel_change(int *,int *,int *);
void inv_mo(int mo[3],int imo[3]);
size_t max_long(size_t a,size_t b);
 int ar3d_cnt(int init,int pack_procs,int sz,int l,int od);
  template <class Type>	void  pack_ar(Type *in,Type *out,int ardims[3],int kst,int ken,int pack_dim,int pack_procs);
int swap0(int new_mo[3],int mo[3],int L,int *next=NULL);

#ifdef CUDA

 typedef cufftResult planResult;
 typedef cufftResult execResult;
 typedef cufftHandle planHandle;
 planResult plan_r2c_s(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_r2c_d(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2r_s(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2r_d(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2c_s(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2c_d(cufftHandle *plan, int *N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);

#elif defined ESSL

  struct essl_plan {
    double *aux1,*aux2;
    int naux1,naux2;
    int N,m,istride,idist,ostride,odist,isign;
  } ;
  
 typedef void planResult;
 typedef void execResult;
 typedef struct essl_plan planHandle;
 planResult plan_r2c_s(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_r2c_d(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2r_s(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2r_d(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist);
 planResult plan_c2c_s(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist, int isign);
  planResult plan_c2c_d(planHandle *plan, int N, int batch, int *inembed,int istride,int idist,int *onembed,int ostride,int odist, int isign);

#elif defined FFTW

typedef void execResult;
typedef long planResult;
typedef long planHandle;
#define plan_r2c_s fftwf_plan_many_dft_r2c
#define plan_r2c_d fftw_plan_many_dft_r2c
#define plan_c2r_s fftwf_plan_many_dft_c2r
#define plan_c2r_d fftw_plan_many_dft_c2r
#define plan_c2c_s fftwf_plan_many_dft
#define plan_c2c_d fftw_plan_many_dft
planResult plan_dct1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


planResult plan_dct2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


planResult plan_dct3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);


planResult plan_dct4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dct4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);
planResult plan_dst4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag = DEF_FFT_FLAGS);

#endif

void scheb_r(planHandle,float *,float *);
void dcheb_r(planHandle,double *,double *);
void scheb_c(planHandle,mycomplex *,mycomplex *);
void dcheb_c(planHandle,complex_double *,complex_double *);

 execResult exec_r2c_s(planHandle,float *,mycomplex *);
 execResult exec_r2c_d(planHandle,double *,complex_double *);
 execResult exec_c2r_s(planHandle,mycomplex *,float *);
 execResult exec_c2r_d(planHandle,complex_double *,double *);
 execResult exec_r2r_s(planHandle,float *,float *);
 execResult exec_r2r_d(planHandle,double *,double *);
 execResult exec_r2r_complex_s(planHandle,float *,float *);
 execResult exec_r2r_complex_d(planHandle,double *,double *);

#if defined FFTW || defined ESSL
#define exec_c2c_forward_s exec_c2c_s
#define exec_c2c_backward_s exec_c2c_s
#define exec_c2c_forward_d exec_c2c_d
#define exec_c2c_backward_d exec_c2c_d
execResult exec_c2c_s(planHandle,mycomplex *,mycomplex *);
execResult exec_c2c_d(planHandle,complex_double *,complex_double *);
#else
execResult exec_c2c_forward_s(planHandle,mycomplex *,mycomplex *);
execResult exec_c2c_forward_d(planHandle,complex_double *,complex_double *);
execResult exec_c2c_backward_d(planHandle,complex_double *,complex_double *);
execResult exec_c2c_backward_s(planHandle,mycomplex *,mycomplex *);
#endif

 template <class Type> void blas_trans(size_t rows,size_t cols,const double alpha,const Type *A,size_t lda,Type *B,size_t ldb);

class gen_trans_type {
 public :
  int isign;
  bool is_set,is_empty;
  int dt1,dt2;   //Datatype before and after
  int prec;
  inline gen_trans_type(int isign_=0)
  {
    //    name = new char[strlen(name_)+1];
    //    strcpy(name,name_);
    is_empty=false;
    is_set = true;
    isign = isign_;
  }
  inline gen_trans_type(int dt1_, int dt2_,int prec_,int isign_=0)
  {
    //    name = new char[strlen(name_)+1];
    //strcpy(name,name_);
    is_set = true;
    is_empty=false;
    isign = isign_;
    dt1 = dt1_;
    dt2 = dt2_;
    prec = prec_;
  }
  ~gen_trans_type() {}
  bool operator==(const gen_trans_type &) const;
  inline bool operator==(const gen_trans_type &rhs) {
    if(rhs.isign == isign && rhs.is_set == is_set && rhs.dt1 == dt1 &&
       rhs.dt2 == dt2)
      return(true);
    else
      return(false);
    
  }

};



template <class Type1,class Type2>  class trans_type1D  : public gen_trans_type{

  int ID;
  char *name;
  public :

  //typedef long (*doplan_type)(const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,...);

 planResult (*doplan)(...);

 // void (*exec)(long,Type1 *,Type2 *);
 void (*exec)(...);

 //trans_type1D(const char *name,  planResult (*doplan_)(...),void (*exec)(...)=NULL,int isign=0);
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
  
// Define a 1D transform type
//template <class Type1,class Type2>  trans_type1D<Type1,Type2>::
 trans_type1D(const char *name_,planResult (*doplan_)(...)=NULL,execResult (*exec_)(...)=NULL,bool is_empty=false,int isign=0) : gen_trans_type(isign) {
  int prec2;

    name = new char[strlen(name_)+1];
    strcpy(name,name_);


  if(typeid(Type1) ==typeid(float)) {
    dt1 = 1;
    prec = 4;
    if(typeid(Type2) ==type_float) {
	dt2 = 1;
	prec2 = 4;
	//	doplan = fftwf_plan_many_r2r;
	if(exec_) exec = (void (*)(...)) exec_;
#ifdef FFTW
	else exec = (void (*)(...)) fftwf_execute_r2r;
#elif defined CUDA
	else if(!is_empty) printf("Error in trans_type1D: CUDA does not support R2R transforms\n");
#endif
    }
    else if(typeid(Type2) ==type_complex) {
	dt2 = 2;
	prec2 = 4;
	//doplan = fftwf_plan_many_r2c;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...)) exec_r2c_s;
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
#ifdef FFTW
	else exec = (void (*)(...))fftw_execute_r2r;

#elif defined CUDA
	else if(!is_empty) printf("Error in trans_type1D: CUFFT does not support R2R transforms\n");
#endif
      }
      else if(typeid(Type2) ==type_complex_double) {
	dt2 = 2;
	prec2 = 8;
	// doplan = fftw_plan_many_r2c;
	if(exec_) exec = (void (*)(...)) exec_;
	else exec = (void (*)(...))exec_r2c_d;
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
	else exec = (void (*)(...))exec_c2r_s;
      }
      else if(typeid(Type2) ==type_complex) {
	dt2 = 2;
	prec2 = 4;
	//doplan = fftwf_plan_many;
	if(exec_) exec = (void (*)(...)) exec_;
	else if(isign == -1)
	  exec = (void (*)(...)) exec_c2c_forward_s;
	else
	  exec = (void (*)(...)) exec_c2c_backward_s;
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
	else exec = (void (*)(...)) exec_c2r_d;
      }
      else if(typeid(Type2) ==type_complex_double) {
	dt2 = 2;
	prec2 = 8;
	//doplan = fftw_plan_many;
	if(exec_) exec = (void (*)(...)) exec_;
	else if(isign == -1)
	  exec = (void (*)(...)) exec_c2c_forward_d;
	else
	  exec = (void (*)(...)) exec_c2c_backward_d;
      }
    }


  doplan = doplan_;

    if(!doplan && !is_empty)
      cout << "Error in trans_type1D: no suitable doplan" << endl;    
    if(prec != prec2)
      cout << "Error in trans_type1D: precisions don't match!" << endl;
  }

};



class ProcGrid;
class DataGrid;
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
  stage *next=NULL;
  int kind;
  void myexec(char *in,char *out,bool OW);
  void myexec_deriv(char *in,char *out, bool OW);

  stage() {}
  virtual ~stage() {}

};


template <class Type> class MPIplan : public stage {
  //    frien class transform3D;
    //    friend class stage;
 protected:
  int prec;
  int numtasks,taskid;
  DataGrid *grid1,*grid2;
  ProcGrid *Pgrid;
  //int mpicomm_ind;
  int *SndCnts,*RcvCnts,*SndStrt,*RcvStrt;
  int d1,d2,du; // Dimension ranks: from local, to local, and unchanged
  int comm_id; //which communicator is this? 0=row, 1=column etc
  int mo1[3],mo2[3];

  void pack_sendbuf(Type *sendbuf,Type *in);
  void unpack_recvbuf(Type *out,Type * recvbuf);
  bool is_set;

 public :

  size_t WorkSpace;

  MPIplan(const DataGrid &gr1,const DataGrid &gr2,int d1,int d2, int prec_);
  //MPIplan() {};
  ~MPIplan();
  void exec(char *in,char *out,char *tmpbuf=NULL);
  template <class Type1,class Type2> friend class trans_MPIplan;
  template <class Type11,class Type2> friend  class transform3D;

  };

 template <class Type1,class Type2>   class trans_MPIplan;


template <class Type1,class Type2>   class transplan : public stage {

 protected:

  int prec;
  DataGrid *grid1,*grid2;
  ProcGrid *Pgrid;
  int N,m,istride,idist,ostride,odist,isign;
  int *inembed,*onembed;
  unsigned fft_flag;
  void compute_deriv_loc(Type2 *in,Type2 *out,int dims[3]);
  //  char *DevBuf,*DevBuf2;

  public :

  size_t WorkSpace=0;

#ifdef CUDA
  //  bool useCuda=false;
  int InLoc=LocHost;
  int OutLoc=LocHost;
  cudaEvent_t EVENT_EXEC,EVENT_H2D,EVENT_D2H;
#endif
  size_t *offset1,*mysize1;
  size_t *offset2,*mysize2;
  bool is_empty=false;
  int trans_dim;  // Rank of dimension of the transform
  int mo1[3],mo2[3];
  bool is_set;
  trans_type1D<Type1,Type2> *trans_type;
  Plantype<Type1,Type2> *plan;
  transplan(const DataGrid &gr1,const DataGrid &gr2,const gen_trans_type *type,int d,int InLoc_=LocHost,int OutLoc_=LocHost); //, bool inplace_);
  transplan(const DataGrid &gr1,const DataGrid &gr2,int type_ID,int d,int InLoc_=LocHost,int OutLoc_=LocHost); //, bool inplace_); 
  void init_tr(const DataGrid &gr1,const DataGrid &gr2, const gen_trans_type *type,int d); // bool inplace_) 
  transplan() {};
  ~transplan() {
    delete grid1,grid2; 
    delete [] offset1,mysize1;
    delete [] offset2,mysize2;
#ifdef CUDA
    cudaEventDestroy(EVENT_EXEC);
    cudaEventDestroy(EVENT_H2D);
    cudaEventDestroy(EVENT_D2H);
#endif
  };

  void reorder_in_slice(Type1 *in,int mo1[3],int mo2[3],int dims1[3], int slice=0,int nslices=1,int pack_dim=-1,int pack_procs=0);
  
  void rot102in_slice(Type1 *in,Type2 *out,bool inplace,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot120in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice=0,int nslices=1,bool deriv=false, char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot210in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot201in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot021_op_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),Plantype<Type1,Type2> *plan, int cache_bl, int slice=0,int nslices=1,bool deriv=false,int pack_dim=-1,int pack_procs=0);

  void rot102out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot120out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot210out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot201out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice=0,int nslices=1,bool deriv=false,char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);

  void rot021_ip(Type1 *in,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl, bool deriv=false,char *tmpbuf=NULL);


#ifdef CUBLAS
    void rot102in_cublas(int lda,int ldb,Type2 *A,Type2 *C);
#endif

#ifdef CUDA
    void ro120in_cu(dim3 gridDim,dim3 blockSize,Type2 *in,Type2 *out,int d[3]);
    void ro120out_cu(dim3 gridDim,dim3 blockSize,Type1 *in,Type1 *out,int d[3]);
#endif

    void exec(char *in,char *out, int dim_deriv,bool OW=false, char *tmpbuf=NULL);
#ifdef CUDA
  void exec_slice(char *in,char *out, int dim_deriv,int slice,int nslices,event_t *event_hold,bool OW=false, char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);
#else
  void exec_slice(char *in,char *out, int dim_deriv,int slice,int nslices,int *event_hold,bool OW=false, char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);
#endif
  //  void reorder_in(Type1 *in,int mo1[3],int mo2[3],int dims1[3]);
  //  void reorder_in(Type1 *in,int mo1[3],int mo2[3],int dims1[3],char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);
  //void reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int dims1[3],int pack_dim=-1,int pack_procs=0);
  //void reorder_out_slice(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int dims1[3],int slice,int nslices,int pack_dim=-1,int pack_procs=0);
  // void reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW, char *tmpbuf=NULL);
  int reorder_trans_slice(Type1 *in,Type2 *out,int *mo1,int *mo2,void (*exec)(...),int dim_deriv,int slice=0,int nslices=1,event_t *event_hold=NULL,bool OW=false, char *tmpbuf=NULL,int pack_dim=-1,int pack_procs=0);
  //  void reorder_deriv(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW);
  void find_plan(trans_type1D<Type1,Type2> *type);
  //void exec_deriv(char *in,char *out, bool OW=false);
  int find_m(int *mo1,int *mo2,int *dims1,int *dims2,int trans_dim);

  //template <class Type1,class Type2>  
  friend class trans_MPIplan<Type1,Type2>;
  template <class Type11,class Type22> friend  class transform3D;
  
  };

template <class Type1,class Type2>   class trans_MPIplan : public stage {
  

 private : 

  int ***sndst,***sndsz,***snden;
  int ***rcvst,***rcvsz,***rcven;
  int **SndCnts,**SndStrt,**RcvCnts,**RcvStrt;   
  int *offset_snd,*offset_rcv;
  
  void pack_sendbuf_trans_slice(Type2 *sendbuf,char *in,int dim_deriv,event_t *event_hold,int slice=0,int nslices=1,char *devbuf=NULL,bool OW=false,char *tmpbuf=NULL);
  void pack_sendbuf_trans(Type2 *sendbuf,char *in,int dim_deriv,event_t *event_hold,int slice=0,int nslices=1,char *devbuf=NULL,bool OW=false,char *tmpbuf=NULL);
  void pack_sendbuf_trans(Type2 *sendbuf,char *in,int dim_deriv,bool OW);
  //  void pack_sendbuf_deriv(Type2 *sendbuf,char *in, bool OW);
  //  void unpack_recv(Type2 *out,Type2 * recvbuf);
  

  public :

  bool is_set;
  transplan<Type1,Type2>* trplan;
  MPIplan<Type2>* mpiplan;
  size_t WorkSpaceHost=0,WorkSpaceDev=0;
  
  trans_MPIplan(const DataGrid &gr1,const DataGrid &intergrid,const DataGrid &gr2,int d1,int d2,const gen_trans_type *type,int trans_dim_,int InLoc_); //,bool inplace_);
  ~trans_MPIplan() {
    int i,j;
    for(j=0;j<3;j++ )
      for(i=0;i<nslices;i++) {
	delete [] sndst[j][i];
	delete [] sndsz[j][i];
	delete [] snden[j][i];
	delete [] rcvst[j][i];
	delete [] rcvsz[j][i];
	delete [] rcven[j][i];
      }
    for(int i=0;i<nslices;i++) {
      delete [] SndCnts[i];
      delete [] SndStrt[i];
      delete [] RcvCnts[i];
      delete [] RcvStrt[i];
    }
    delete [] SndCnts,SndStrt,RcvCnts,RcvStrt,offset_snd,offset_rcv;
    for(j=0;j<3;j++) {
      delete [] sndst[j];
      delete [] snden[j];
      delete [] sndsz[j];
      delete [] rcvst[j];
      delete [] rcven[j];
      delete [] rcvsz[j];
    }
    delete [] sndst,sndsz,snden,rcvst,rcvsz,rcven;
  };
#ifdef CUDA
  int InLoc=LocHost;
#endif
  void exec(char *in,char *out,  int dim_deriv,event_t *event_hold,bool OW=false,char *tmpbuf=NULL,char *devbuf=NULL,double *tmpi=NULL);
  void exec_nb(char *in,char *out,  int dim_deriv,event_t *event_hold,bool OW=false,char *tmpbuf=NULL,char *devbuf=NULL,double *tmpi=NULL);
  void unpack_recvbuf_slice(Type2 *out,Type2 * recvbuf,int slice=0,int nslices=1);
#ifdef P2P
  void unpack_recvbuf_slice_p2p(Type2 *out,Type2 * recvbuf,int rank,int slice=0,int nslices=1,int **rcvstrt=NULL);
#endif
  //  void exec_deriv(char *in,char *out, bool OW);

  template <class TypeIn1,class TypeOut1> friend class transplan;
  template <class Type> friend class MPIplan;

  };

template <class Type>  void write_buf(Type *buf,char *label,int sz[3],int mo[3]);

 class ProcGrid {

  bool is_set;
 public:
  int taskid,numtasks;
  int nd;  //number of dimensions the volume is split over
  //  int proc_order[3];   //Ordering of tasks in processor grid, e.g. (1,2,3) : first dimension - adjacent tasks,then second, then third dimension
  int ProcDims[3];  //Processor grid size (in inverse order of split dimensions, i.e. rows first, then columns etc
  int grid_id_cart[3];
  MPI_Comm mpi_comm_glob; // Global MPi communicator we are starting from
  MPI_Comm mpi_comm_cart;
  MPI_Comm mpicomm[3]; //MPI communicators for each dimension 

  ProcGrid(int procdims[3],MPI_Comm mpi_comm_init);
  ProcGrid(const ProcGrid &rhs);
  ~ProcGrid();
  bool operator==( const ProcGrid &P1) const;
  inline bool operator==( const ProcGrid &P) {
    //    if(!(this->is_set))      return(false);
    //if (!(P.is_set))     return(false);

    if(this->nd != P.nd)
      return(false);
    if(this->taskid != P.taskid || this->numtasks != P.numtasks)
      return(false);
    int res;
    MPI_Comm_compare(this->mpi_comm_glob,P.mpi_comm_glob,&res);
    if(res != MPI_IDENT && res != MPI_CONGRUENT)
	return(false);
    MPI_Comm_compare(this->mpi_comm_cart,P.mpi_comm_cart,&res);
    if(res != MPI_IDENT && res != MPI_CONGRUENT)
	return(false);
    for(int i=0;i<3;i++) {
      if(this->ProcDims[i] != P.ProcDims[i])
	return(false);
      if(this->grid_id_cart[i] != P.grid_id_cart[i])
	return(false);
    }
    return(true);

  }

  friend class DataGrid;
  template<class Type>  friend  class MPIplan;
  template<class Type1,class Type2>  friend  class trans_plan;
  template<class Type1,class Type2>  friend  class trans_MPIplan;
};  


class DataGrid {
 private :

  //  void InitPencil2D(int myid[2]);
  void InitPencil();
  int *st[3],*sz[3],*en[3];  // Lowest, size and uppermost location in 3D, for each processor in subcommunicator

 public :
  //  int datatype;  //Datatype of each element (1 - real, 2 - complex)
  //  int prec;  //Precision (4/8)
  int nd;
  int Gdims[3];  //Global dimensions
  int dim_conj_sym; // Dimension of conjugate symmetry, where we store N/2+1 of the data after R2C transform due to conjugate symmety; =-1 if none
  int MemOrder[3];  //Memory ordering inside the data volume
  int Ldims[3];  //Local dimensions on THIS processor
  ProcGrid *Pgrid;
  int Pdims[3];  //Processor grid dimensions as mapped onto data dimensions
  int Dmap[3]; // How this data grid maps onto processor grid dimensions 
  int L[3];  //Rank of Local dimension (p=1)
  int D[3]; //Ranks of Dimensions of physical grid split over rows and columns correspondingly
  int GlobStart[3]; // Starting coords of this cube in the global grid
  int grid_id[3];  //Position of this pencil/cube in the processor grid
  //  int (*st)[3],(*sz)[3],(*en)[3];  // Lowest, size and uppermost location in 3D, for each processor in subcommunicator
  bool is_set;
  //bool IsLocal(int);
  inline  bool IsLocal(int dim) 
{
  if(dim < 0 || dim > 2)
    return(false);
  return(Pdims[dim] == 1);
    
};
  DataGrid(int *gdims_,int dim_conj_sym_,ProcGrid *pgrid,int *dmap,int *mem_order);
  DataGrid(const DataGrid &rhs);
  DataGrid() {};
  ~DataGrid();
  inline void set_gdims(int gdims[3]) {
    for(int i=0;i<3;i++) {
      Gdims[i] = gdims[i];
      //      ldims[i] = div_proc(gdims[i],pgrid[i],grid_id[i]);
      //      ldims[i] = gdims[i] / pgrid[i];
    }
    InitPencil();
  };


  void get_gdims(int gdims[3]) {for(int i=0;i<3;i++) gdims[i] = Gdims[i];};

  void set_mo(int mo[3]) {for(int i=0;i<3;i++) MemOrder[i] = mo[i];};

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
  planHandle *libplan_in,*libplan_out,*libplan_inout;
  Plan() {};
  ~Plan() {};
  };

template <class Type1,class Type2> class Plantype : public Plan
{
  int dt1,dt2;
  int prec;
  int N,m,d1,d2;
  //  bool inplace;
  int istride,idist,ostride,odist;
  int *inembed,*onembed;
  int isign,fft_flag;
  int typeID1,typeID2;
protected:
  int *mysize;

 public:
  planResult (*doplan)(...);
		 //int rank,const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,...);
  //  void (*exec)(long,Type1 *in,Type2 *out);
  void (*exec)(...);

  //int rank,const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,

  inline Plantype(planResult (*doplan_)(...),void (*exec_)(...),int N_,int m_,int istride_,int idist_,int ostride_,int odist_,int *inembed_=NULL,int *onembed_=NULL,int isign_=0,unsigned fft_flag_=DEF_FFT_FLAGS,int d1_=1, int d2_=1) 
  { doplan = doplan_;
    exec = exec_; 
    N = N_;m=m_;istride = istride_;istride = istride_;idist = idist_;
    ostride = ostride_;odist = odist_;isign = isign_;fft_flag = fft_flag_;
    inembed = inembed_;
    onembed = onembed_;
    d1 = d1_;
    d2 = d2_;
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
  friend void cleanup();
};


class trans_type3D {
 public :
  char *name;
  //  int dt1,dt2;
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

  bool find_order(int L[3],const trans_type3D *tp,const DataGrid *gr1,const DataGrid *gr2,bool *return_steps);

/*
class variable { // Do we need this?
  int ext,datatype,prec;
  void *data;
  grid *decomp;<Type1,Type2>
  bool alloc;
  variable(const variable &rhs);
  variable(const DataGrid &gr);
  variable~();
};
*/

  stage *init_transplan(const DataGrid &gr1,const DataGrid &gr2,const gen_trans_type *type,int d,int prec);
  stage *init_MPIplan(const DataGrid &gr1,const DataGrid &gr2,int d1,int d2, int dt,int prec);
  stage *init_trans_MPIplan(const DataGrid &gr1,const DataGrid &gr2,int d1,int d2, const gen_trans_type *type,int d, int prec);

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
  DataGrid *grid1,*grid2;
  ProcGrid *Pgrid;
  //  int trans_type[7];
  friend class stage;
#ifdef CUDA
  int InLoc,OutLoc;
#endif
  //  size_t *offset,*mysize;

 public:

  size_t WorkSpaceHost=0,WorkSpaceDev=0;
  void exec(Type1 *in,Type2 *out, bool OW=false, char *work_host=NULL, char *work_dev=NULL);
  void exec_deriv(Type1 *in,Type2 *out,int idir, bool OW=false, char *work_host=NULL, char *work_dev=NULL);

  transform3D(const DataGrid &grid1_, const DataGrid &grid2_,const trans_type3D *type,const int InLoc_=0, const int OutLoc_=0);
  ~transform3D();
};

 template <class Type> stage *final_seq(const DataGrid &grid1, const DataGrid &grid2, size_t *workspace, int loc1=0, int loc2=0);
 template <class Type> DataGrid *final_trans(DataGrid *grid1, const DataGrid &grid2, stage *curr,int prec,size_t *workspace);
//extern int ntrans;
//extern int npl;
//static Ntrans_types=0;
//const int MAX_TYPES=50;

//trans_type1D types1D[MAX_TYPES];


//transform3D stored_trans3D[MAX_NS];

extern int padd;
const int gblock=1;

//template <class Type> void reorder_out(Type *in,Type *out,int mo1[3],int mo2[3],int dims1[3]);
//template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int dims1[3]);

//template <class Type> void compute_deriv_loc(Type *in,Type *out,int dims[3],bool r2c); 
template <class Type> void compute_deriv(Type *in,Type *out,DataGrid *gr,int idir); 
// template <class Type> void compute_deriv(Type *in,Type *out,grid *gr,int idir);


extern vector<Plan *> Plans;
extern vector<gen_trans_type *> types1D;
extern vector<gen_transform3D *> stored_trans3D;
extern vector<stage *> stored_trans1D;
extern vector<trans_type3D> types3D;
extern vector<ProcGrid *> stored_proc_grids;
extern vector<DataGrid *> stored_data_grids;

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
    double gpu_transfer;
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

  //#include "templ.C"
  //#include "exec.C"
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
