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

#include "p3dfft.h"
#include <string.h>

int P3DFFT_EMPTY_TYPE_SINGLE,P3DFFT_EMPTY_TYPE_DOUBLE,P3DFFT_EMPTY_TYPE_SINGLE_COMPLEX,P3DFFT_EMPTY_TYPE_DOUBLE_COMPLEX;
int P3DFFT_R2CFFT_S,P3DFFT_R2CFFT_D,P3DFFT_C2RFFT_S,P3DFFT_C2RFFT_D,P3DFFT_CFFT_FORWARD_S,P3DFFT_CFFT_FORWARD_D,P3DFFT_CFFT_BACKWARD_S,P3DFFT_CFFT_BACKWARD_D;
int P3DFFT_DCT1_REAL_S,P3DFFT_DCT1_REAL_D,P3DFFT_DST1_REAL_S,P3DFFT_DST1_REAL_D,P3DFFT_DCT1_COMPLEX_S,P3DFFT_DCT1_COMPLEX_D,P3DFFT_DST1_COMPLEX_S,P3DFFT_DST1_COMPLEX_D;
int P3DFFT_DCT2_REAL_S,P3DFFT_DCT2_REAL_D,P3DFFT_DST2_REAL_S,P3DFFT_DST2_REAL_D,P3DFFT_DCT2_COMPLEX_S,P3DFFT_DCT2_COMPLEX_D,P3DFFT_DST2_COMPLEX_S,P3DFFT_DST2_COMPLEX_D;
int P3DFFT_DCT3_REAL_S,P3DFFT_DCT3_REAL_D,P3DFFT_DST3_REAL_S,P3DFFT_DST3_REAL_D,P3DFFT_DCT3_COMPLEX_S,P3DFFT_DCT3_COMPLEX_D,P3DFFT_DST3_COMPLEX_S,P3DFFT_DST3_COMPLEX_D;
int P3DFFT_DCT4_REAL_S,P3DFFT_DCT4_REAL_D,P3DFFT_DST4_REAL_S,P3DFFT_DST4_REAL_D,P3DFFT_DCT4_COMPLEX_S,P3DFFT_DCT4_COMPLEX_D,P3DFFT_DST4_COMPLEX_S,P3DFFT_DST4_COMPLEX_D;

//P3DFFT_CHEB_REAL_S,P3DFFT_CHEB_REAL_D,P3DFFT_CHEB_COMPLEX_S,P3DFFT_CHEB_COMPLEX_D;

namespace p3dfft {

using namespace std;

  // External vectors 

vector<Plan *> Plans;  // Defined 1D transform plans
vector<gen_trans_type *> types1D; // Defined 1D transform types
vector<trans_type3D> types3D; // Defined 3D transform types
vector<stage *> stored_trans1D; // Initialized and planned 1D transforms
vector<gen_transform3D *> stored_trans3D; // Initialized and planned 3D transforms
vector<ProcGrid*> stored_proc_grids; // Defined proc. grids (used in Fortran and C wrappers)
vector<DataGrid*> stored_data_grids; // Defined data grids (used in Fortran and C wrappers)


  //  extern "C" {
int EMPTY_TYPE_SINGLE,EMPTY_TYPE_DOUBLE,EMPTY_TYPE_SINGLE_COMPLEX,EMPTY_TYPE_DOUBLE_COMPLEX;
  int R2CFFT_S,R2CFFT_D,C2RFFT_S,C2RFFT_D,CFFT_FORWARD_S,CFFT_FORWARD_D,CFFT_BACKWARD_S,CFFT_BACKWARD_D;
  int DCT1_REAL_S,DCT1_REAL_D,DST1_REAL_S,DST1_REAL_D,DCT1_COMPLEX_S,DCT1_COMPLEX_D,DST1_COMPLEX_S,DST1_COMPLEX_D;
  int DCT2_REAL_S,DCT2_REAL_D,DST2_REAL_S,DST2_REAL_D,DCT2_COMPLEX_S,DCT2_COMPLEX_D,DST2_COMPLEX_S,DST2_COMPLEX_D;
  int DCT3_REAL_S,DCT3_REAL_D,DST3_REAL_S,DST3_REAL_D,DCT3_COMPLEX_S,DCT3_COMPLEX_D,DST3_COMPLEX_S,DST3_COMPLEX_D;
  int DCT4_REAL_S,DCT4_REAL_D,DST4_REAL_S,DST4_REAL_D,DCT4_COMPLEX_S,DCT4_COMPLEX_D,DST4_COMPLEX_S,DST4_COMPLEX_D;
  //CHEB_REAL_S,CHEB_REAL_D,CHEB_COMPLEX_S,CHEB_COMPLEX_D;

  // }

int padd;

void setup(int nslices_)
//#else
//void setup()
//#endif
{
  const char *name;
  int isign;
  gen_trans_type *p;

  int types_count=0;
  /*
#ifdef CUDA
  TILE_DIM=32;
  SHMEM_SIZE=16384;
#endif
  */
  /*
#ifdef ESSL
 FN iusadr;
 int ierno,inoal,inomes,itrace,irange,irc,dummy;
 int naux;

 iusadr = enotrm;

 dummy = 0;
 einfo (0,dummy,dummy);

 ierno = 2015;
 inoal = 0;
 inomes = -1;
 itrace = 0;
 irange = 2015;
 errset (ierno,inoal,inomes,itrace,&iusadr,irange);
#endif
  */
#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Single" << endl;
#endif
  name = "Empty Type Single";
  p = new trans_type1D<float,float>(name,NULL,NULL,true);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_SINGLE = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Double" << endl;
#endif
  name = "Empty Type Double";
  p = new trans_type1D<double,double>(name,NULL,NULL,true);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_DOUBLE = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Single Complex" << endl;
#endif
  name = "Empty Type Single Complex";
  p = new trans_type1D<mycomplex,mycomplex>(name,NULL,NULL,true);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_SINGLE_COMPLEX = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty TypeDouble Complex" << endl;
#endif
  name = "Empty Type Double Complex";
  p = new trans_type1D<complex_double,complex_double>(name,NULL,NULL,true);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_DOUBLE_COMPLEX = types_count;
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

  p = new trans_type1D<float,mycomplex>(name,(planResult (*)(...) ) plan_r2c_s,(execResult (*)(...)) exec_r2c_s);
  types1D.push_back(p);
  R2CFFT_S =  types_count;
  types_count++;
  

#ifdef DEBUG
  cout << "p3dft_setup: adding R2C double type" << endl;
#endif

  name ="Real-to-complex Fourier Transform double precision";

  p = new trans_type1D<double,complex_double>(name,(planResult (*)(...) ) plan_r2c_d,(execResult (*)(...)) exec_r2c_d);

  types1D.push_back(p);
  R2CFFT_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2R single type" << endl;
#endif
  name = "Complex-to-real Fourier Transform, single precision";

  p = new trans_type1D<mycomplex,float>(name,(planResult (*)(...) ) plan_c2r_s,(execResult (*)(...)) exec_c2r_s);

  types1D.push_back(p);
  C2RFFT_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2R double type" << endl;
#endif

  name = "Complex-to-real Fourier Transform, double precision";

  p = new trans_type1D<complex_double,double>(name,(planResult (*)(...) ) plan_c2r_d,(execResult (*)(...)) exec_c2r_d);

  types1D.push_back(p);
  C2RFFT_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C forward single type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, forward transform, singple precision";

#ifdef FFTW
  isign = FFTW_FORWARD;
#elif defined ESSL
  isign = 1;
#else
  isign = 0;
#endif
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_c2c_s,(execResult (*)(...) )exec_c2c_forward_s,false,isign);

  types1D.push_back(p);
  CFFT_FORWARD_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C forward double type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, forward transform, double precision";

#ifdef FFTW
  isign = FFTW_FORWARD;
#elif defined ESSL
  isign = 1;
#else
  isign = 0;
#endif
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_c2c_d,(execResult (*)(...)) exec_c2c_forward_d,false,isign);


  types1D.push_back(p);
  CFFT_FORWARD_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C backward single type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, backward transform, single precision";

#ifdef FFTW
  isign = FFTW_BACKWARD;
#elif defined ESSL
  isign = -1;
#else
  isign = 0;
#endif
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_c2c_s,(execResult (*)(...)) exec_c2c_backward_s,false,isign);

  types1D.push_back(p);
  CFFT_BACKWARD_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding C2C double backward type" << endl;
#endif

  name = "Complex-to-complex Fourier Transform, backward transform, double precision";

#ifdef FFTW
  isign = FFTW_BACKWARD;
#elif defined ESSL
  isign = -1;
#else
  isign = 0;
#endif
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_c2c_d,(execResult (*)(...)) exec_c2c_backward_d,false,isign);

  types1D.push_back(p);
  CFFT_BACKWARD_D = types_count;
  types_count++;

  /*
  name = "Real-valued Chebyshev Transform, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_cos_s,(void (*)(...)) scheb_r);
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

#ifndef CUDA
  // CUFFT does not support real-to-real transforms

#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT1 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT1, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dct1_s, (execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT1_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT1 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT1, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dct1_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT1_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT1 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT1, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(planResult (*)(...) ) plan_dct1_complex_s, (execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT1_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT1 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT1, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dct1_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT1_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST1 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST1, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dst1_s,(execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST1_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST1 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST1, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dst1_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST1_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST1 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST1, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_dst1_complex_s,(execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST1_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST1 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST1, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dst1_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST1_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT2 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT2, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dct2_s, (execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT2_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT2 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT2, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dct2_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT2_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT2 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT2, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(planResult (*)(...) ) plan_dct2_complex_s, (execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT2_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT2 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT2, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dct2_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT2_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST2 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST2, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dst2_s,(execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST2_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST2 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST2, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dst2_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST2_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST2 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST2, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_dst2_complex_s,(execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST2_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST2 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST2, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dst2_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST2_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT3 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT3, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dct3_s, (execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT3_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT3 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT3, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dct3_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT3_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT3 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT3, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(planResult (*)(...) ) plan_dct3_complex_s, (execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT3_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT3 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT3, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dct3_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT3_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST3 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST3, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dst3_s,(execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST3_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST3 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST3, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dst3_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST3_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST3 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST3, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_dst3_complex_s,(execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST3_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST3 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST3, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dst3_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST3_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT4 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT4, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dct1_s, (execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT4_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT4 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT4, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dct1_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT4_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT4 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT4, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(planResult (*)(...) ) plan_dct1_complex_s, (execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT4_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT4 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT4, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dct1_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT4_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST4 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST4, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(planResult (*)(...) ) plan_dst4_s,(execResult (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST4_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST4 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST4, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(planResult (*)(...) ) plan_dst4_d,(execResult (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST4_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST4 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST4, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(planResult (*)(...) ) plan_dst4_complex_s,(execResult (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST4_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST4 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST4, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(planResult (*)(...) ) plan_dst4_complex_d,(execResult (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST4_COMPLEX_D = types_count;
  types_count++;
#endif

  nslices = nslices_;
#ifdef CUDA
  streams = new cudaStream_t[nslices];
  for (int i = 0; i < nslices; i++)
    {
      checkCudaErrors(cudaStreamCreate(&streams[i]));
    }
#endif

}



void cleanup()
{
  
  //  vector<grid>::iterator it1=stored_grids.begin();
  //  stored_grids.erase(stored_grids.begin(),stored_grids.end());
  //  for(vector<trans_type3D>::iterator it=types3D.begin();it != types3D.end();it++) {
  //    types3D.erase(types3D.begin(),types3D.end());
  stored_proc_grids.clear();
  stored_data_grids.clear();
  types3D.clear();

  // Since these are vectors of pointers, simply erasing is not enough; must delete each referenced class. 
  //  printf("Clearing types1D\n");
  vector<gen_trans_type *>::iterator it1=types1D.begin();
  while(it1 != types1D.end()) {
    if((*it1)->prec == 4) {
      if((*it1)->dt1 == 1 && (*it1)->dt2 == 1) {
	trans_type1D<float,float> *t = reinterpret_cast<trans_type1D<float,float>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 1 && (*it1)->dt2 == 2) {
	trans_type1D<float,mycomplex> *t = reinterpret_cast<trans_type1D<float,mycomplex>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 2 && (*it1)->dt2 == 1) {
	trans_type1D<mycomplex,float> *t = reinterpret_cast<trans_type1D<mycomplex,float>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 2 && (*it1)->dt2 == 2) {
	trans_type1D<mycomplex,mycomplex> *t = reinterpret_cast<trans_type1D<mycomplex,mycomplex>*>( *it1);
	delete t;
      }
    }
    else if((*it1)->prec == 8) {
      if((*it1)->dt1 == 1 && (*it1)->dt2 == 1) {
	trans_type1D<double,double> *t = reinterpret_cast<trans_type1D<double,double>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 1 && (*it1)->dt2 == 2) {
	trans_type1D<double,complex_double> *t = reinterpret_cast<trans_type1D<double,complex_double>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 2 && (*it1)->dt2 == 1) {
	trans_type1D<complex_double,double> *t = reinterpret_cast<trans_type1D<complex_double,double>*>( *it1);
	delete t;
      }
      else      if((*it1)->dt1 == 2 && (*it1)->dt2 == 2) {
	trans_type1D<complex_double,complex_double> *t = reinterpret_cast<trans_type1D<complex_double,complex_double>*>( *it1);
	delete t;
      }
    }
    
    it1 = types1D.erase(it1);
  }
  
  //  printf("Clearing Plans\n");
  //  for(vector<Plan *>::iterator it=Plans.begin();it != Plans.end();it++) {
  vector<Plan *>::iterator it=Plans.begin();
  while(it != Plans.end()) {
    for(int i=0;i<nslices;i++) {
#ifdef FFTW
      if((*it)->libplan_in[i] != NULL) 
	fftw_destroy_plan((fftw_plan) (*it)->libplan_in[i]);
      if((*it)->libplan_out[i] != NULL) 
	fftw_destroy_plan((fftw_plan) (*it)->libplan_out[i]);
      if((*it)->libplan_inout[i] != NULL) 
	fftw_destroy_plan((fftw_plan) (*it)->libplan_inout[i]);
#elif defined ESSL
      //      if((*it)->libplan_in[i] != NULL) 
      delete [] (*it)->libplan_in[i].aux1;
	//if((*it)->libplan_out[i] != NULL) 
      //	delete [] (*it)->libplan_out[i].aux1;
	// if((*it)->libplan_inout[i] != NULL) 
      //	delete [] (*it)->libplan_inout[i].aux1;
#elif defined CUDA
      if((*it)->libplan_in[i] != NULL) 
	cufftDestroy((cufftHandle) (*it)->libplan_in[i]);
      if((*it)->libplan_out[i] != NULL) 
	cufftDestroy((cufftHandle) (*it)->libplan_out[i]);
      if((*it)->libplan_inout[i] != NULL) 
	cufftDestroy((cufftHandle) (*it)->libplan_inout[i]);
#endif
    }
    delete [] (*it)->libplan_in,(*it)->libplan_out,(*it)->libplan_inout;


    /*
    if((*it)->prec == 4) {
      if((*it)->dt1 == 1 && (*it)->dt2 == 1) {
	Plantype<float,float> *t = reinterpret_cast<Plantype<float,float>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 1 && (*it)->dt2 == 2) {
	Plantype<float,mycomplex> *t = reinterpret_cast<Plantype<float,mycomplex>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 2 && (*it)->dt2 == 1) {
	Plantype<mycomplex,float> *t = reinterpret_cast<Plantype<mycomplex,float>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 2 && (*it)->dt2 == 2) {
	Plantype<mycomplex,mycomplex> *t = reinterpret_cast<Plantype<mycomplex,mycomplex>*>( *it);
	delete t;
      }
    }
    else if((*it)->prec == 8) {
      if((*it)->dt1 == 1 && (*it)->dt2 == 1) {
	Plantype<double,double> *t = reinterpret_cast<Plantype<double,double>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 1 && (*it)->dt2 == 2) {
	Plantype<double,complex_double> *t = reinterpret_cast<Plantype<double,complex_double>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 2 && (*it)->dt2 == 1) {
	Plantype<complex_double,double> *t = reinterpret_cast<Plantype<complex_double,double>*>( *it);
	delete t;
      }
      else      if((*it)->dt1 == 2 && (*it)->dt2 == 2) {
	Plantype<complex_double,complex_double> *t = reinterpret_cast<Plantype<complex_double,complex_double>*>( *it);
	delete t;
      }
    }
    */
    it = Plans.erase(it);
  }


  //  printf("Clearing stored_trans1D\n");
  vector<stage *>::iterator it2=stored_trans1D.begin();
  while(it2 != stored_trans1D.end()) {
    switch((*it2)->kind) {
    case 1:
    
      if((*it2)->stage_prec == 4) {
	if((*it2)->dt1 == 1 && (*it2)->dt2 == 1) {
	  transplan<float,float> *t = reinterpret_cast<transplan<float,float>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 1 && (*it2)->dt2 == 2) {
	  transplan<float,mycomplex> *t = reinterpret_cast<transplan<float,mycomplex>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 1) {
	  transplan<mycomplex,float> *t = reinterpret_cast<transplan<mycomplex,float>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 2) {
	  transplan<mycomplex,mycomplex> *t = reinterpret_cast<transplan<mycomplex,mycomplex>*>( *it2);
	  delete t;
	}
      }
      else if((*it2)->stage_prec == 8) {
	if((*it2)->dt1 == 1 && (*it2)->dt2 == 1) {
	  transplan<double,double> *t = reinterpret_cast<transplan<double,double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 1 && (*it2)->dt2 == 2) {
	  transplan<double,complex_double> *t = reinterpret_cast<transplan<double,complex_double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 1) {
	  transplan<complex_double,double> *t = reinterpret_cast<transplan<complex_double,double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 2) {
	  transplan<complex_double,complex_double> *t = reinterpret_cast<transplan<complex_double,complex_double>*>( *it2);
	  delete t;
	}
      }
      
      break;
    case2:
      
      if((*it2)->stage_prec == 4) {
	if((*it2)->dt1 == 1) {
	  MPIplan<float> *t = reinterpret_cast<MPIplan<float>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt2 == 2) {
	  MPIplan<mycomplex> *t = reinterpret_cast<MPIplan<mycomplex>*>( *it2);
	  delete t;
	}
      }
      else if((*it2)->stage_prec == 8) {
	if((*it2)->dt1 == 1) {
	  MPIplan<double> *t = reinterpret_cast<MPIplan<double> *>( *it2);
	  delete t;
	}
	else      if((*it2)->dt2 == 2) {
	  MPIplan<complex_double> *t = reinterpret_cast<MPIplan<complex_double>*>( *it2);
	  delete t;
	}
      }
      
      break;
    case 3:
      
      if((*it2)->stage_prec == 4) {
	if((*it2)->dt1 == 1 && (*it2)->dt2 == 1) {
	  trans_MPIplan<float,float> *t = reinterpret_cast<trans_MPIplan<float,float>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 1 && (*it2)->dt2 == 2) {
	  trans_MPIplan<float,mycomplex> *t = reinterpret_cast<trans_MPIplan<float,mycomplex>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 1) {
	  trans_MPIplan<mycomplex,float> *t = reinterpret_cast<trans_MPIplan<mycomplex,float>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 2) {
	  trans_MPIplan<mycomplex,mycomplex> *t = reinterpret_cast<trans_MPIplan<mycomplex,mycomplex>*>( *it2);
	  delete t;
	}
      }
      else if((*it2)->stage_prec == 8) {
	if((*it2)->dt1 == 1 && (*it2)->dt2 == 1) {
	  trans_MPIplan<double,double> *t = reinterpret_cast<trans_MPIplan<double,double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 1 && (*it2)->dt2 == 2) {
	  trans_MPIplan<double,complex_double> *t = reinterpret_cast<trans_MPIplan<double,complex_double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 1) {
	  trans_MPIplan<complex_double,double> *t = reinterpret_cast<trans_MPIplan<complex_double,double>*>( *it2);
	  delete t;
	}
	else      if((*it2)->dt1 == 2 && (*it2)->dt2 == 2) {
	  trans_MPIplan<complex_double,complex_double> *t = reinterpret_cast<trans_MPIplan<complex_double,complex_double>*>( *it2);
	  delete t;
	}
      }
      
      break;
    }
    
    it2 = stored_trans1D.erase(it2);
  }

  //  printf("Clearing stored_trans3D\n");
  vector<gen_transform3D *>::iterator it3=stored_trans3D.begin();
  while(it3 != stored_trans3D.end()) {

    if((*it3)->prec == 4) {
      if((*it3)->dt1 == 1 && (*it3)->dt2 == 1) {
	transform3D<float,float> *t = reinterpret_cast<transform3D<float,float>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 1 && (*it3)->dt2 == 2) {
	transform3D<float,mycomplex> *t = reinterpret_cast<transform3D<float,mycomplex>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 2 && (*it3)->dt2 == 1) {
	transform3D<mycomplex,float> *t = reinterpret_cast<transform3D<mycomplex,float>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 2 && (*it3)->dt2 == 2) {
	transform3D<mycomplex,mycomplex> *t = reinterpret_cast<transform3D<mycomplex,mycomplex>*>( *it3);
	delete t;
      }
    }
    else if((*it3)->prec == 8) {
      if((*it3)->dt1 == 1 && (*it3)->dt2 == 1) {
	transform3D<double,double> *t = reinterpret_cast<transform3D<double,double>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 1 && (*it3)->dt2 == 2) {
	transform3D<double,complex_double> *t = reinterpret_cast<transform3D<double,complex_double>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 2 && (*it3)->dt2 == 1) {
	transform3D<complex_double,double> *t = reinterpret_cast<transform3D<complex_double,double>*>( *it3);
	delete t;
      }
      else      if((*it3)->dt1 == 2 && (*it3)->dt2 == 2) {
	transform3D<complex_double,complex_double> *t = reinterpret_cast<transform3D<complex_double,complex_double>*>( *it3);
	delete t;
      }
    }

    it3 = stored_trans3D.erase(it3);
  }

#ifdef FFTW
  //  printf("Cleaning FFTW\n");
  fftw_cleanup();
#endif    
#ifdef CUDA
  for(int i=0;i<nslices;i++)
    cudaStreamDestroy(streams[i]);
  delete [] streams;
#endif

  //  printf("Done Cleaning\n");

}

  
/*
planResult plan_r2c_s(const int *n,int howmany,float *in,const int *inembed,int istride,int idist,mycomplex *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftwf_plan_many_dft_r2c(1,n,howmany,in,inembed,istride,idist,(fftwf_complex *) out,onembed,ostride,odist,fft_flag));
#endif
}

planResult plan_r2c_d(const int *n,int howmany,double *in,const int *inembed,int istride,int idist,complex_double *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftw_plan_many_dft_r2c(1,n,howmany,in,inembed,istride,idist,(fftw_complex *) out,onembed,ostride,odist,fft_flag));
#endif
}
planResult plan_c2r_s(const int *n,int howmany,mycomplex *in,const int *inembed,int istride,int idist,float *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftwf_plan_many_dft_c2r(1,n,howmany,(fftwf_complex *) in,inembed,istride,idist,out,onembed,ostride,odist,fft_flag));
#endif
}

planResult plan_c2r_d(const int *n,int howmany,complex_double *in,const int *inembed,int istride,int idist,double *out,const int *onembed,int ostride,int odist,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftw_plan_many_dft_c2r(1,n,howmany,(fftw_complex *)in,inembed,istride,idist,out,onembed,ostride,odist,fft_flag));
#endif
}
planResult plan_c2c_s(const int *n,int howmany,mycomplex *in,const int *inembed,int istride,int idist,mycomplex *out,const int *onembed,int ostride,int odist,int isign,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftwf_plan_many_dft(1,n,howmany,(fftwf_complex *) in,inembed,istride,idist,(fftwf_complex *) out,onembed,ostride,odist,isign,fft_flag));
#endif
}

planResult plan_c2c_d(const int *n,int howmany,complex_double *in,const int *inembed,int istride,int idist,complex_double *out,const int *onembed,int ostride,int odist,int isign,unsigned fft_flag=DEF_FFT_FLAGS)
{
#ifdef FFTW
  return((planResult) fftw_plan_many_dft(1,n,howmany,(fftw_complex *) in,inembed,istride,idist,(fftw_complex *) out,onembed,ostride,odist,isign,fft_flag));
#endif
}
*/


#ifdef CUDA
  planResult plan_r2c_s(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_R2C,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT R2C_s\n");
  return(res);
}
  planResult plan_r2c_d(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_D2Z,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT R2C_d\n");
  return(res);
}

  planResult plan_c2r_s(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_C2R,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT C2R_s\n");
  return(res);
}
  planResult plan_c2r_d(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_Z2D,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT C2R_d\n");
  return(res);
}

  planResult plan_c2c_s(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_C2C,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT c2c_s\n");
  return(res);
}
  planResult plan_c2c_d(cufftHandle *plan, int *N, int batch,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  planResult res;
  res = cufftPlanMany(plan,1,N,inembed,istride,idist,onembed,ostride,odist,CUFFT_Z2Z,batch);
  if(res != CUFFT_SUCCESS)
    printf("Error in CUFFT C2C_d\n");
  return(res);
}



  execResult exec_r2c_s(planHandle plan,float *in,mycomplex *out)
{
  cufftExecR2C(plan,(cufftReal *) in,(cufftComplex *) out);
}
execResult exec_r2c_d(planHandle plan,double *in,complex_double *out)
{
  cufftExecD2Z(plan,(cufftDoubleReal *) in,(cufftDoubleComplex *) out);
}
execResult exec_c2r_s(planHandle plan,mycomplex *in,float *out)
{
  cufftExecC2R(plan,(cufftComplex *) in,(cufftReal *) out);
}
execResult exec_c2r_d(planHandle plan,complex_double *in,double *out)
{
  cufftExecZ2D(plan,(cufftDoubleComplex *) in,(cufftDoubleReal *) out);
}
execResult exec_c2c_forward_s(planHandle plan,mycomplex *in,mycomplex *out)
{
  cufftExecC2C(plan,(cufftComplex *) in,(cufftComplex *) out,CUFFT_FORWARD);
}
execResult exec_c2c_backward_s(planHandle plan,mycomplex *in,mycomplex *out)
{
  cufftExecC2C(plan,(cufftComplex *) in,(cufftComplex *) out,CUFFT_INVERSE);
}
execResult exec_c2c_forward_d(planHandle plan,complex_double *in,complex_double *out)
{
  cufftExecZ2Z(plan,(cufftDoubleComplex *) in,(cufftDoubleComplex *) out,CUFFT_FORWARD);
}
execResult exec_c2c_backward_d(planHandle plan,complex_double *in,complex_double *out)
{
  cufftExecZ2Z(plan,(cufftDoubleComplex *) in,(cufftDoubleComplex *) out,CUFFT_INVERSE);
}
#endif

#ifdef ESSL

planResult plan_c2c_d(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist, int isign) 
{
  if(N <= 8192)
    plan->naux1 = 30000;
  else
    plan->naux1 = 30000 + 1.14 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = isign;
  plan->N = N;
  plan->m = m;
  //  printf("Calling DCFT plan: N=%d,m=%d,isign=%d,istride=%d,idist=%d\n",N,m,isign,istride,idist);
  dcft(1,NULL,istride,idist,NULL,ostride,odist,N,m,isign,1.0,plan->aux1,plan->naux1,NULL,0);
}

planResult plan_c2c_s(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist, int isign) 
{
  if(N <= 8192)
    plan->naux1 = 30000;
  else
    plan->naux1 = 30000 + 1.14 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = isign;
  plan->N = N;
  plan->m = m;
  scft(1,NULL,istride,idist,NULL,ostride,odist,N,m,isign,1.0,plan->aux1,plan->naux1,NULL,0);
}

planResult plan_r2c_d(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  if(N <= 16384)
    plan->naux1 = 35000;
  else
    plan->naux1 = 30000 + 0.82 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = 1;
  plan->N = N;
  plan->m = m;
  drcft(1,NULL,idist,NULL,odist,N,m,plan->isign,1.0,plan->aux1,plan->naux1,NULL,0);
}

planResult plan_r2c_s(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  if(N <= 16384)
    plan->naux1 = 35000;
  else
    plan->naux1 = 30000 +0.82 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = 1;
  plan->N = N;
  plan->m = m;
  srcft(1,NULL,idist,NULL,odist,N,m,plan->isign,1.0,plan->aux1,plan->naux1,NULL,0,NULL,0);
}

planResult plan_c2r_d(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  if(N <= 16384)
    plan->naux1 = 35000;
  else
    plan->naux1 = 30000 + 0.82 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = -1;
  plan->N = N;
  plan->m = m;
  dcrft(1,NULL,idist,NULL,odist,N,m,plan->isign,1.0,plan->aux1,plan->naux1,NULL,0);
}

planResult plan_c2r_s(planHandle *plan, int N, int m,int *inembed,int istride,int idist,int *onembed,int ostride,int odist) 
{
  if(N <= 16384)
    plan->naux1 = 35000;
  else
    plan->naux1 = 30000 +0.82 * N;

  plan->aux1 = new double[plan->naux1];
  plan->istride = istride;
  plan->idist = idist;
  plan->ostride = ostride;
  plan->odist = odist;
  plan->isign = -1;
  plan->N = N;
  plan->m = m;
  scrft(1,NULL,idist,NULL,odist,N,m,plan->isign,1.0,plan->aux1,plan->naux1,NULL,0,NULL,0);
}



execResult exec_c2c_d(planHandle plan,complex_double *in,complex_double *out)
{
  //  printf("Calling DCFT exec: N=%d,m=%d,isign=%d,istride=%d,idist=%d\n",plan.N,plan.m,plan.isign,plan.istride,plan.idist);
  dcft(0,in,plan.istride,plan.idist,out,plan.ostride,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0);
}
execResult exec_c2c_s(planHandle plan,mycomplex *in,mycomplex *out)
{
  scft(0,in,plan.istride,plan.idist,out,plan.ostride,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0);
}

execResult exec_r2c_d(planHandle plan,double *in,complex_double *out)
{
  drcft(0,in,plan.idist,out,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0);
}
execResult exec_r2c_s(planHandle plan,float *in,mycomplex *out)
{
  srcft(0,in,plan.idist,out,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0,NULL,0);
}

execResult exec_c2r_d(planHandle plan,complex_double *in,double *out)
{
  dcrft(0,in,plan.idist,out,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0);
}
execResult exec_c2r_s(planHandle plan,mycomplex *in,float *out)
{
  scrft(0,in,plan.idist,out,plan.odist,plan.N,plan.m,plan.isign,1.0,plan.aux1,plan.naux1,NULL,0,NULL,0);
}
#endif

#ifdef FFTW
void exec_r2c_s(planHandle plan,float *in,mycomplex *out)
{
  fftwf_execute_dft_r2c((fftwf_plan) plan,in,(fftwf_complex *) out);
}
void exec_r2c_d(planHandle plan,double *in,complex_double *out)
{
  fftw_execute_dft_r2c((fftw_plan) plan,in,(fftw_complex *) out);
}
void exec_c2r_s(planHandle plan,mycomplex *in,float *out)
{
  fftwf_execute_dft_c2r((fftwf_plan) plan,(fftwf_complex *) in,out);
}
void exec_c2r_d(planHandle plan,complex_double *in,double *out)
{
  fftw_execute_dft_c2r((fftw_plan) plan,(fftw_complex *) in, out);
}
void exec_c2c_s(planHandle plan,mycomplex *in,mycomplex *out)
{
  fftwf_execute_dft((fftwf_plan) plan,(fftwf_complex *) in,(fftwf_complex *) out);
}
void exec_c2c_d(planHandle plan,complex_double *in,complex_double *out)
{
  fftw_execute_dft((fftw_plan) plan,(fftw_complex *) in,(fftw_complex *) out);
}
void exec_r2r_s(planHandle plan,float *in,float *out)
{
  fftwf_execute_r2r((fftwf_plan) plan,in,out);
}
void exec_r2r_d(planHandle plan,double *in,double *out)
{
  fftw_execute_r2r((fftw_plan) plan,in,out);
}

void exec_r2r_complex_s(planHandle plan,float *in,float *out)
{
  fftwf_execute_r2r((fftwf_plan) plan,in,out);
  fftwf_execute_r2r((fftwf_plan) plan,in+1,out+1);
}
void exec_r2r_complex_d(planHandle plan,double *in,double *out)
{
  fftw_execute_r2r((fftw_plan) plan,in,out);
  fftw_execute_r2r((fftw_plan) plan,in+1,out+1);
}


planResult plan_dct1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT00;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


planResult plan_dct1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT00;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst1_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT00;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst1_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT00;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dct1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT00;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dct1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT00;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT00;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT00;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



planResult plan_dct2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT10;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


planResult plan_dct2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT10;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT10;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT10;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dct2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT10;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dct2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT10;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT10;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT10;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



planResult plan_dct3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT01;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


planResult plan_dct3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT01;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT01;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT01;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dct3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT01;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dct3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT01;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT01;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT01;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



planResult plan_dct4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT11;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


planResult plan_dct4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT11;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT11;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dst4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT11;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


planResult plan_dct4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT11;
   return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dct4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT11;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT11;
  return((planResult) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


planResult plan_dst4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT11;
  return((planResult) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}
#endif



#ifdef TIMERS

void timer::init()
{

  reorder_trans = 0.0;
  reorder_out = 0.0;
  reorder_in = 0.0;
  trans_exec = 0.0;
  packsend = 0.0;
  packsend_trans = 0.0;
  unpackrecv = 0.0;
  alltoall = 0.0;
  gpu_transfer = 0.0;
}

void timer::print(MPI_Comm comm)
{

  timer gtimers_avg,gtimers_min,gtimers_max;
  gtimers_avg.init();
  gtimers_min.init();
  gtimers_max.init();
  int nproc,taskid;
  MPI_Comm_rank(comm,&taskid);
  MPI_Comm_size(comm,&nproc);

  MPI_Reduce(&reorder_trans,&gtimers_avg.reorder_trans,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&reorder_trans,&gtimers_min.reorder_trans,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&reorder_trans,&gtimers_max.reorder_trans,1,MPI_DOUBLE,MPI_MAX,0,comm);

  MPI_Reduce(&reorder_out,&gtimers_avg.reorder_out,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&reorder_out,&gtimers_min.reorder_out,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&reorder_out,&gtimers_max.reorder_out,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&reorder_in,&gtimers_avg.reorder_in,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&reorder_in,&gtimers_min.reorder_in,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&reorder_in,&gtimers_max.reorder_in,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&trans_exec,&gtimers_avg.trans_exec,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&trans_exec,&gtimers_min.trans_exec,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&trans_exec,&gtimers_max.trans_exec,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&packsend_trans,&gtimers_avg.packsend_trans,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&packsend_trans,&gtimers_min.packsend_trans,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&packsend_trans,&gtimers_max.packsend_trans,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&packsend,&gtimers_avg.packsend,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&packsend,&gtimers_min.packsend,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&packsend,&gtimers_max.packsend,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&unpackrecv,&gtimers_avg.unpackrecv,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&unpackrecv,&gtimers_min.unpackrecv,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&unpackrecv,&gtimers_max.unpackrecv,1,MPI_DOUBLE,MPI_MAX,0,comm);
  MPI_Reduce(&alltoall,&gtimers_avg.alltoall,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&alltoall,&gtimers_min.alltoall,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&alltoall,&gtimers_max.alltoall,1,MPI_DOUBLE,MPI_MAX,0,comm);
#ifdef CUDA
  MPI_Reduce(&gpu_transfer,&gtimers_avg.gpu_transfer,1,MPI_DOUBLE,MPI_SUM,0,comm);
  MPI_Reduce(&gpu_transfer,&gtimers_min.gpu_transfer,1,MPI_DOUBLE,MPI_MIN,0,comm);
  MPI_Reduce(&gpu_transfer,&gtimers_max.gpu_transfer,1,MPI_DOUBLE,MPI_MAX,0,comm);
#endif

  if(taskid == 0) {
    printf("TIMERS (avg/min/max): \nReorder_trans (%lf %lf %lf)\n",gtimers_avg.reorder_trans/nproc,gtimers_min.reorder_trans,gtimers_max.reorder_trans);
    printf("Reorder_out   (%lf %lf %lf)\n",gtimers_avg.reorder_out/nproc,gtimers_min.reorder_out,gtimers_max.reorder_out);
    printf("Reorder_in    (%lf %lf %lf)\n",gtimers_avg.reorder_in/nproc,gtimers_min.reorder_in,gtimers_max.reorder_in);
    printf("Trans_exec    (%lf %lf %lf)\n",gtimers_avg.trans_exec/nproc,gtimers_min.trans_exec,gtimers_max.trans_exec);
    printf("Packsend      (%lf %lf %lf)\n",gtimers_avg.packsend/nproc,gtimers_min.packsend,gtimers_max.packsend);
    printf("Packsend_trans(%lf %lf %lf)\n",gtimers_avg.packsend_trans/nproc,gtimers_min.packsend_trans,gtimers_max.packsend_trans);
    printf("Unpackrecv    (%lf %lf %lf)\n",gtimers_avg.unpackrecv/nproc,gtimers_min.unpackrecv,gtimers_max.unpackrecv);
    printf("Alltoall      (%lf %lf %lf)\n",gtimers_avg.alltoall/nproc,gtimers_min.alltoall,gtimers_max.alltoall);
#ifdef CUDA
    printf("GPU Transfer      (%lf %lf %lf)\n",gtimers_avg.gpu_transfer/nproc,gtimers_min.gpu_transfer,gtimers_max.gpu_transfer);
#endif


  }

}

  extern  timer timers;

#endif


template <class Type> MPIplan<Type>::~MPIplan()
{
  delete [] SndCnts,SndStrt,RcvCnts,RcvStrt;
  delete grid1,grid2;
}

ProcGrid::ProcGrid(int procdims[3],MPI_Comm mpi_comm_init)
{
  int i,j,k;

  nd = 0;
  MPI_Comm_dup(mpi_comm_init,&mpi_comm_glob);
  for(i=0;i <3; i++) {
    if((ProcDims[i] = procdims[i]) > 1)
      nd++;
  }
  if(nd == 0)
    nd = 1;

  MPI_Comm_rank(mpi_comm_glob,&taskid);
  MPI_Comm_size(mpi_comm_glob,&numtasks);

  // Set up communicators for pencils or 3D decomposition
  // if(nd >1) {
  int periodic[3];
  int reorder=0;
  
  for(i=0;i < 3;i++)
    periodic[i]=1;
  MPI_Cart_create(mpi_comm_glob,3,ProcDims,periodic,reorder,&mpi_comm_cart);
  MPI_Cart_coords(mpi_comm_cart,taskid,3,grid_id_cart);
  
  MPI_Comm mpi_comm_tmp;
  int remain_dims[] = {1,0,0};
  for(i=0;i < 3;i++) {
    
    // Create COLUMN sub comm.
    MPI_Cart_sub(mpi_comm_cart,remain_dims,mpicomm+i);
    MPI_Comm_rank(mpicomm[i],grid_id_cart+i);
    
#ifdef DEBUG
    printf("%d: myid=%d %d\n",taskid,grid_id_cart[i]);
#endif
    remain_dims[(i+1)%3] = 1;
    remain_dims[i] = 0;
  }

}

ProcGrid::ProcGrid(const ProcGrid &rhs) 
  {
    
    //  if(rhs.is_set) {
    // is_set = true;

    //    prec = rhs.prec;
    int i,j,m,l;
    nd = rhs.nd;
    //    mpi_comm_glob = rhs.mpi_comm_glob;
    MPI_Comm_dup(rhs.mpi_comm_glob,&mpi_comm_glob);
    MPI_Comm_dup(rhs.mpi_comm_cart,&mpi_comm_cart);
    for(i=0;i<3;i++)
      MPI_Comm_dup(rhs.mpicomm[i],&mpicomm[i]);

    for(i=0;i<3; i++) {
      ProcDims[i] = rhs.ProcDims[i];
      grid_id_cart[i] = rhs.grid_id_cart[i];
    }
    taskid = rhs.taskid;
    numtasks = rhs.numtasks;
    //}
  }

 // grid constructor: initialize grid description and setup up MPI structures
  DataGrid::DataGrid(int gdims[3],int dim_conj_sym_,ProcGrid *Pgrid_,int dmap[3],int mem_order[3])
{
  int i,j;

  dim_conj_sym = dim_conj_sym_;
  Pgrid = Pgrid_;
  nd = Pgrid->nd;

  // Find dimension of processor grid (1 to 3 non-unit values in pgrid)
  for(i=0;i <3; i++) {
    Gdims[i] = gdims[i];
    Dmap[i] = dmap[i];
    Pdims[i] = Pgrid->ProcDims[Dmap[i]];
    MemOrder[i] = mem_order[i];
  }

  // Find grid id 3D coordinates of the local grid withint the global grid
  for(i=0;i<3;i++)
    grid_id[i] = Pgrid->grid_id_cart[Dmap[i]];

  // Allocate structures for local data sizes

  for(i=0;i<3;i++) {
    st[i] = new int[Pdims[i]];
    sz[i] = new int[Pdims[i]];
    en[i] = new int[Pdims[i]];
  }

  InitPencil(); // Initialize pencils/slabs etc
  is_set = true;
}

  // Copy constructor
DataGrid::DataGrid(const DataGrid &rhs)
{
  
  if(rhs.is_set) {
    is_set = true;

    //    prec = rhs.prec;
    int i,j,m,l;
    dim_conj_sym = rhs.dim_conj_sym;
    nd = rhs.nd;
    Pgrid = rhs.Pgrid;

    for(i=0;i<3; i++) {
      Gdims[i] = rhs.Gdims[i];
      Ldims[i] = rhs.Ldims[i];
      Pdims[i] = rhs.Pdims[i];
      MemOrder[i] = rhs.MemOrder[i];
      L[i] = rhs.L[i];
      D[i] = rhs.D[i];
      Dmap[i] = rhs.Dmap[i];
      grid_id[i] = rhs.grid_id[i];
      GlobStart[i] = rhs.GlobStart[i];
    }

    for(i=0;i<3;i++) {
      st[i] = new int[Pdims[i]];
      sz[i] = new int[Pdims[i]];
      en[i] = new int[Pdims[i]];
      
      for(j=0; j < Pdims[i]; j++) {
	//	st[i][j] = new int[3];
	//sz[i][j] = new int[3];
	//en[i][j] = new int[3];
	//for(int k=0;k<3; k++) {
	st[i][j] = rhs.st[i][j];
	sz[i][j] = rhs.sz[i][j];
	en[i][j] = rhs.en[i][j];
	
      }
    }
    
  }
  else
    is_set = false;
}

DataGrid::~DataGrid() 
{
  if(is_set) {
    for(int i=0;i < 3;i++) {
      //      for(int j=0;j<ProcGrid[i];j++) {
	delete [] st[i];
	delete [] sz[i];
	delete [] en[i];
    }
    
    /*
    delete [] st;
    delete [] sz;
    delete [] en;
    */
  }

}

ProcGrid::~ProcGrid()
{
  for(int i=0;i<3;i++)
    MPI_Comm_free(&mpicomm[i]);
  MPI_Comm_free(&mpi_comm_cart);
  MPI_Comm_free(&mpi_comm_glob);
}


void DataGrid::InitPencil()
{
  int i,j,k,n,pm,l,size,nl,nu,data,proc,li,myproc,pdim,ploc,p;

  // First determine local dimension(s)
  ploc=0;n=0;
  for(k=0; k < 3; k++) 
    if(Pdims[k] == 1) 
      L[ploc++] = k;
    else
      D[n++] = k;
  
  /*
    for(i=0;i<nd;i++)
      if(D[i] == k)
	break;
    if(i >= nd) {
      L[ploc++] = k;
      i = -1;
    }
    }*/
      //      i = proc_order[pdim++];

  // Distribute the dimension each communicator spans
  for(i=0;i<3;i++) {
    //    j = D[i];
    //    for(j=0;j<3;j++) {
    
    data = Gdims[i]; //[li];
    proc = Pdims[i];
    size = data/proc;
    nu = data%proc;
    nl = proc-nu;
    
    st[i][0] = 0;
    sz[i][0] = size;
    en[i][0] = size;
    for(p=1; p < nl; p++) {
      st[i][p] = st[i][p-1] +size;
      sz[i][p] = size;
      en[i][p] = en[i][p-1] +size;
    }
    size++;
    for(;p < proc; p++) {
      st[i][p] = st[i][p-1] +sz[i][p-1];
      sz[i][p] = size;
      en[i][p] = en[i][p-1] +size;
    }
    //st[i][p][j] = st[i][p-1][j] + size;
    //sz[i][p][j] = size;
    en[i][p-1] = data;
    sz[i][p-1] = data - st[i][p-1];
    
    // Assign local size for this spanning dimension
    myproc = grid_id[i];
    Ldims[i] = sz[i][myproc];
    GlobStart[i] = st[i][myproc];
  }
}


// Initialize a 3D transform type based on IDs of 1D transforms
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

  //  dt1 = types1D[types[0]]->dt1;
  //dt2 = types1D[types[2]]->dt2;
  prec = types1D[types[0]]->prec;
  is_set = true;
}

// Initialize a 3D transform type based on 1D transforms
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
  //  dt1 = types1D[types[0]]->dt1;
  //dt2 = types1D[types[2]]->dt2;
  prec = types1D[types[0]]->prec;
  is_set = true;

}

// Copy constructor
trans_type3D::trans_type3D(const trans_type3D &tr) {
  //  name = new char[sizeof(tr.name)+1];
  //strcpy(name,tr.name);
  for(int i=0;i<3;i++) 
    types[i] = tr.types[i];

  //  dt1 = tr.dt1;
  //dt2 = tr.dt2;
  is_set = true;
  prec = tr.prec;
} 


trans_type3D::~trans_type3D()
{
  //delete [] name;
}


}
