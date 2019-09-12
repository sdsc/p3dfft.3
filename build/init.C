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
vector<grid> stored_grids; // Defined grids (used in Fortran and C wrappers)


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

void setup()
{
  const char *name;
  int isign;
  gen_trans_type *p;

  int types_count=0;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Single" << endl;
#endif
  name = "Empty Type Single";
  p = new gen_trans_type(name,1,1,4,0);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_SINGLE = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Double" << endl;
#endif
  name = "Empty Type Double";
  p = new gen_trans_type(name,1,1,8,0);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_DOUBLE = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty Type Single Complex" << endl;
#endif
  name = "Empty Type Single Complex";
  p = new gen_trans_type(name,2,2,4,0);
  p->is_empty = true;
  types1D.push_back(p);


  EMPTY_TYPE_SINGLE_COMPLEX = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dfft_setup: adding Empty TypeDouble Complex" << endl;
#endif
  name = "Empty Type Double Complex";
  p = new gen_trans_type(name,2,2,8,0);
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
  cout << "p3dft_setup: adding Cosine R2R DCT1 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT1, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dct1_s, (void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT1_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT1 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT1, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dct1_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT1_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT1 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT1, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(long (*)(...) ) plan_dct1_complex_s, (void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT1_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT1 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT1, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dct1_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT1_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST1 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST1, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dst1_s,(void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST1_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST1 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST1, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dst1_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST1_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST1 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST1, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_dst1_complex_s,(void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST1_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST1 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST1, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dst1_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST1_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT2 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT2, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dct2_s, (void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT2_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT2 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT2, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dct2_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT2_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT2 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT2, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(long (*)(...) ) plan_dct2_complex_s, (void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT2_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT2 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT2, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dct2_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT2_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST2 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST2, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dst2_s,(void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST2_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST2 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST2, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dst2_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST2_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST2 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST2, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_dst2_complex_s,(void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST2_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST2 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST2, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dst2_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST2_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT3 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT3, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dct3_s, (void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT3_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT3 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT3, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dct3_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT3_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT3 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT3, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(long (*)(...) ) plan_dct3_complex_s, (void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT3_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT3 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT3, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dct3_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT3_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST3 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST3, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dst3_s,(void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST3_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST3 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST3, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dst3_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST3_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST3 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST3, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_dst3_complex_s,(void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST3_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST3 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST3, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dst3_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST3_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding Cosine R2R DCT4 single type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT4, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dct1_s, (void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DCT4_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding cosine R2R DCT4 double type" << endl;
#endif

  name = "Real-valued Cosine Transform DCT4, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dct1_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DCT4_REAL_D = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex Cosine R2R DCT4 single type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT4, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex >(name,(long (*)(...) ) plan_dct1_complex_s, (void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DCT4_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex cosine R2R DCT4 double type" << endl;
#endif

  name = "Complex-valued Cosine Transform DCT4, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dct1_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DCT4_COMPLEX_D = types_count;
  types_count++;



#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST4 single type" << endl;
#endif

  name = "Real-valued Sine Transform DST4, single precision";
#ifdef FFTW
  p = new trans_type1D<float,float>(name,(long (*)(...) ) plan_dst4_s,(void (*)(...)) exec_r2r_s);
#endif
  types1D.push_back(p);
  DST4_REAL_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding sine R2R DST4 double type" << endl;
#endif

  name = "Real-valued Sine Transform DST4, double precision";
#ifdef FFTW
  p = new trans_type1D<double,double>(name,(long (*)(...) ) plan_dst4_d,(void (*)(...)) exec_r2r_d);
#endif
  types1D.push_back(p);
  DST4_REAL_D = types_count;
  types_count++;


#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST4 single type" << endl;
#endif

  name = "Complex-valued Sine Transform DST4, single precision";
#ifdef FFTW
  p = new trans_type1D<mycomplex,mycomplex>(name,(long (*)(...) ) plan_dst4_complex_s,(void (*)(...)) exec_r2r_complex_s);
#endif
  types1D.push_back(p);
  DST4_COMPLEX_S = types_count;
  types_count++;

#ifdef DEBUG
  cout << "p3dft_setup: adding complex sine R2R DST4 double type" << endl;
#endif

  name = "Complex-valued Sine Transform DST4, double precision";
#ifdef FFTW
  p = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) plan_dst4_complex_d,(void (*)(...)) exec_r2r_complex_d);
#endif
  types1D.push_back(p);
  DST4_COMPLEX_D = types_count;
  types_count++;


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



}

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

  if(taskid == 0) {
    printf("TIMERS (avg/min/max): \nReorder_trans (%lf %lf %lf)\n",gtimers_avg.reorder_trans/nproc,gtimers_min.reorder_trans,gtimers_max.reorder_trans);
    printf("Reorder_out   (%lf %lf %lf)\n",gtimers_avg.reorder_out/nproc,gtimers_min.reorder_out,gtimers_max.reorder_out);
    printf("Reorder_in    (%lf %lf %lf)\n",gtimers_avg.reorder_in/nproc,gtimers_min.reorder_in,gtimers_max.reorder_in);
    printf("Trans_exec    (%lf %lf %lf)\n",gtimers_avg.trans_exec/nproc,gtimers_min.trans_exec,gtimers_max.trans_exec);
    printf("Packsend      (%lf %lf %lf)\n",gtimers_avg.packsend/nproc,gtimers_min.packsend,gtimers_max.packsend);
    printf("Packsend_trans(%lf %lf %lf)\n",gtimers_avg.packsend_trans/nproc,gtimers_min.packsend_trans,gtimers_max.packsend_trans);
    printf("Unpackrecv    (%lf %lf %lf)\n",gtimers_avg.unpackrecv/nproc,gtimers_min.unpackrecv,gtimers_max.unpackrecv);
    printf("Alltoall      (%lf %lf %lf)\n",gtimers_avg.alltoall/nproc,gtimers_min.alltoall,gtimers_max.alltoall);


  }

}

  extern  timer timers;

#endif



void cleanup()
{
  
  //  vector<grid>::iterator it1=stored_grids.begin();
  //  stored_grids.erase(stored_grids.begin(),stored_grids.end());
  //  for(vector<trans_type3D>::iterator it=types3D.begin();it != types3D.end();it++) {
  //    types3D.erase(types3D.begin(),types3D.end());
  stored_grids.clear();
  types3D.clear();

  // Since these are vectors of pointers, simply erasing is not enough; must delete each referenced class. 
  //  printf("Clearing types1D\n");
  vector<gen_trans_type *>::iterator it1=types1D.begin();
  while(it1 != types1D.end()) {
    delete *it1;
    it1 = types1D.erase(it1);
  }

  //  printf("Clearing Plans\n");
  //  for(vector<Plan *>::iterator it=Plans.begin();it != Plans.end();it++) {
  vector<Plan *>::iterator it=Plans.begin();
  while(it != Plans.end()) {
#ifdef FFTW
    if((*it)->libplan_in != NULL) 
      fftw_destroy_plan((fftw_plan) (*it)->libplan_in);
    if((*it)->libplan_out != NULL) 
      fftw_destroy_plan((fftw_plan) (*it)->libplan_out);
#endif
    delete *it;
    it = Plans.erase(it);
  }

  //  printf("Clearing stored_trans1D\n");
  vector<stage *>::iterator it2=stored_trans1D.begin();
  while(it2 != stored_trans1D.end()) {
    delete *it2;
    it2 = stored_trans1D.erase(it2);
  }

  //  printf("Clearing stored_trans3D\n");
  vector<gen_transform3D *>::iterator it3=stored_trans3D.begin();
  while(it3 != stored_trans3D.end()) {
    delete *it3;
    it3 = stored_trans3D.erase(it3);
  }

#ifdef FFTW
  //  printf("Cleaning FFTW\n");
  fftw_cleanup();
#endif    

  //  printf("Done Cleaning\n");

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
void exec_r2r_s(long plan,float *in,float *out)
{
  fftwf_execute_r2r((fftwf_plan) plan,in,out);
}
void exec_r2r_d(long plan,double *in,double *out)
{
  fftw_execute_r2r((fftw_plan) plan,in,out);
}

void exec_r2r_complex_s(long plan,float *in,float *out)
{
  fftwf_execute_r2r((fftwf_plan) plan,in,out);
  fftwf_execute_r2r((fftwf_plan) plan,in+1,out+1);
}
void exec_r2r_complex_d(long plan,double *in,double *out)
{
  fftw_execute_r2r((fftw_plan) plan,in,out);
  fftw_execute_r2r((fftw_plan) plan,in+1,out+1);
}


long plan_dct1_s(int rank, const int *n,		   
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


long plan_dct1_d(int rank, const int *n,		   
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


long plan_dst1_s(int rank, const int *n,		   
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


long plan_dst1_d(int rank, const int *n,		   
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


long plan_dct1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT00;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dct1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT00;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst1_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT00;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst1_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT00;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



long plan_dct2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT10;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


long plan_dct2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT10;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst2_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT10;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst2_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT10;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dct2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT10;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dct2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT10;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst2_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT10;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst2_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT10;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



long plan_dct3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT01;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


long plan_dct3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT01;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst3_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT01;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst3_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT01;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dct3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT01;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dct3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT01;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst3_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT01;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst3_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT01;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}



long plan_dct4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT11;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed, \
    ostride,odist,&type,fft_flag));
}


long plan_dct4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT11;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst4_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT11;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dst4_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT11;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride,idist,out,onembed,\
    ostride,odist,&type,fft_flag));
}


long plan_dct4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		int ostride, int odist,unsigned fft_flag)
{
   fftwf_r2r_kind type = FFTW_REDFT11;
   return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed, \
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dct4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_REDFT11;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst4_complex_s(int rank, const int *n,		   
                         int howmany,					   
                         float *in, const int *inembed,			   
                         int istride, int idist,			   
                         float *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftwf_r2r_kind type = FFTW_RODFT11;
  return((long) fftwf_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
}


long plan_dst4_complex_d(int rank, const int *n,		   
                         int howmany,					   
                         double *in, const int *inembed,			   
                         int istride, int idist,			   
                         double *out, const int *onembed,			
		  int ostride, int odist,unsigned fft_flag)
{
  fftw_r2r_kind type = FFTW_RODFT11;
  return((long) fftw_plan_many_r2r(1,n,howmany,in,inembed,istride*2,idist*2,out,onembed,\
    ostride*2,odist*2,&type,fft_flag));
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


 // grid constructor: initialize grid description and setup up MPI structures
grid::grid(int gdims_[3],int dim_conj_sym_,int pgrid_[3],int proc_order_[3],int mem_order_[3],
	   MPI_Comm mpicomm_)
	   //,int prec_,int dt_)
{
  int i,j;
  int myid[2];
  MPI_Comm mpi_comm_tmp;

  dim_conj_sym = dim_conj_sym_;
  MPI_Comm_dup(mpicomm_,&mpi_comm_glob);
  nd=0;
  P[0]=P[1]=P[2]=1;
  D[0]=D[1]=D[2]=L[0]=L[1]=L[2]=-1;
  // Find dimension of processor grid (1 to 3 non-unit values in pgrid)

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

  // Set up communicators for pencils or 3D decomposition
  if(nd >1) {
    int periodic[nd];
    int reorder=0;

    for(i=0;i < nd;i++)
      periodic[i]=1;
    MPI_Cart_create(mpicomm_,nd,P,periodic,reorder,&mpi_comm_cart);
    MPI_Cart_coords(mpi_comm_cart,taskid,nd,grid_id_cart);

    int remain_dims[3];
    if(nd == 2) { // 2D (pencils) decomposition
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
    }
    else if(nd == 3) { // 3D (volumetric) decomsposition 
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
  }
  else { //nd=1 Slab (1D) decomposition
    grid_id_cart[0] = taskid;
    *mpicomm = mpicomm_;
  }

  // Find grid id 3D coordinates of the local grid withint the global grid
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

  // Allocate structures for local data sizes

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

  InitPencil(); // Initialize pencils/slabs etc
  is_set = true;
}

  // Copy constructor
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
    dim_conj_sym = rhs.dim_conj_sym;

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
    numtasks = rhs.numtasks;
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
  dt1 = types1D[types[0]]->dt1;
  dt2 = types1D[types[2]]->dt2;
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
  dt1 = types1D[types[0]]->dt1;
  dt2 = types1D[types[2]]->dt2;
  prec = types1D[types[0]]->prec;
  is_set = true;

}

// Copy constructor
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

// Define a 1D transform type
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

