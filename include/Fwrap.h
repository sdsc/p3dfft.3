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

  void p3dfft_init_3Dtype_f(int *,int[3]); //,char *);
  void p3dfft_plan_1Dtrans_f(int *,int *,int *,int *,int *);
  void p3dfft_plan_3Dtrans_f(int *,int *,int *,Type3D *);
void p3dfft_init_grid_f(int *,int *,int *,int *,int *,int *,int *,int *);
void p3dfft_exec_1Dtrans_double_f(int *,double *,double *, int *);
void p3dfft_exec_1Dtrans_single_f(int *,float *,float *, int *);
void p3dfft_exec_3Dtrans_double_f(Plan3D *,double *,double *, int *);
void p3dfft_exec_3Dtrans_single_f(Plan3D *,float *,float *, int *);
void p3dfft_exec_3Dderiv_double_f(Plan3D *,double *,double *,int *, int *);
void p3dfft_exec_3Dderiv_single_f(Plan3D *,float *,float *,int *, int *);
void p3dfft_compute_deriv_single_f(float *,float *,int *,int *);
void p3dfft_compute_deriv_double_f(double *,double *,int *,int *);
