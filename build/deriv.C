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

#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

#include "p3dfft.h"

namespace p3dfft {

template <class Type> void  compute_deriv(Type *IN,Type *OUT,grid *gr, int idir)
{
  int ldir,i,j,k,g,mid,sdims[3];

  // Find local storage dimension to be differentiated and storage array dimensions
  for(i=0;i<3;i++) {
    sdims[gr->mem_order[i]] = gr->ldims[i];
    if(gr->mem_order[i] == idir)
      ldir = i;
  }

  // Adjust for reduced X space after real-to-complex transform
  if(gr->dim_conj_sym == idir)
    g = (gr->gdims[idir]-1)*2;
  else
    g = gr->gdims[idir];

  //  Compute the middle point in the spectrum (Nyquist frequency)
  mid = g/2 -gr->glob_start[idir];

  //  if(typeid(Type) == type_Type) {
    Type *p1,*p2,mult;

    switch(ldir) {
    case 0:
      p1 = IN;
      p2 = OUT;
      for(k=0;k<sdims[2];k++)
	for(j=0;j<sdims[1];j++) {
// Lower half: complex-multiply by i*k
	  for(i=0;i < min(sdims[0],mid);i++) 
	    *p2++ = Type(0.0,i+gr->glob_start[idir]) * *p1++;

// Nyquist frequency: zero
	  if(mid >= 0 && mid < sdims[ldir]) {
	    *p2++=0; p1++;
	  }

// Upper half: complex-multiply by i*(k-N)
	  for(i=max(0,mid+1);i < sdims[0];i++)
	    *p2++ = Type(0.0,i+gr->glob_start[idir] - g) * *p1++;
	}
	  
      break;
 
    case 1:
      p1 = IN;
      p2 = OUT;
      for(k=0;k<sdims[2];k++) {
// Lower half: complex-multiply by i*k
	for(j=0;j < min(sdims[1],mid);j++) { 
	  mult = Type(0.0,j+gr->glob_start[idir]);
	  for(i=0;i<sdims[0];i++) 
	    *p2++ = mult * *p1++;
	}
    // Nyquist frequency: zero
	if(mid >= 0 && mid < sdims[ldir]) 
	  for(i=0;i<sdims[0];i++) {
	    *p2++=0; p1++;
	  }

      // Upper half: complex-multiply by i*(k-N)
	for(j=max(0,mid+1);j < sdims[1];j++){
	  mult = Type(0.0,j+gr->glob_start[idir] -g);
	  for(i=0;i<sdims[0];i++) 
	    *p2++ = mult * *p1++;
	}
      }
    
    break;

    case 2:
      p1 = IN;
      p2 = OUT;
      // Lower half: complex-multiply by i*k
      for(k=0;k < min(sdims[2],mid);k++) { 
	mult = Type(0.0,k+gr->glob_start[idir]);
	for(j=0;j<sdims[1];j++)
	  for(i=0;i<sdims[0];i++) 
	    *p2++ = mult * *p1++;
      }
// Nyquist frequency: zero
      if(mid >= 0 && mid < sdims[ldir]) 
	for(j=0;j<sdims[1];j++)
	  for(i=0;i<sdims[0];i++) {
	    *p2++=0; p1++;
	  }

// Upper half: complex-multiply by i*(k-N)
      for(k=max(0,mid+1);k < sdims[2];k++){
	mult = Type(0.0,k+gr->glob_start[idir] -g);
	for(j=0;j<sdims[1];j++)
	  for(i=0;i<sdims[0];i++) 
	    *p2++ = mult * *p1++;
      } 
  	  
    break;
    
    }
  
}

  template void compute_deriv<mycomplex>(mycomplex *,mycomplex *,grid *,int);
  template void compute_deriv<complex_double>(complex_double *,complex_double *,grid *,int);

}

