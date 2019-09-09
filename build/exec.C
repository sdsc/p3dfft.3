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

namespace p3dfft {

#ifdef TIMERS
  timer timers;
#endif

#ifdef DEBUG
static int cnt_pack=0;
static int cnt_trans=0;
void printbuf(char *,int[3],int,int);
#endif

/* Execute 3D transform, optionally followed by spectral derivative in logical dimension idir.
   Input/output: buffers containing 3D arrays, contiguously stored according to memory mappings (mo1 and mo2) expected by transform3D class. Can be same buffer, otherwise need to be non-overlapping.
   OW: if != 0 input buffer can be overwritten (obviously the case if in==out).
*/

  template <class Type1,class Type2> void transform3D<Type1,Type2>::exec_deriv(Type1 *in,Type2 *out,int idir, bool OW) 
{
  int ldir,i,j,k,g,mid;
  complex_double *p1,*p2,mult;

  int nvar=0;
  char *buf[2],*var[2],*buf_in,*buf_out;
  int next,curr,dt_1,dt_2,nextvar,stage_cnt=0;

  if((void *) in == (void *) out && !OW) {
    printf("Warning in transform3D:exec_deriv: input and output are the same, but overwrite priviledge is not set\n");
    //    OW = 1;
  }

#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(grid1->mpi_comm_glob,&taskid);
#endif

  next = 1; curr = 0;
  buf[curr] = buf_in = (char *) in;
  buf_out = (char *) out;
  dt_1 = dt1;
  nextvar = 0;
  int cnt=0;

  for(stage *curr_stage=Stages;curr_stage != NULL;curr_stage = curr_stage->next) {
#ifdef DEBUG
    printf("stage %d, kind=%d\n",stage_cnt++,curr_stage->kind);
#endif  

	//	transplan<Type1,Type2> *st = (transplan<Type1,Type2> *) curr_stage;
    stage *st = curr_stage;
    int size1 = st->dims1[0] * st->dims1[1] * st->dims1[2]; 
    int size2 = st->dims2[0] * st->dims2[1] * st->dims2[2]; 
    dt_1 = curr_stage->dt1;
    dt_2 = curr_stage->dt2;
    if(!curr_stage->next) { 
      buf[next] = buf_out; // Last stage always writes to out destination
    }
	/*
	  if(buf[curr] == buf[next] && !st->inplace) {
	    printf("Error in transform3D::exec_deriv: stage %d was intended for out-place transform but is used in-place\n",stage_cnt);
	    MPI_Abort(grid1->mpi_comm_glob,0);
	  }
	  }*/
    else 
      if((!OW && buf[curr] == buf_in) || size2*dt_2 > size1 *dt_1) {
	//  out-place
	// check if we can use the destination array as the output buffer
	if(dt_2 * size2 <= dt2 * grid2->ldims[0] *grid2->ldims[1] *grid2->ldims[2]) {
	  buf[next] = buf_out;
#ifdef DEBUG
	  printf("%d: Using output buffer in out-of-place transform\n",taskid);
#endif
	}
	else { //  allocate a new buffer
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
	  var[nextvar] = new char[size2*dt_2*st->stage_prec];
	  buf[next] = var[nextvar];
	  nvar++;
	  nextvar = 1-nextvar;
	}
      }
      else // In-place
	buf[next] = buf[curr];
    
    if(curr_stage->kind == TRANS_ONLY) {
      // Only transform, no exchange
      
      transplan<Type1,Type2> *tr = (transplan<Type1,Type2> *) st;
     
      if(idir == tr->trans_dim)
      	st->myexec_deriv(buf[curr],buf[next], OW || buf[curr] != buf_in);
      else
        st->myexec(buf[curr],buf[next],OW || buf[curr] != buf_in);
      
#ifdef DEBUG
      int imo[3];
      char str[80];
      for(i=0;i<3;i++)
	imo[tr->mo2[i]] = i; 
      sprintf(str,"exec-out.%d.%d",stage_cnt,taskid);
      write_buf<Type2>((Type2 *) buf[next],str,curr_stage->dims2,imo);
#endif
    }
    else if(curr_stage->kind == MPI_ONLY) { // Only MPI plan (exchange, no transform)
      MPIplan<Type1> *tr = (MPIplan<Type1> *) curr_stage; 
      curr_stage->myexec(buf[curr],buf[next],true);
    }
    else { // MPI and transform combined

      trans_MPIplan<Type1,Type2> *tr = (trans_MPIplan<Type1,Type2> *) (transplan<Type1,Type2> *) curr_stage; 
      if(idir == tr->trplan->trans_dim)
      	curr_stage->myexec_deriv(buf[curr],buf[next],OW  || buf[curr] != buf_in);
      else
      	curr_stage->myexec(buf[curr],buf[next], OW  || buf[curr] != buf_in);

#ifdef DEBUG
      int imo[3];
      char str[80];
      for(i=0;i<3;i++)
	imo[tr->trplan->mo2[i]] = i; 
      sprintf(str,"exec-out.%d.%d",stage_cnt,taskid);
      write_buf<Type2>((Type2 *) buf[next],str,curr_stage->dims2,imo);
#endif
    }
	
    dt_1 = dt_2;
	
    if(nvar > 1) { // Keep track and delete buffers that are no longer used
      delete [] var[nextvar];
      nvar--;
    }

    next = 1-next;
    curr = 1-curr;
  }
  if(nvar > 0)
    delete [] var[1-nextvar];
}


  // Execute local spectral derivative
// dims are storage dimensions; differentiate in the first dimension, assumed to be local
  template <class Type1,class Type2> void transplan<Type1,Type2>::compute_deriv_loc(Type2 *in,Type2 *out,int dims[3]) 
  {
    int g,mid,i,j,k;
    bool r2c;

    if(trans_type->dt1 == 1 && trans_type->dt2 == 2)
      r2c = true;
    else
      r2c = false;


  // Adjust for reduced X space after real-to-complex transform
  if(r2c)
    g = (dims[0]-1)*2;
  else
    g = dims[0];

  //  Compute the middle point in the spectrum (Nyquist frequency)
  mid = g/2;

  if(typeid(Type2) == type_complex_double) {

    complex_double *p1 = (complex_double *) in;
  complex_double *p2 = (complex_double *) out;

  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++) {
// Lower half: complex-multiply by i*k
      for(i=0;i < mid;i++) 
	*p2++ = complex_double(0.0,i) * *p1++;

// Nyquist frequency: zero
      *p2++=0; p1++;

// Upper half: complex-multiply by i*(k-N)
      for(i=mid+1;i < dims[0];i++)
	*p2++ = complex_double(0.0,i - g) * *p1++;
    }  
  }

  else if(typeid(Type2) == type_complex) {

  mycomplex *p1 = (mycomplex *) in;
  mycomplex *p2 = (mycomplex *) out;

  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++) {
// Lower half: complex-multiply by i*k
      for(i=0;i < mid;i++) 
	*p2++ = mycomplex(0.0,i) * *p1++;

// Nyquist frequency: zero
      *p2++=0; p1++;

// Upper half: complex-multiply by i*(k-N)
      for(i=mid+1;i < dims[0];i++)
	*p2++ = mycomplex(0.0,i - g) * *p1++;
    }  
  }
}

  // Execute 3D transform
  template <class Type1,class Type2> void transform3D<Type1,Type2>::exec(Type1 *in,Type2 *out, bool OW)
//void transform3D::exec(char *in,char *out,int OW)
{
  // Call exec_deriv with -1 (void) for derivative dimension
  exec_deriv(in,out,-1, OW);
}

  void stage::myexec(char *in,char *out, bool OW)
{
  switch(kind) {
  case TRANS_ONLY: 
    if(dt1 == REAL)
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<float,float> *st=(transplan<float,float> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  transplan<double,double> *st=(transplan<double,double> *) this;
	  st->exec(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  transplan<float,mycomplex> *st=(transplan<float,mycomplex> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  transplan<double,complex_double> *st=(transplan<double,complex_double> *) this;
	  st->exec(in,out,OW);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<mycomplex,float> *st=(transplan<mycomplex,float> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  transplan<complex_double,double> *st=(transplan<complex_double,double> *) this;
	  st->exec(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  transplan<mycomplex,mycomplex> *st=(transplan<mycomplex,mycomplex> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  transplan<complex_double,complex_double> *st=(transplan<complex_double,complex_double> *) this;
	  st->exec(in,out,OW);
	}
    
    break;
    
  case MPI_ONLY:
    if(dt1 == REAL)
      if(stage_prec == 4) {
	MPIplan<float> *st=(MPIplan<float> *) this;
	st->exec(in,out);
      }
      else {
	MPIplan<double> *st=(MPIplan<double> *) this;
	st->exec(in,out);
      }
    else
      if(stage_prec == 4) {
	MPIplan<mycomplex> *st=(MPIplan<mycomplex> *) this;
	st->exec(in,out);
      }
      else {
	MPIplan<complex_double> *st=(MPIplan<complex_double> *) this;
	st->exec(in,out);
      }
    
    break;
  case TRANSMPI:
    if(dt1 == REAL)
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  trans_MPIplan<float,float> *st=(trans_MPIplan<float,float> *) (transplan<float,float> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  trans_MPIplan<double,double> *st=(trans_MPIplan<double,double> *) (transplan<double,double> *) this;
	  st->exec(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<float,mycomplex> *st=(trans_MPIplan<float,mycomplex> *) (transplan<float,mycomplex> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  trans_MPIplan<double,complex_double> *st=(trans_MPIplan<double,complex_double> *) (transplan<double,complex_double> *) this;
	  st->exec(in,out,OW);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,float> *st=(trans_MPIplan<mycomplex,float> *) (transplan<mycomplex,float> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  trans_MPIplan<complex_double,double> *st=(trans_MPIplan<complex_double,double> *) (transplan<complex_double,double> *) this;
	  st->exec(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,mycomplex> *st=(trans_MPIplan<mycomplex,mycomplex> *)  (transplan<mycomplex,mycomplex> *) this;
	  st->exec(in,out,OW);
	}
	else {
	  trans_MPIplan<complex_double,complex_double> *st=(trans_MPIplan<complex_double,complex_double> *) (transplan<complex_double,complex_double> *) this;
	  st->exec(in,out,OW);
	}
    break;
  }

}

  void stage::myexec_deriv(char *in,char *out, bool OW)
{
  switch(kind) {
  case TRANS_ONLY: 
    if(dt1 == REAL)
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<float,float> *st=(transplan<float,float> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  transplan<double,double> *st=(transplan<double,double> *) this;
	  st->exec_deriv(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  transplan<float,mycomplex> *st=(transplan<float,mycomplex> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  transplan<double,complex_double> *st=(transplan<double,complex_double> *) this;
	  st->exec_deriv(in,out,OW);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<mycomplex,float> *st=(transplan<mycomplex,float> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  transplan<complex_double,double> *st=(transplan<complex_double,double> *) this;
	  st->exec_deriv(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  transplan<mycomplex,mycomplex> *st=(transplan<mycomplex,mycomplex> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  transplan<complex_double,complex_double> *st=(transplan<complex_double,complex_double> *) this;
	  st->exec_deriv(in,out,OW);
	}
    
    break;
    
  case MPI_ONLY:
    if(dt1 == REAL)
      if(stage_prec == 4) {
	MPIplan<float> *st=(MPIplan<float> *) this;
	st->exec(in,out);
      }
      else {
	MPIplan<double> *st=(MPIplan<double> *) this;
	st->exec(in,out);
      }
    else
      if(stage_prec == 4) {
	MPIplan<mycomplex> *st=(MPIplan<mycomplex> *) this;
	st->exec(in,out);
      }
      else {
	MPIplan<complex_double> *st=(MPIplan<complex_double> *) this;
	st->exec(in,out);
      }
    
    break;
  case TRANSMPI:
    if(dt1 == REAL)
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  trans_MPIplan<float,float> *st=(trans_MPIplan<float,float> *) (transplan<float,float> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  trans_MPIplan<double,double> *st=(trans_MPIplan<double,double> *) (transplan<double,double> *) this;
	  st->exec_deriv(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<float,mycomplex> *st=(trans_MPIplan<float,mycomplex> *) (transplan<float,mycomplex> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  trans_MPIplan<double,complex_double> *st=(trans_MPIplan<double,complex_double> *) (transplan<double,complex_double> *) this;
	  st->exec_deriv(in,out,OW);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,float> *st=(trans_MPIplan<mycomplex,float> *) (transplan<mycomplex,float> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  trans_MPIplan<complex_double,double> *st=(trans_MPIplan<complex_double,double> *) (transplan<complex_double,double> *) this;
	  st->exec_deriv(in,out,OW);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,mycomplex> *st=(trans_MPIplan<mycomplex,mycomplex> *)  (transplan<mycomplex,mycomplex> *) this;
	  st->exec_deriv(in,out,OW);
	}
	else {
	  trans_MPIplan<complex_double,complex_double> *st=(trans_MPIplan<complex_double,complex_double> *) (transplan<complex_double,complex_double> *) this;
	  st->exec_deriv(in,out,OW);
	}
    break;
  }

}


  // Execute 1D transform, combined with local transpose as needed
// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]

  template <class Type1,class Type2> void transplan<Type1,Type2>::exec(char *in_,char *out_, bool OW)
{
  int L,N,m,mocurr[3],mc[3];
  Type1 *in;
  Type2 *out;
  Type2 *buf=NULL;
  void swap0(int newmo[3],int mo[3],int L);

  in = (Type1 *) in_;
  out = (Type2 *) out_;

  L = trans_dim;

  

#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf("%d: In transplan::exec, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  if(mo1[L] == 0 || mo2[L] == 0) {

    if(!arcmp(mo1,mo2,3)) { // If there is no change in memrory mapping, there is no need to do transpose. Just call 1D transform.
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      if((void *) in != (void *) out)
	if(trans_type->is_empty) 
	  memcpy(out,in,prec*dims1[0]*dims1[1]*dims1[2]);
	
	else 
	  /*
#ifdef FFTW
	  if(dt1 > dt2 && !OW) { // C2R
	    Type1 *buf = new Type1[dims1[0]*dims1[1]*dims1[2]];
	    memcpy(buf,in,prec*dims1[0]*dims1[1]*dims1[2]);
	    (*(trans_type->exec))(plan->libplan_out,buf,out);
	    delete [] buf;
	  }
	  else
#endif
	  */
	  (*(trans_type->exec))(plan->libplan_out,in,out);
      else if(!trans_type->is_empty) 

	//      if(!OW || dt2* dims2[0]*dims2[1]*dims2[2] > dt1 *dims1[0]*dims1[1]*dims1[2]) { // out-of-place
      
	  (*(trans_type->exec))(plan->libplan_in,in,out);

#ifdef TIMERS
      timers.trans_exec += MPI_Wtime() -t1;
#endif
    }
    else { // Otherwise, combine transpose with transform
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      if(trans_type->is_empty) {
	if((void *) in == (void *) out)
	  reorder_in<Type1>(in,mo1,mo2,dims1);   
	else
	  reorder_out((Type2 *)in,out,mo1,mo2,dims1);   
      }
      else
	reorder_trans(in,out,mo1,mo2,dims1,OW);   
#ifdef TIMERS
      timers.reorder_trans += MPI_Wtime() -t1;
#endif
    }	
  }

  else  //If first dimension is not the one to be transformed, need to reorder first, combined with transform

      if(trans_type->is_empty) {
	if((void *) in == (void *) out)
	  reorder_in<Type1>(in,mo1,mo2,dims1);   
	else
	  reorder_out((Type2 *) in,out,mo1,mo2,dims1);   
      }
      else {

    swap0(mocurr,mo1,L);
    buf = new Type2[dims2[0]*dims2[1]*dims2[2]];
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      reorder_trans(in,buf,mo1,mocurr,dims1,OW);   
#ifdef TIMERS
      timers.reorder_trans += MPI_Wtime() -t1;
#endif
#ifdef TIMERS
      t1=MPI_Wtime();
#endif
    reorder_out(buf,out,mocurr,mo2,dims2);
#ifdef TIMERS
      timers.reorder_out += MPI_Wtime() -t1;
#endif
    delete [] buf;
  }
}

// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]
// Transform in 1D and take spectral derivative in the dimension of transform
template <class Type1,class Type2> void transplan<Type1,Type2>::exec_deriv(char *in_,char *out_, bool OW)
{
  int L,N,m,mocurr[3],mc[3];
  Type1 *in;
  Type2 *out;
  Type2 *buf=NULL;
  bool alloc=false;
  void swap0(int newmo[3],int mo[3],int L);

  in = (Type1 *) in_;
  out = (Type2 *) out_;

  L = trans_dim;

  if(dt2 != 2)
    printf("Error in exec_deriv: expected complex output type\n");

#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf("%d: In transplan::exec, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  if(mo1[L] == 0 || mo2[L] == 0) {

    if(!arcmp(mo1,mo2,3)) {  // If there is no change in memory ordering
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif

      if((void *) in == (void *) out) 
  // Do 1D transform directly on input, into 'out'
	(*(trans_type->exec))(plan->libplan_in,in,out); 
      else
	(*(trans_type->exec))(plan->libplan_out,in,out); 
      int sdims[3]; // Find storage dimensions
      for(int i=0;i<3;i++)
	sdims[mo1[i]] = grid2->ldims[i];
      //      bool r2c = (grid2->dim_conj_sym >= 0);
      // Now compute local derivative
      compute_deriv_loc(out,out,sdims);
#ifdef TIMERS
      timers.trans_deriv += MPI_Wtime() -t1;
#endif
    }
    else { // Combine reordering with 1D transform and derivative
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      reorder_deriv(in,out,mo1,mo2,dims1,OW);   
#ifdef TIMERS
      timers.reorder_deriv += MPI_Wtime() -t1;
#endif
    }	
  }

  else { //If first dimension is not the one to be transformed, need to reorder first, combined with transform
    swap0(mocurr,mo1,L);
  //    Try to organize so out of place transpose is used, if possible
    //    if(mocurr[0] != mo2[0] || mocurr[1] != mo2[1] || mocurr[2] != mo2[2]) {
    buf = new Type2[dims2[0]*dims2[1]*dims2[2]];
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      reorder_deriv(in,buf,mo1,mocurr,dims1,OW);   
#ifdef TIMERS
      timers.reorder_deriv += MPI_Wtime() -t1;
      t1=MPI_Wtime();
#endif
    reorder_out(buf,out,mocurr,mo2,dims2);
#ifdef TIMERS
      timers.reorder_out += MPI_Wtime() -t1;
#endif
    delete [] buf;
  }
}

#ifdef DEBUG
  void printbuf(char *buf,int dims[3],int dt,int taskid)
{
  int i,j,k;
  complex_double x;

  for(k=0;k<dims[2];k++)
    for(j=0;j<dims[1];j++)
      for(i=0;i<dims[0];i++) {
	x = *((complex_double *) (buf+8*dt*(i+dims[0]*j+dims[0]*dims[1]*k)));
	if(abs(x) > 0.0001)
	  printf("%d: (%d %d %d) %lg %lg\n",taskid,i,j,k,x.real(),x.imag());
      }

}
#endif
  
int arcmp(int *A,int *B,int N)
{
  int i;
  for(i=0; i < N; i++)
    if(A[i] != B[i])
      return(1);

  return(0);
}

#define TRANS_IN 0
#define TRANS_OUT 1

  // Reorder (local transpose), followed by 1D transform
// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);

  int imo1[3],imo2[3],d1[3],d2[3];
  int scheme,nb13,nb31,nb23,nb32;
  
#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf("%d: In transplan::reorder_trans, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  // Find inverse memory mappings
  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  // Find local storage dimensions
  for(i=0;i<3;i++) {
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }
 
  /* 
    if(dims1[trans_dim] != *inembed) {
      printf("Error in reorder_trans: leading dimension on input %d doesn't match transform type %d\n",dims1[trans_dim],*inembed);
      return;
    }
    if(dims2[trans_dim] != *onembed) {
      printf("Error in reorder_trans: leading dimension on output %d doesnt match transform type %d\n", dims2[trans_dim],*onembed);
      return;
    }
  */

    if(mo1[trans_dim] == 0) 
    scheme = TRANS_IN;
  else if(mo2[trans_dim] == 0) 
    scheme = TRANS_OUT;
  else {
    printf("Error in reorder_trans: expected dimension %d to be the leading dimension for input or output\n",trans_dim);
    return;
   }

  rel_change(imo1,imo2,mc);

#ifdef DEBUG
  printf("In reorder_trans, mc=%d %d %d, scheme=%s\n",mc[0],mc[1],mc[2],(scheme == TRANS_IN) ? "IN" : "OUT"); 
  char str[80];
  static int cnt_reorder_trans=0;
#endif

  bool inplace=false;

  if(scheme == TRANS_IN) { // Scheme is transform before reordering

    Type1 *pin_t1,*pIN;
    Type2 *tmp,*pin2,*pout2;
    float *pin,*pin1,*pout,*pout1,*ptran2; 
    int ds=sizeof(Type2)/4;
#ifdef FFTW
    if(dt1 > dt2 && !OW) {
      pIN = new Type1[d1[0]*d1[1]*d1[2]];
      memcpy(pIN,in,sizeof(Type1)*d1[0]*d1[1]*d1[2]);
    }
    else
#endif
      pIN = in;

   switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
	
	if(OW && dt2 <= dt1)
	  inplace=true;
	else
	  tmp = new Type2[d2[0]*d2[1]];
      	
	for(k=0;k <d1[2];k++) {
	  if(inplace) {
	    tmp = (Type2 *) (in+k*d1[0]*d1[1]);
	    (*(trans_type->exec))(plan->libplan_in,pIN+k*d1[0]*d1[1],tmp);
	  }
	  else
	    (*(trans_type->exec))(plan->libplan_out,pIN+k*d1[0]*d1[1],tmp);

	  pout = (float *) out+ds*k*d2[0]*d2[1];
#ifdef MKL_BLAS
#ifdef DEBUG
	  if(taskid == 0) {
	    printf("A=\n");
	    for(i=0;i<d2[0];i++) {
	      for(j=0;j<d2[1];j++)
		printf("%lf ",*(tmp+i*d2[1]+j));
	      printf("\n");
	    }
	  }
#endif
	  blas_trans<Type2>(d2[1],d2[0],1.0,tmp,d2[1],pout,d2[0]);
#ifdef DEBUG
	  if(taskid == 0) {
	    printf("B=\n");
	    for(i=0;i<d2[1];i++) {
	      for(j=0;j<d2[0];j++)
		printf("%lf ",*(pout+i*d2[0]+j));
	      printf("\n");
	    }
	  }
#endif
	
#else

	  pout = (float *) out +ds*k*d2[0]*d2[1];
	  for(j=0;j < d2[1];j++) {
	    pin = (float *) tmp + ds*j;
	    for(i=0;i < d2[0];i++) {
	      for(int m=0;m<ds;m++)
		*pout++ = *pin++;
	      pin += ds*(d2[1]-1);
	    }	
	 
	  }
#endif
	  //      if(!trans_type->is_empty)
	}
      	
	if(!inplace)
	  delete [] tmp;
	break;
      	
      case 2: //1,2,0
	if(d1[0]*d1[1] >0)
	  nb31 = nb13 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = nb13 = 1;
	if(nb31 < 1) nb31 = 1;
	if(nb13 < 1) nb13 = 1;

	if(OW && dt2 <= dt1)
	  inplace = true;
	else
	  tmp = new Type2[d2[2]*d2[1]*nb31];

	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  if(inplace) {
	    tmp = (Type2 *) (in+k*d1[0]*d1[1]);
	    (*(trans_type->exec))(plan->libplan_in,in+k*d1[0]*d1[1],tmp);
	  }
	  else
	    (*(trans_type->exec))(plan->libplan_out,pIN+k*d1[0]*d1[1],tmp);

	  for(i=0;i < d2[1];i+=nb13) {
	    i2 = min(i+nb13,d2[1]);
	    pout2 =  out + i*d2[0];
	    pin2 =  tmp + i;
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) pin2;
	      pout1 = (float *) pout2+ds*kk;
	      for(j=0;j < d1[1];j++) {
		pin = pin1  ;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d2[0]-1);
		}
		pin1 += ds*d2[1];
		pout1+= ds*d2[0]*d2[1];
	      }
	      pin2 += d2[1]*d2[2];
	    }
	  }
	}
	if(!inplace)
	  delete [] tmp;
	
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	if(d1[0]*d1[1] >0)
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d2[0]*d2[1] >0)
	  nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else
	  nb13 = 1;
	if(nb13 < 1) nb13 = 1;

	tmp = new Type2[d2[2]*d2[1]*d2[0]];
	(*(trans_type->exec))(plan->libplan_out,pIN,tmp);
	
	for(k=0;k <d2[0];k+=nb31) {
	  k2 = min(k+nb31,d2[0]);
	  for(i=0;i < d2[2];i+=nb13) {
	    i2 = min(i+nb13,d2[2]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) tmp + ds*(kk*d2[2]*d2[1] +i);
	      pout1 = (float *) out + ds*(kk + i *d2[1]*d2[0]);
	      for(j=0;j < d2[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d2[0]*d2[1]-1);
		}
		pin1 += ds*d2[2];
		pout1 += ds*d2[0];
	      }
	    }
	  }
	}

	delete [] tmp;
	break;

      
      case 0: //2,0,1
	if(d2[0]*d2[1] > 0)
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else
	  nb23 = 1;
	if(nb23 < 1) nb23 = 1;
	if(d2[0]>0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d2[0]*nb23);
	else
	  nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	
	tmp = new Type2[d2[0]*d2[1]*d2[2]];
	(*(trans_type->exec))(plan->libplan_out,pIN,tmp);
	
	for(k=0;k <d1[1];k+=nb32) {
	  k2 = min(k+nb32,d1[1]);
	  for(j=0;j < d1[2];j+=nb23) {
	    j2 = min(j+nb23,d1[2]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = (float *) tmp +ds*(kk*d2[2] +j*d2[0]*d2[2]);
	      pout1 = (float *) out +ds*(kk +j*d2[0]);
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d2[2];i++) {
		  for(int m=0;m<ds;m++)
		    *pout++ =  *pin++;
		  pout += ds*(d2[0]*d2[1]-1);
		}
		pin1 += ds*d2[2]*d2[0];
		pout1+= ds*d2[0];
	      }
	    }
	  }
	}
	
	delete [] tmp;
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {
	if(d1[0]*d1[1] > 0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d1[0]*d1[1] > 0)
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;

	if((void *) in != (void *) out) {

	  for(k=0;k <d1[2];k+=nb32) {
	    k2 = min(k+nb32,d1[2]);
	    for(j=0;j < d1[1];j+=nb23) {
	      j2 = min(j+nb23,d1[1]);
	      for(kk=k; kk < k2; kk++) {
		for(jj=j; jj < j2; jj++) {
		  pin_t1 = pIN +kk*d1[0]*d1[1] +jj*d1[0];
		  pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
		  (*(trans_type->exec))(plan->libplan_out,pin_t1,pout2);
		}
	      }
	    }
	  }
	}
	else { // In-place
	  tmp = new Type2[d2[0]*d2[1]*d2[2]];
	  (*(trans_type->exec))(plan->libplan_out,pIN,tmp);

	  for(k=0;k <d2[1];k+=nb32) {
	    k2 = min(k+nb32,d2[1]);
	    for(j=0;j < d2[2];j+=nb23) {
	      j2 = min(j+nb23,d2[2]);
	      for(kk=k; kk < k2; kk++) {
		for(jj=j; jj < j2; jj++) {
		  pin2 = tmp +kk*d2[0]*d2[2] +jj*d2[0];
		  pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
		  for(i=0;i<d2[0];i++)
		    *pout2++ = *pin2++;
		}
	      }
	    }
	  }
	  delete [] tmp;
	  
	}
      }	
      break;
   }
    
#ifdef FFTW
   if(dt1 > dt2 && !OW)
     delete [] pIN;
#endif

  }
  else { // Scheme is transform after reordering

    Type1 *tmp,*pin2,*pin3,*pout3,*pout4;
    float *pin,*pout,*pin1,*pout1,*ptran;
    Type2 *pout2,*ptran2;
    int ds=sizeof(Type1)/4;

    switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2

	if((void *) in == (void *) out && dt2 > dt1) {

	tmp = new Type1[d1[0]*d1[1]*d1[2]];
	for(k=0;k <d1[2];k++) {
	  pin = (float *) in + ds*k*d1[0]*d1[1];
#ifdef MKL_BLAS
	  blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmp,d1[1]);
#else
	  pout1 = (float *) tmp +ds*k*d1[0]*d1[1];
	  for(j=0;j < d1[1];j++) {
	    pout = pout1 + ds*j;
	    for(i=0;i < d1[0];i++) {
	      for(int m=0;m<ds;m++)
		*pout++ = *pin++;
	      pout += ds*(d1[1]-1);
	    }	
	  }
#endif
	  pout2 = out + d2[0]*d2[1]*k;
	  (*(trans_type->exec))(plan->libplan_out,(Type1 *) pout1, pout2);
	}

	}
	else {
	tmp = new Type1[d1[0]*d1[1]];
	for(k=0;k <d1[2];k++) {
	  pin = (float *) in + ds*k*d1[0]*d1[1];
#ifdef MKL_BLAS
	  blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmp,d1[1]);
#else
	  for(j=0;j < d1[1];j++) {
	    pout = (float *) tmp +ds*j;
	    for(i=0;i < d1[0];i++) {
	      for(int m=0;m<ds;m++)
		*pout++ = *pin++;
	      pout += ds*(d1[1]-1);
	    }	
	  }
#endif
	  pout2 = out + k*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out,tmp,pout2);
	}

	}
	delete [] tmp;

	break;

      case 2: //1,2,0
	if(d1[0]*d1[1] > 0)
	  nb31 = nb13 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = nb13 = 1;
	if(nb31 < 1){  nb31 = 1; nb13 = 1; }

	if((void *) in == (void *) out) {

	tmp = new Type1[d1[0]*d1[1]*d1[2]]; // tmp[k][i][j] = in[i][j][k]

	for(j=0;j < d1[1];j++) {
	  pin1 = (float *) in +ds*j*d1[0];
	  pout1 =(float *)  tmp + ds*j*d1[2]*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + ds*(i + kk*d1[0]*d1[1]);
		pout = pout1 + ds*(i*d1[2] +kk);
		
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d1[2]-1);
		}
	      }
	    }
	  }
	}
	for(j=0;j < d1[1];j++) {
	  pout1 =(float *)  tmp + ds*j*d1[2]*d1[0];
	  pout2 = out + j*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out,(Type1 *) pout1,pout2);
	}
	}
	else {
	tmp = new Type1[d1[0]*d1[2]]; // tmp[k][i] = in[i][j][k]

	for(j=0;j < d1[1];j++) {
	  pin1 = (float *) in +ds*j*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + ds*(i + kk*d1[0]*d1[1]);
		pout = (float *)  tmp + ds*(i*d1[2] +kk);
		
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d1[2]-1);
		}
	      }
	    }
	  }

	  pout2 = out + j*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out,tmp,pout2);
	}
	
	}	
	delete [] tmp;

       	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	if(d1[0]*d1[1] >0)
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d2[0]*d2[1] >0)
	  nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb13 = 1;
	if(nb13 < 1) nb13 = 1;
	
	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  for(i=0;i < d1[0];i+=nb13) {
	    i2 = min(i+nb13,d1[0]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) in + ds*(kk*d1[0]*d1[1] +i);
	      pout1 = (float *) tmp + ds*(kk +i*d1[2]*d1[1]);
	      for(j=0;j < d1[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d1[2]*d1[1]-1);
		}
		pin1 += ds*d1[0];
		pout1 += ds*d1[2];
	      }
	    }
	  }
	}
	//      if(!trans_type->is_empty)
	(*(trans_type->exec))(plan->libplan_out,tmp,out);
	delete [] tmp;
	
	break;
      case 0: //2,0,1
	if(d1[0]*d1[1] >0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb32 < 1) nb32 = 1;
	nb23 = nb32;
	
	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[1];k+=nb32) {
	  k2 = min(k+nb32,d1[1]);
	  for(j=0;j < d1[2];j+=nb23) {
	    j2 = min(j+nb23,d1[2]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = (float *) in + ds*(kk*d1[0] +j*d1[0]*d1[1]);
	      pout1 = (float *) tmp + ds*(kk +j*d1[1]);
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d1[0];i++) {
		  for(int m=0;m<ds;m++)
		    *pout++ =  *pin++;
		  pout += ds*(d1[1]*d1[2]-1);
		}
		pin1 += ds*d1[0]*d1[1];
		pout1 += ds*d1[1];
	      }
	    }
	  }
	}
	
	(*(trans_type->exec))(plan->libplan_out,tmp,out);
	delete [] tmp;
	
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {

	if(d1[0]*d1[1] >0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d2[0]*d2[1] >0)
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;
	
	if((void *) in == (void *) out) {

	tmp = new Type1[d1[0]*d1[1]*d1[2]]; // tmp[k][i] = in[i][j][k]
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      pin2 = in + kk*d1[0]*d1[1] +j*d1[0];
	      pout3 = tmp +kk*d2[0] +j*d2[0]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		pout4 = pout3;
		pin3 = pin2;
		for(i=0;i<d1[0];i++) 
		  *pout4++ = *pin3++;
		pin2 += d1[0];
		pout3 += d2[0]*d2[1];
	      }
	    }
	  }
	}

	for(k=0;k<d2[2];k++)
	  for(j=0;j<d2[1];j++) {
	    pin2 = tmp+d1[0]*(j+k*d2[1]);
	    pout2 = out + d2[0]*(j+k*d2[1]);
	    if(trans_type->exec != NULL)
	      (*(trans_type->exec))(plan->libplan_out,pin2,pout2);
	  }

	}
	else {
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      pin2 = in + kk*d1[0]*d1[1] +j*d1[0];
	      ptran2 = out +kk*d2[0] +j*d2[0]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		if(trans_type->exec != NULL)
		  (*(trans_type->exec))(plan->libplan_out,pin2,ptran2);
		pin2 += d1[0];
		ptran2 += d2[0]*d2[1];
	      }
	    }
	  }
	}
	}	
	break;
      }
      
    }
  }

#ifdef DEBUG
  sprintf(str,"reorder_trans.out%d",cnt_reorder_trans++);
  write_buf<Type2>(out,str,dims2,imo2);
#endif  
}

// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]
// Reorder (local transpose), perform 1D transform followed by spectral derivative
  template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_deriv(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1, bool OW)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2,dims[3];
  void rel_change(int *,int *,int *);

  int imo1[3],imo2[3],d1[3],d2[3];
  int scheme,nb13,nb31,nb23,nb32;
  
#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  printf("%d: In transplan::reorder_deriv, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  inv_mo(mo1,imo1); // Find inverse storage dimensions mappings
  inv_mo(mo2,imo2);
//Fins local storage dimensions, i.e. in[d1[0]][d1[1]]d1[2]], out[d2[0]][d2[1]][d2[2]]
  for(i=0;i<3;i++) { 
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }
  
  /*
    if(dims1[trans_dim] != *inembed) {
      printf("Error in reorder_deriv: leading dimension on input %d doesn't match transform type %d\n",dims1[trans_dim],*inembed);
      return;
    }
    if(dims2[trans_dim] != *onembed) {
      printf("Error in reorder_deriv: leading dimension on output %d doesnt match transform type %d\n", dims2[trans_dim],*onembed);
      return;
    }
  */

  if(mo1[trans_dim] == 0) 
    scheme = TRANS_IN;
  else if(mo2[trans_dim] == 0) 
    scheme = TRANS_OUT;
  else {
    printf("Error in reorder_deriv: expected dimension %d to be the leading dimension for input or output\n",trans_dim);
    return;
   }

  rel_change(imo1,imo2,mc);

#ifdef DEBUG
  printf("In reorder_deriv, mc=%d %d %d, scheme=%s\n",mc[0],mc[1],mc[2],(scheme == TRANS_IN) ? "IN" : "OUT"); 
  char str[80];
  static int cnt_reorder_trans=0;
  /*
  sprintf(str,"reorder_deriv.in%d",cnt_reorder_trans);
  write_buf<Type1>(in,str,dims1,imo1);
  */
#endif

  if(scheme == TRANS_IN) { // Scheme is transform before reordering

    Type1 *pin_t1;
    Type2 *tmp,*pin2,*pout2;
    float *pin,*pin1,*pout,*pout1,*ptran2; 
    int ds=sizeof(Type2)/4;

   switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
;
	tmp = new Type2[d2[0]*d2[1]];
	dims[0] = d2[1];
	dims[1] = d2[0];
	dims[2] = 1;

	for(k=0;k <d1[2];k++) {
	  (*(trans_type->exec))(plan->libplan_out,in+k*d1[0]*d1[1],tmp);
	  compute_deriv_loc(tmp,tmp,dims);
	  pout = (float *) out+ds*k*d2[0]*d2[1];
#ifdef MKL_BLAS
#ifdef DEBUG
	  if(taskid == 0) {
	    printf("A=\n");
	    for(i=0;i<d2[0];i++) {
	      for(j=0;j<d2[1];j++)
		printf("%lf ",*(tmp+i*d2[1]+j));
	      printf("\n");
	    }
	  }
#endif
	  blas_trans<Type2>(d2[1],d2[0],1.0,tmp,d2[1],pout,d2[0]);
#ifdef DEBUG
	  if(taskid == 0) {
	    printf("B=\n");
	    for(i=0;i<d2[1];i++) {
	      for(j=0;j<d2[0];j++)
		printf("%lf ",*(pout+i*d2[0]+j));
	      printf("\n");
	    }
	  }
#endif
	
#else

	  pout = (float *) out +ds*k*d2[0]*d2[1];
	  for(j=0;j < d2[1];j++) {
	    pin = (float *) tmp + ds*j;
	    for(i=0;i < d2[0];i++) {
	      for(int m=0;m<ds;m++)
		*pout++ = *pin++;
	      pin += ds*(d2[1]-1);
	    }	
	 
	  }
#endif
	  //      if(!trans_type->is_empty)
	}
	
	delete [] tmp;
	break;
	
      case 2: //1,2,0
	if(d1[0]*d1[1] >0)
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	nb13 = nb31;

	// d1[0] -> d2[1]
	// d1[1] -> d2[2]
	// d1[2] -> d2[0]
	tmp = new Type2[d2[2]*d2[1]*nb31];
	dims[0] = d2[1];
	dims[1] = d2[2];
	dims[2] = nb13;

	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  (*(trans_type->exec))(plan->libplan_out,in+k*d1[0]*d1[1],tmp);
	  compute_deriv_loc(tmp,tmp,dims);
	  for(i=0;i < d2[1];i+=nb13) {
	    i2 = min(i+nb13,d2[1]);
	    pout2 =  out + i*d2[0];
	    pin2 =  tmp + i;
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) pin2;
	      pout1 = (float *) pout2+ds*kk;
	      for(j=0;j < d1[1];j++) {
		pin = pin1  ;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d2[0]-1);
		}
		pin1 += ds*d2[1];
		pout1+= ds*d2[0]*d2[1];
	      }
	      pin2 += d2[1]*d2[2];
	    }
	  }
	}
	delete [] tmp;
	
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	if(d1[0]*d1[1] >0)
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d1[0]*d1[1] >0)
	  nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb13 = 1;
	if(nb13 < 1) nb13 = 1;

	// d1[0] -> d2[2]
	// d1[1] -> d2[1]
	// d1[2] -> d2[0]
	dims[0] = d2[2];
	dims[1] = d2[1];
	dims[2] = d2[0];

	tmp = new Type2[d2[0]*d2[1]*d2[2]];

	(*(trans_type->exec))(plan->libplan_out,in,tmp);
	compute_deriv_loc(tmp,tmp,dims);
	
	for(k=0;k <d2[0];k+=nb31) {
	  k2 = min(k+nb31,d2[0]);
	  for(i=0;i < d2[2];i+=nb13) {
	    i2 = min(i+nb13,d2[2]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) tmp + ds*(kk*d2[2]*d2[1] +i);
	      pout1 = (float *) out + ds*(kk + i *d2[1]*d2[0]);
	      for(j=0;j < d2[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d2[0]*d2[1]-1);
		}
		pin1 += ds*d2[2];
		pout1 += ds*d2[0];
	      }
	    }
	  }
	}
	
	delete [] tmp;
	break;
      case 0: //2,0,1
	if(d2[0]*d2[1] >0)	
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;
	if(d2[0] >0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d2[0]*nb23);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	// d1[0] -> d2[2]
	// d1[1] -> d2[0]
	// d1[2] -> d2[1]
	dims[0] = d2[2];
	dims[1] = d2[0];
	dims[2] = d2[1];

	tmp = new Type2[d2[0]*d2[1]*d2[2]];

	(*(trans_type->exec))(plan->libplan_out,in,tmp);
	compute_deriv_loc(tmp,tmp,dims);
	
	for(k=0;k <d1[1];k+=nb32) {
	  k2 = min(k+nb32,d1[1]);
	  for(j=0;j < d1[2];j+=nb23) {
	    j2 = min(j+nb23,d1[2]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = (float *) tmp +ds*(kk*d2[2] +j*d2[0]*d2[2]);
	      pout1 = (float *) out +ds*(kk +j*d2[0]);
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d2[2];i++) {
		  for(int m=0;m<ds;m++)
		    *pout++ =  *pin++;
		  pout += ds*(d2[0]*d2[1]-1);
		}
		pin1 += ds*d2[2]*d2[0];
		pout1+= ds*d2[0];
	      }
	    }
	  }
	}
	
	delete [] tmp;
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {
	if(d1[0]*d1[1] >0)
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d2[0]*d2[1] >0)
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;

	// d1[0] -> d2[0]
	// d1[1] -> d2[2]
	// d1[2] -> d2[1]
	dims[0] = d2[0];
	dims[1] = 1;
	dims[2] = 1;
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      for(jj=j; jj < j2; jj++) {
		pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
		pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
		(*(trans_type->exec))(plan->libplan_out,pin_t1,pout2);
		compute_deriv_loc(pout2,pout2,dims);
	      }
	    }
	  }
	}    
      }
      
      break;
    }
    
  }
  else { // Scheme is transform after reordering

    Type1 *tmp,*pin2;
    float *pin,*pout,*pin1,*pout1,*ptran;
    Type2 *pout2,*ptran2;
    int ds=sizeof(Type1)/4;

    switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
	  
	// d1[0] -> d2[1]
	// d1[1] -> d2[0]
	// d1[2] -> d2[2]
	dims[0] = d2[0];
	dims[1] = d2[1];
	dims[2] = 1;

	tmp = new Type1[d1[0]*d1[1]];
	for(k=0;k <d1[2];k++) {
	  pin = (float *) in + ds*k*d1[0]*d1[1];
#ifdef MKL_BLAS
	  blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmp,d1[1]);
#else
	  for(j=0;j < d1[1];j++) {
	    pout = (float *) tmp +ds*j;
	    for(i=0;i < d1[0];i++) {
	      for(int m=0;m<ds;m++)
		*pout++ = *pin++;
	      pout += ds*(d1[1]-1);
	    }	
	  }
#endif
	  pout2 = out + k*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out,tmp,pout2);
	  compute_deriv_loc(pout2,pout2,dims);
	}
	delete [] tmp;

	break;

      case 2: //1,2,0
	if(d1[0]*d1[1] >0)	
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	nb13 = nb31;

	// d1[0] -> d2[1]
	// d1[1] -> d2[2]
	// d1[2] -> d2[0]
	dims[0] = d2[0];
	dims[1] = d2[1];
	dims[2] = 1;

	tmp = new Type1[d1[0]*d1[2]]; // tmp[k][i] = in[i][j][k]

	for(j=0;j < d1[1];j++) {
	  pin1 = (float *) in +ds*j*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + ds*(i + kk*d1[0]*d1[1]);
		pout = (float *)  tmp + ds*(i*d1[2] +kk);
		
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d1[2]-1);
		}
	      }
	    }
	  }

	  pout2 = out + j*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out,tmp,pout2);
          compute_deriv_loc(pout2,pout2,dims);
	}
	
	delete [] tmp;
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	if(d1[0]*d1[1] >0)	
	  nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d2[0]*d2[1] >0)	
	  nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb13 = 1;
	if(nb13 < 1) nb13 = 1;
	
	// d1[0] -> d2[2]
	// d1[1] -> d2[1]
	// d1[2] -> d2[0]
	dims[0] = d2[0];
	dims[1] = d2[1];
	dims[2] = d2[2];

	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  for(i=0;i < d1[0];i+=nb13) {
	    i2 = min(i+nb13,d1[0]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = (float *) in + ds*(kk*d1[0]*d1[1] +i);
	      pout1 = (float *) tmp + ds*(kk +i*d1[2]*d1[1]);
	      for(j=0;j < d1[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  for(int m=0;m<ds;m++)
		    *pout++ = *pin++;
		  pout += ds*(d1[2]*d1[1]-1);
		}
		pin1 += ds*d1[0];
		pout1 += ds*d1[2];
	      }
	    }
	  }
	}
	//      if(!trans_type->is_empty)
	(*(trans_type->exec))(plan->libplan_out,tmp,out);
	compute_deriv_loc(out,out,dims);
	delete [] tmp;
	
	break;
      case 0: //2,0,1
	if(d1[0]*d1[1] >0)	
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	nb23 = nb32;

	// d1[0] -> d2[2]
	// d1[1] -> d2[0]
	// d1[2] -> d2[1]
	dims[0] = d2[0];
	dims[1] = d2[1];
	dims[2] = d2[2];
	
	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[1];k+=nb32) {
	  k2 = min(k+nb32,d1[1]);
	  for(j=0;j < d1[2];j+=nb23) {
	    j2 = min(j+nb23,d1[2]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = (float *) in + ds*(kk*d1[0] +j*d1[0]*d1[1]);
	      pout1 = (float *) tmp + ds*(kk +j*d1[1]);
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d1[0];i++) {
		  for(int m=0;m<ds;m++)
		    *pout++ =  *pin++;
		  pout += ds*(d1[1]*d1[2]-1);
		}
		pin1 += ds*d1[0]*d1[1];
		pout1 += ds*d1[1];
	      }
	    }
	  }
	}
	
	(*(trans_type->exec))(plan->libplan_out,tmp,out);
	compute_deriv_loc(out,out,dims);
	delete [] tmp;
	
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {
	if(d1[0]*d1[1] >0)	
	  nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d2[0]*d2[1] >0)	
	  nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;

	// d1[0] -> d2[0]
	// d1[1] -> d2[2]
	// d1[2] -> d2[1]
	dims[0] = d2[0];
	dims[1] = 1;
	dims[2] = 1;
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      pin2 = in + kk*d1[0]*d1[1] +j*d1[0];
	      ptran2 = out +kk*d2[0] +j*d2[0]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		if(trans_type->exec != NULL) {
		  (*(trans_type->exec))(plan->libplan_out,pin2,ptran2);
		  compute_deriv_loc(ptran2,ptran2,dims);
		}
		pin2 += d1[0];
		ptran2 += d2[0]*d2[1];
	      }
	    }
	  }
	}
	
	break;
      }
      
    }
  }

#ifdef DEBUG
  sprintf(str,"reorder_deriv.out%d",cnt_reorder_trans++);
  write_buf<Type2>(out,str,dims2,imo2);
#endif  
}

  // Reorder (local transpose) out-of-place
// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// dims2[mo2[i]] = dims1[mo1[i]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int *dims_init)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);
  float *pin,*pout,*pin1,*pout1;
  int imo1[3],imo2[3],d1[3],d2[3],nb13,nb31,nb23,nb32;
  int ds=sizeof(Type2) /4;

  // Inverse storage mapping
  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  // Find local storage dimensions
  for(i=0;i<3;i++) {
    d1[i] = dims_init[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }

  // Find relative change in the memory mapping from input to output
  rel_change(imo1,imo2,mc);
  pin = (float *) in;
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2
      for(k=0;k <d1[2];k++) {
	pout1 = (float *) out + ds*k*d1[0]*d1[1];
	for(j=0;j < d1[1];j++)
	  pout = pout1 + j;
	  for(i=0;i < d1[0];i++) {
	    for(int m=0;m<ds;m++)
	      *(pout++)  = *(pin++);
	    pout += ds*(d2[0]-1);
	  }
      }
      break;

    case 2: //1,2,0
      if(d1[0]*d1[1] >0)	
	nb31 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      else nb31 = 1;
      if(nb31 < 1) nb31 = 1;
      nb13 = nb31;

      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = (float *) in +ds*(kk*d1[0]*d1[1] +i);
	    pout1 = (float *) out + ds*(kk*d2[1] +i*d2[0]*d2[1]);
	    for(j=0;j < d1[1];j++) {
	      pin = pin1;
	      pout = pout1;
	      for(ii=i; ii < i2; ii++) {
		for(int m=0;m<ds;m++)
		  *pout++ = *pin++;
		pout += ds*(d2[0]*d2[1]-1);
	      }
	      pin1 += ds*d1[0];
	      pout1+=ds;
	    }
	  }
	}
      }

      break;
    }
  
    break;
  case 2:
    switch(mc[1]) {
    case 1: //2,1,0
      if(d1[0]*d1[1] >0)	
	nb31 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      else nb31 = 1;
      if(nb31 < 1) nb31 = 1;
      nb13 = nb31;

      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = (float *) in + ds*(kk*d1[0]*d1[1] +i);
	    pout1 = (float *) out + ds*(kk + i * d2[0]*d2[1]);
	    for(j=0;j < d1[1];j++) {
	      pin = pin1;
	      pout = pout1;
	      for(ii=i; ii < i2; ii++) {
		for(int m=0;m<ds;m++)
		  *pout++ = *pin++;
		pout += ds*(d2[0]*d2[1]-1);
	      }
	      pin1 += ds*d1[0];
	      pout1 += ds*d2[0];
	    }
	  }
	}
      }
      
      break;
    case 0: //2,0,1
      if(d1[0]*d1[1] >0)	
	nb32 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      else nb32 = 1;
      if(nb32 < 1) nb32 = 1;
      if(d2[0]*d2[1] >0)	
	nb23 = CACHE_BL / (sizeof(Type2)*d2[0]*d2[1]);
      else nb23 = 1;
      if(nb23 < 1) nb23 = 1;

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = (float *) in +ds*kk*d1[0]*d1[1];
	    pout1 = (float *) out +ds*kk;
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1 +ds*jj*d1[0];
	      pout = pout1 +ds*jj*d2[0]*d2[1];
	      for(i=0;i < d1[0];i++) {
		for(int m=0;m<ds;m++)
		  *pout++ = *pin++;
		pout += ds*(d2[0]-1);
	      }
	    }
	  }
	}
      }

      break;
    }
    break;
  case 0: //0,2,1
    if(mc[1] == 2) {
      if(d1[0]*d1[1] >0)	
	nb32 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      else nb32 = 1;
      if(nb32 < 1) nb32 = 1;
      if(d2[0]*d2[1] >0)	
	nb23 = CACHE_BL / (sizeof(Type2)*d2[0]*d2[1]);
      else nb23 = 1;
      if(nb23 < 1) nb23 = 1;

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = (float *) in + ds*(kk*d1[0]*d1[1] +j*d1[0]);
	    pout1 = (float *) out +ds*(kk*d2[0] +j*d2[0]*d2[1]);
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) 
		for(int m=0;m<ds;m++)
		  *pout++ = *pin++;
	      pin += ds*d1[0];
	      pout += ds*d2[0]*d2[1];
	    }
	  }
	}
      }      
    }
    else
      for(k=0;k<d1[2];k++)
	for(j=0;j < d1[1];j++)
	  for(i=0;i < d1[0];i++)
	    *out++ = *in++;
    break;
  }


}

  // Reorder (local transpose) in-place (out=in)
// Input: in[d1[imo1[0]]][d1[imo1[1]]][d1[imo1[2]]]
// dims2[mo2[i]] = d1[mo1[i]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]
// Assume the 'in' buffer is large enough for both input and output
template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int dims[3])
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  Type *pin,*pin1,*pout,*pout1,*tmp;
  void rel_change(int *,int *,int *);
  int d1[3];

  for(i=0;i<3;i++)
    d1[mo1[i]] = dims[i];

  // Find relative change in memory mapping from input to output
  rel_change(mo1,mo2,mc);

  int nb13,nb31,nb23,nb32;
  int pad = CACHEPAD/sizeof(Type);
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2
      tmp = new Type[(d1[0]+pad)*d1[1]];
      pin1 = in;
      for(k=0;k <d1[2];k++) {
	pout = tmp;
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout++  = *pin1++;
	  pout+= pad;//Cache shift
	}
	pin = tmp;
	pout1 = in +k*d1[0]*d1[1];
	for(j=0;j < d1[1];j++) {
	  pout = pout1;
	  for(i=0;i < d1[0];i++) {	    
	    *pout  = *pin++;
	    pout += d1[1];
	  }
	  pin += pad;
	  pout1++;
	}
      }
      delete [] tmp;

      break;

    case 2: //1,2,0
      if(d1[0]*d1[1] >0)	
	nb31 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      else nb31 = 1;
      if(nb31 < 1) nb31 = 1;
      if(d1[2]*d1[1] >0)	
	nb13 = CACHE_BL / (sizeof(Type)*d1[2]*d1[1]);
      else nb13 = 1;
      if(nb13 < 1) nb13 = 1;

      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pin = in;
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout++  = *pin++;
	  pout++;//Cache shift
	}
      
      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*(d1[0]+1)*d1[1];
	    pout1 = in + kk*d1[1];
	    for(j=0;j < d1[1];j++) {
	      pin = pin1+i;
	      pout = pout1+i*d1[2]*d1[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += d1[2]*d1[1];
	      }
	      pin1 += d1[0]+1;
	      pout1++;
	    }
	  }
	}
      }
      delete [] tmp;

      break;
    }
  
    break;
  case 2:
    switch(mc[1]) {
    case 1: //2,1,0
      if(d1[0]*d1[1] >0)	
	nb31 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      else nb31 = 1;
      if(nb31 < 1) nb31 = 1;
      if(d1[2]*d1[1] >0)	
	nb13 = CACHE_BL / (sizeof(Type)*d1[2]*d1[1]);
      else nb13 = 1;
      if(nb13 < 1) nb13 = 1;

      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pin = in;
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout++  = *pin++;
	  pout++;//Cache shift
	}

      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp + kk*(d1[0]+1)*d1[1];
	    pout1 = in + kk;
	    for(j=0;j < d1[1];j++) {
	      pin = pin1 + i;
	      pout = pout1 + i * d1[1]*d1[2];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += d1[2]*d1[1];
	      }
	      pin1 += d1[0]+1;
	      pout1 += d1[2];
	    }
	  }
	}
      }

      delete [] tmp;
      
      break;
    case 0: //2,0,1
      if(d1[0]*d1[1] >0)	
	nb32 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      else nb32 = 1;
      if(nb32 < 1) nb32 = 1;
      if(d1[0]*d1[2] >0)	
	nb23 = CACHE_BL / (sizeof(Type)*d1[0]*d1[2]);
      else nb23 = 1;
      if(nb23 < 1) nb23 = 1;

      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pin = in;
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++){
	  for(i=0;i < d1[0];i++) 
	    *pout++  = *pin++;
	  pout++;//Cache shift
	}

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = tmp +kk*d1[0]*d1[1] +j*d1[0];
	    pout1 = in +kk +j*d1[0]*d1[2];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) {
		*pout = *pin++;
		pout += d1[2];
	      }
	      pin1 += d1[0]+1;
	      pout1 += d1[0]*d1[2];
	    }
	  }
	}
      }

      delete [] tmp;
      break;
    }
    break;
  case 0: //0,2,1
    if(mc[1] == 2) {
      if(d1[0]*d1[1] >0)	
	nb32 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      else nb32 = 1;
      if(nb32 < 1) nb32 = 1;
      if(d1[0]*d1[2] >0)	
	nb23 = CACHE_BL / (sizeof(Type)*d1[0]*d1[2]);
      else nb23 = 1;
      if(nb23 < 1) nb23 = 1;

      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pin = in;
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout++  = *pin++;
	  pout++; //Cache shift
	}

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*d1[0]*d1[1] +j*d1[0];
	    pout1 = in +kk*d1[0] +j*d1[0]*d1[2];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) 
		*pout++ = *pin++;
	      pin += d1[0]+1;
	      pout += d1[0]*d1[2];
	    }
	  }
	}
      }      
      delete [] tmp;

    }
    break;
  }
}


// mo2[mc[i]] = mo1[i]
void  rel_change(int mo1[3],int mo2[3],int mc[3])
{
  if(mo2[2] !=mo1[2]) 
    if(mo2[1] == mo1[2])
      mc[2]=1;
    else
      mc[2]=0;
  else
    mc[2]=2;

  if(mo2[1] != mo1[1])
    if(mo2[0] == mo1[1])
      mc[1]=0;
    else
      mc[1]=2;
  else
    mc[1]=1;
  
  if(mo2[0] != mo1[0])
    if(mo2[1] == mo1[0])
      mc[0]=1;
    else
      mc[0] = 2;
  else
    mc[0] = 0;
  
}

// Perform MPI exchange
template <class Type> void MPIplan<Type>::exec(char *in_,char *out_) {
  Type *in,*out;
  in = (Type *) in_;
  out = (Type *) out_;
  Type *sendbuf = new Type[dims1[0]*dims1[1]*dims1[2]];
#ifdef TIMERS
  double t1=MPI_Wtime();
#endif
  pack_sendbuf(sendbuf,(Type *) in);
#ifdef TIMERS
      timers.packsend += MPI_Wtime() -t1;
#endif
  Type *recvbuf = new Type[dims2[0]*dims2[1]*dims2[2]];
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  MPI_Alltoallv(sendbuf,SndCnts,SndStrt,MPI_REAL,recvbuf,RcvCnts,RcvStrt,MPI_REAL,mpicomm);
#ifdef TIMERS
      timers.alltoall += MPI_Wtime() -t1;
#endif
  delete [] sendbuf;
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  unpack_recvbuf((Type *) out,recvbuf);
#ifdef TIMERS
      timers.unpackrecv += MPI_Wtime() -t1;
#endif

  delete [] recvbuf;
}

template <class Type> void MPIplan<Type>::pack_sendbuf(Type *sendbuf,Type *src)
{
  int i,x,y,z;
  float *p1,*p0;

  int j,istart[numtasks][3],iend[numtasks][3],d[3];
  int imo1[3],imo2[3];
  int ds=sizeof(Type)/4;

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);

  for(i=0;i<numtasks;i++)
    for(j=0;j<3;j++) {
      istart[i][j] = 0;
      iend[i][j] = dims1[j];
      //      istart[i][j] = grid1->st[comm_id][i][j]; iend[i][j] = grid1->en[comm_id][i][j];
    }
  for(j=0;j<numtasks;j++) {
    istart[j][d2] = grid2->st[comm_id][j][d2]; iend[j][d2] = grid2->en[comm_id][j][d2];
  }    


  for(i=0;i<3;i++)
    d[i] = dims1[imo1[i]];
  
  for(i=0;i < numtasks;i++) {
    p1 = (float *) sendbuf + *(SndStrt+i);
    for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++)
      for(y=istart[i][imo1[1]];y< iend[i][imo1[1]];y++) {
	p0 = (float *) src+ds*(d[0]*(z*d[1]+y) +istart[i][imo1[0]]);
	for(x=ds*istart[i][imo1[0]];x < ds*iend[i][imo1[0]];x++)
	  *p1++ = *p0++;;
      }
  }

}

template <class Type> void MPIplan<Type>::unpack_recvbuf(Type *dest,Type *recvbuf)
{
  int i,ii,x,y,z,k,x2,y2,z2;
  int ds=sizeof(Type)/4;
  float *p1,*p0,*pin,*pin1,*pout,*pout1;

  int j,istart[numtasks][3],iend[numtasks][3],isize[numtasks][3],d[3];
  int xstart,xend,ystart,yend,zstart,zend;
  int xsize,ysize,zsize;
  int imo1[3],imo2[3],nb13,nb31,nb23,nb32;

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);

  for(i=0;i<numtasks;i++) 
    for(j=0;j<3;j++) {
      istart[i][j] = 0;
      isize[i][j] = iend[i][j] = dims2[j];
      //       isize[i][j] = grid2->sz[comm_id][i][j];
      //istart[i][j] = grid2->st[comm_id][i][j]; iend[i][j] = grid2->en[comm_id][i][j];
    }
  for(j=0;j<numtasks;j++) {
    isize[j][d1] = grid1->sz[comm_id][j][d1];
    istart[j][d1] = grid1->st[comm_id][j][d1]; iend[j][d1] = grid1->en[comm_id][j][d1];
  }    
  /*
  for(i=0;i<numtasks;i++) 
    for(j=0;j<3;j++) {
      printf("%d: istart[%d][%d]=%d\n",taskid,i,j,istart[i][j]);
      printf("%d: iend[%d][%d]=%d\n",taskid,i,j,iend[i][j]);
      printf("%d: grid1.en[%d][%d]=%d\n",taskid,i,j,grid1->en[comm_id][i][j]);
      printf("%d: grid2.en[%d][%d]=%d\n",taskid,i,j,grid2->en[comm_id][i][j]);
      printf("%d: dims2[%d]=%d\n",taskid,i,dims2[i]);
    }
  */

  for(i=0;i<3;i++)
    d[i] = dims2[imo2[i]];

  if(mo1[0] == mo2[0] && mo1[1] == mo2[1] && mo1[2] == mo2[2]) 

    for(i=0;i < numtasks;i++) {
      p1 = (float *) recvbuf + *(RcvStrt+i);
      for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++)
	for(y=istart[i][imo1[1]];y < iend[i][imo1[1]];y++) {
	  p0 = (float *) dest + ds*(d[0]*(z*d[1]+y) + istart[i][imo1[0]]);
	  for(x=istart[i][imo1[0]];x < iend[i][imo1[0]];x++)
	    for(j=0;j<ds;j++)
	      *p0++ = *p1++;
	}
    }
  else {
    int mc[3];
    rel_change(imo1,imo2,mc);

    //    printf("mo1=%d %d %d, mo2=%d %d %d, mc=%d %d %d\n",mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2],mc[0],mc[1],mc[2]);
    
    switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
	for(i=0;i < numtasks;i++) {
	  p1 = (float *) recvbuf + *(RcvStrt+i);
	  p0 = (float *) dest + ds*(istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]])); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];
	  for(z=0;z < zsize;z++) 
	    for(y=0;y < ysize;y++) {
	      pout = p0 + ds*(d[0]*d[1]*z +y);
	      for(x=0;x < xsize;x++) {
		for(j=0;j<ds;j++)
		  *pout++ = *p1++;
		pout += ds*(d[0]-1);
	      }
	    }
	}
	
	break;
	
      case 2: //1,2,0
	if(dims1[imo1[0]]*dims1[imo1[1]] > 0)
	  nb31 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d[0]*d[1]>0)
	  nb13 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	else nb13 = 1;
	if(nb13 < 1) nb13 = 1;

	for(i=0;i < numtasks;i++) {
	  p1 = (float *) recvbuf + *(RcvStrt+i);
	  p0 = (float *) dest + ds*(istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]])); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];

	  for(z=0;z<zsize;z+=nb31) { 
	    z2 = min(zsize,z+nb31);

	    for(x=0;x < xsize;x+=nb13) {
	      x2 = min(xsize,x+nb13);
	      
	      for(k=z;k<z2;k++) {
		pin1 = p1 + ds*(xsize*ysize*k+ x);
		pout1 = p0 + ds*(k +x*d[0]);
		for(y=0;y < ysize;y++) {
		  pin = pin1;
		  pout = pout1;
		  //		  printf("%d: y,z=%d %d, x from %d to %d\n",taskid,y,k,x,x2);
		  for(ii=x;ii<x2;ii++) {
		    for(j=0;j<ds;j++)
		      *pout++ = *pin++;
		    pout += ds*(d[0]-1);
		  }
		  pin1 += ds*xsize;
		  pout1+= ds*d[0]*d[1];
		}
	      }
	    }
	  }
	}
	
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	if(dims1[imo1[0]]*dims1[imo1[1]] > 0)
	  nb31 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	else nb31 = 1;
	if(nb31 < 1) nb31 = 1;
	if(d[0]*d[1]>0)
	  nb13 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	else nb13 = 1;
	if(nb13 < 1) nb13 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = (float *) recvbuf + *(RcvStrt+i);
	  p0 = (float *) dest + ds*(istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]])); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];

	  for(z=0;z<zsize;z+=nb31) { 
	    z2 = min(zsize,z+nb31);
	    for(x=0;x < xsize;x+=nb13) {
	      x2 = min(xsize,x+nb13);
	      
	      for(k=z;k<z2;k++) {
		pin1 = p1 + ds*(xsize*ysize*k + x);
		pout1 = p0 + ds*(k + x*d[0]*d[1]);
		for(y=0;y < ysize;y++) {
		  pin = pin1;
		  pout = pout1 + ds*y*d[0];
		  for(ii=x;ii<x2;ii++) {
		    for(j=0;j<ds;j++)
		      *pout++ = *pin++;
		    pout += ds*(d[0]*d[1]-1);
		  }
		  pin1 += ds*xsize;
		}
	      }
	    }
	  }
	}
	
	break;
      case 0: //2,0,1
	if(dims1[imo1[0]]*dims1[imo1[1]] > 0)
	  nb32 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d[0]*d[1]>0)
	  nb23 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = (float *) recvbuf + *(RcvStrt+i);
	  p0 = (float *) dest + ds*(istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]])); 
	  /*
	  xstart = istart[i][imo1[0]];
	  xend = iend[i][imo1[0]];
	  ystart = istart[i][imo1[1]];
	  yend = iend[i][imo1[1]];
	  zstart = istart[i][imo1[2]];
	  zend = iend[i][imo1[2]];
	  */
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];
	  for(z=0;z<zsize;z+=nb32) { 
	    z2 = min(zsize,z+nb32);
	    for(y=0;y < ysize;y+=nb23) {
	      y2 = min(ysize,y+nb23);
	      for(k=z;k<z2;k++) {
		pin = p1 + ds*xsize*(k*ysize + y);
		pout1 = p0 + ds*(k + y*d[0]*d[1]);
		for(j=y;j<y2;j++) {
		  // pin = pin1;
		  pout = pout1;
		  for(x=0;x < xsize;x++) {
		    for(j=0;j<ds;j++)
		      *pout++ = *pin++;
		    pout += ds*(d[0]-1);
		  }
		  pout1+=ds*d[0]*d[1];
		  //pin1 += xsize;
		}
	      }
	    }
	  }
	}
	
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {
	if(dims1[imo1[0]]*dims1[imo1[1]] > 0)
	  nb32 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	else nb32 = 1;
	if(nb32 < 1) nb32 = 1;
	if(d[0]*d[1]>0)
	  nb23 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	else nb23 = 1;
	if(nb23 < 1) nb23 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = (float *) recvbuf + *(RcvStrt+i);
	  p0 = (float *) dest + ds*(istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]])); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];
	  /*
	  xstart = istart[i][imo1[0]];
	  xend = iend[i][imo1[0]];
	  ystart = istart[i][imo1[1]];
	  yend = iend[i][imo1[1]];
	  zstart = istart[i][imo1[2]];
	  zend = iend[i][imo1[2]];
	  */
	  for(z=0;z<zsize;z+=nb32) { 
	    z2 = min(zsize,z+nb32);
	    for(y=0;y < ysize;y+=nb23) {
	      y2 = min(ysize,y+nb23);
	      for(k=z;k<z2;k++) {
		pin = p1 + ds*xsize*(k*ysize+ y);
		pout1 = p0 + ds*d[0]*(k + y*d[1]);
		for(j=y;j<y2;j++) {
		  //pin = pin1;
		  pout = pout1;
		  for(x=0;x < xsize*ds;x++) {
		    *pout++ = *pin++;
		  }
		  pout1 += ds*d[0]*d[1];
		  //pin1 += xsize;
		}
	      }
	    }
	  }
	}
      }        
    }
  }
}


// Perform transform followed by transpose (MPI exchange)
 template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::exec(char *in,char *out, bool OW) {
  //  Type1 *in;
  //Type2 *out;
  //in = (Type1 *) in_;
  //out = (Type2 *) out_;

   if(trplan->trans_type->is_empty) {
     mpiplan->exec(in,out);
     return;
   }

  int *tmpdims;

  tmpdims = trplan->grid2->ldims;

  Type2 *sendbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
#ifdef TIMERS
  double t1=MPI_Wtime();
#endif
  pack_sendbuf_trans(sendbuf,in,OW);
#ifdef TIMERS
      timers.packsend_trans += MPI_Wtime() -t1;
#endif
  tmpdims = mpiplan->grid2->ldims;

  Type2 *recvbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  MPI_Alltoallv(sendbuf,mpiplan->SndCnts,mpiplan->SndStrt,MPI_REAL,recvbuf,mpiplan->RcvCnts,mpiplan->RcvStrt,MPI_REAL,mpiplan->mpicomm);
#ifdef TIMERS
      timers.alltoall += MPI_Wtime() -t1;
#endif

  delete [] sendbuf;
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  mpiplan->unpack_recvbuf((Type2 *) out,recvbuf);
#ifdef TIMERS
      timers.unpackrecv += MPI_Wtime() -t1;
#endif
  delete [] recvbuf;

#ifdef DEBUG
  char str[80];
  sprintf(str,"transmpi.out%d.%d",cnt_trans++,mpiplan->taskid);
  int imo2[3];
  inv_mo(mpiplan->mo2,imo2);
  write_buf<Type2>((Type2 *) out,str,tmpdims,imo2);
#endif
}

// Perform transform followed by transpose (MPI exchange)
 template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::exec_deriv(char *in,char *out, bool OW) {
  //  Type1 *in;
  //Type2 *out;
  //in = (Type1 *) in_;
  //out = (Type2 *) out_;
  int *tmpdims;

  tmpdims = trplan->grid2->ldims;

  Type2 *sendbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
#ifdef TIMERS
  double t1=MPI_Wtime();
#endif
  pack_sendbuf_deriv(sendbuf,in,OW);
#ifdef TIMERS
      timers.packsend_deriv += MPI_Wtime() -t1;
#endif
  tmpdims = mpiplan->grid2->ldims;
  Type2 *recvbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  MPI_Alltoallv(sendbuf,mpiplan->SndCnts,mpiplan->SndStrt,MPI_REAL,recvbuf,mpiplan->RcvCnts,mpiplan->RcvStrt,MPI_REAL,mpiplan->mpicomm);
#ifdef TIMERS
      timers.alltoall += MPI_Wtime() -t1;
#endif

  delete [] sendbuf;
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  mpiplan->unpack_recvbuf((Type2 *) out,recvbuf);
#ifdef TIMERS
      timers.unpackrecv += MPI_Wtime() -t1;
#endif
  delete [] recvbuf;

#ifdef DEBUG
  char filename[80],str[80];
  sprintf(filename,"transmpi.out%d",cnt_trans++);
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  sprintf(str,".%d",taskid);
  strcat(filename,str);
  int imo2[3];
  inv_mo(mpiplan->mo2,imo2);
  write_buf<Type2>((Type2 *) out,filename,tmpdims,imo2);
#endif
}

template <class Type> void write_buf(Type *buf,char *filename,int sz[3],int mo[3]) {
  int i,j,k;
  FILE *fp;
  Type *p=buf;
  complex_double *p1;

  fp=fopen(filename,"w");
  fprintf(fp,"Size %d %d %d\n",sz[mo[0]],sz[mo[1]],sz[mo[2]]);
  for(k=0;k<sz[mo[2]];k++)
    for(j=0;j<sz[mo[1]];j++)
      for(i=0;i<sz[mo[0]];i++) {
	if(abs(*p) > 1.e-7) {
	  p1 = (complex_double *) p;
	  fprintf(fp,"(%d %d %d) %lg %lg\n",i,j,k,p1->real(),p1->imag());
	}
	p++;
      }
  fclose(fp); 
}


 template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::pack_sendbuf_trans(Type2 *sendbuf,char *src, bool OW)
{
  int i,nt,x,y,z,l;

  int d1 = mpiplan->d1;
  int d2 = mpiplan->d2;
  int *dims1 = mpiplan->dims1;
  int *mo1 = mpiplan->mo1;
  int *mo2 = mpiplan->mo2;
  int d[3],*tmpdims;
  int ds=sizeof(Type2)/4;
  float *p1,*p0,*pz,*p2,*buf;
  int imo1[3],imo2[3],xs,xe,ys,ye;

  nt = mpiplan->numtasks;
  int *SndStrt = mpiplan->SndStrt;
  int j,istart[nt][3],iend[nt][3];
  int comm_id = mpiplan-> comm_id;

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);

  for(i=0;i<nt;i++) 
    for(j=0;j<3;j++) {
      istart[i][j] = 0;
      iend[i][j] = dims1[j];
      //      istart[i][j] = mpiplan->grid1->st[comm_id][i][j]; 
      //iend[i][j] = mpiplan->grid1->en[comm_id][i][j];
  }
  for(j=0;j<nt;j++) {
    istart[j][d2] = mpiplan->grid2->st[comm_id][j][d2]; iend[j][d2] = mpiplan->grid2->en[comm_id][j][d2];
  }    

  tmpdims = trplan->grid2->ldims;
  //tmpdims[0]*tmpdims[1]*tmpdims[2]*
  if(OW && trplan->dt2 <= trplan->dt1)  // CC or C2R 
    buf = (float *) src;
  else
    buf = new float[ds*tmpdims[0]*tmpdims[1]*tmpdims[2]];

#ifdef DEBUG
  char str[80];
  sprintf(str,"pack_send_trans.in%d.%d",cnt_pack,mpiplan->taskid);
  //  inv_mo(mo2,imo2);
  // Optimize in the future
  write_buf<Type2>((Type2 *)src,str,dims1,imo1);
#endif

  /*
#ifdef FFTW
  if(trplan->dt1 > trplan->dt2 && !OW) {
    char *pin = new char[sizeof(Type1)*dims1[0]*dims1[1]*dims1[2]];
    trplan->exec(pin,(char *) buf,OW);
    delete [] pin;
  }
  else
#endif
  */
  trplan->exec(src,(char *) buf,OW);

#ifdef DEBUG
  sprintf(str,"pack_send_trans.out%d.%d",cnt_pack++,mpiplan->taskid);
  write_buf<Type2>((Type2 *) buf,str,tmpdims,imo2);
#endif
  
  for(i=0;i<3;i++)
    d[i] = dims1[imo1[i]];

for(i=0;i < nt;i++) {
    xs = istart[i][imo1[0]]*ds;
    xe = iend[i][imo1[0]]*ds;
    ys = istart[i][imo1[1]];
    ye = iend[i][imo1[1]];
    p1 = (float *) sendbuf + *(SndStrt+i);
    p0 = buf + ds * istart[i][imo1[0]];
    for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++) {
      pz = p0 + ds*d[0]*d[1]*z;
      for(y=ys;y< ye;y++) {
	p2 = pz +ds*d[0]*y;
	for(x=xs;x < xe;x++)
	  *p1++ = *p2++;;
      }
    }
}

 if((void *) buf != (void *) src) 
    delete [] buf;
}

 template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::pack_sendbuf_deriv(Type2 *sendbuf,char *src, bool OW)
{
  int i,nt,x,y,z,l;

  int d1 = mpiplan->d1;
  int d2 = mpiplan->d2;
  int *dims1 = mpiplan->dims1;
  int *mo1 = mpiplan->mo1;
  int *mo2 = mpiplan->mo2;
  int d[3],*tmpdims;
  int ds=sizeof(Type2)/4;
  float *p1,*p0,*pz,*p2,*buf;
  int imo1[3],imo2[3],xs,xe,ys,ye;

  nt = mpiplan->numtasks;
  int *SndStrt = mpiplan->SndStrt;
  int j,istart[nt][3],iend[nt][3];
  int comm_id = mpiplan-> comm_id;

  inv_mo(mo1,imo1);

  for(i=0;i<nt;i++) 
    for(j=0;j<3;j++) {
      istart[i][j] = 0;
      iend[i][j] = dims1[j];
      //      istart[i][j] = mpiplan->grid1->st[comm_id][i][j]; 
      //iend[i][j] = mpiplan->grid1->en[comm_id][i][j];
  }
  for(j=0;j<nt;j++) {
    istart[j][d2] = mpiplan->grid2->st[comm_id][j][d2]; iend[j][d2] = mpiplan->grid2->en[comm_id][j][d2];
  }    

  tmpdims = trplan->grid2->ldims;

  if(OW && trplan->dt2 <= trplan->dt1) 
    buf = (float *) src;
  else
    buf = new float[ds*tmpdims[0]*tmpdims[1]*tmpdims[2]];

#ifdef DEBUG
  char str[80];
  sprintf(str,"pack_send_trans.in%d",cnt_pack);
  //  inv_mo(mo2,imo2);
  // Optimize in the future
  write_buf<Type2>((Type2 *)src,str,dims1,imo1);
#endif

  trplan->exec_deriv(src,(char *) buf,OW);

#ifdef DEBUG
  sprintf(str,"pack_send_trans.out%d",cnt_pack++);
  write_buf<Type2>((Type2 *) buf,str,tmpdims,imo1);
#endif
  
  for(i=0;i<3;i++)
    d[i] = dims1[imo1[i]];

for(i=0;i < nt;i++) {
    xs = istart[i][imo1[0]]*ds;
    xe = iend[i][imo1[0]]*ds;
    ys = istart[i][imo1[1]];
    ye = iend[i][imo1[1]];
    p1 = (float *) sendbuf + *(SndStrt+i);
    p0 = buf + ds * istart[i][imo1[0]];
    for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++) {
      pz = p0 + ds*d[0]*d[1]*z;
      for(y=ys;y< ye;y++) {
	p2 = pz +ds*d[0]*y;
	for(x=xs;x < xe;x++)
	  *p1++ = *p2++;;
      }
    }
}

 if((void *) buf != (void *) src) 
    delete [] buf;
}

void inv_mo(int mo[3],int imo[3])
{
  int i,j;
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) {
      if(mo[i] == j)
	imo[j] = i;
    }
}

#ifdef MKL_BLAS

 template <class Type> void blas_trans(size_t rows,size_t cols,const double alpha,const Type *A,size_t lda,Type *B,size_t ldb)
 {
   if(typeid(Type) == type_float)      
     mkl_somatcopy('r','t',rows,cols,alpha,(const float *) A,lda,(float *) B,ldb);
   else if(typeid(Type) == type_double) {
     //     printf("Calling mkl_domatcopy, rows=%d, cols=%d\n",rows,cols);
     mkl_domatcopy('c','t',rows,cols,alpha,(const double *) A,lda,(double *) B,ldb);
   }
   else if(typeid(Type) == type_complex)
     mkl_comatcopy('r','t',rows,cols,alpha,(const mycomplex *) A,lda,(mycomplex *) B,ldb);
   else if(typeid(Type) == type_complex_double) {
     //printf("Calling mkl_zomatcopy, rows=%d, cols=%d\n",rows,cols);
     mkl_zomatcopy('r','t',rows,cols,alpha,(const complex_double *) A,lda,(complex_double *) B,ldb);
   }
 }
#endif


}
