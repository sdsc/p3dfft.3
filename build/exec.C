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
#include <unistd.h>

namespace p3dfft {

  extern int nslices;
  int nslices;
#ifdef CUDA
  extern cudaStream_t *streams;
  cudaStream_t *streams;
#endif

#include "reorder.C"

#ifdef TIMERS
  timer timers;
#endif

#ifdef DEBUG
static int cnt_pack=0;
static int cnt_trans=0;
static int cnt_mpi=0;
void printbuf(char *,int[3],int,int);
#endif

/* Execute 3D transform, optionally followed by spectral derivative in logical dimension idir.
   Input/output: buffers containing 3D arrays, contiguously stored according to memory mappings (mo1 and mo2) expected by transform3D class. Can be same buffer, otherwise need to be non-overlapping.
   OW: if != 0 input buffer can be overwritten (obviously the case if in==out).
*/

  template <class Type1,class Type2> void transform3D<Type1,Type2>::exec_deriv(Type1 *in,Type2 *out,int idir, bool OW, char *work_host, char *work_dev) 
{
  int ldir,i,j,k,g,mid;
  complex_double *p1,*p2,mult;

  char *buf[2],*buf_in,*buf_out;
  int next,curr,dt_1,dt_2,stage_cnt=0;

  if((void *) in == (void *) out && !OW) {
    printf("Warning in transform3D:exec_deriv: input and output are the same, but overwrite priviledge is not set\n");
    //    OW = 1;
  }


#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(grid1->Pgrid->mpi_comm_glob,&taskid);
  int *mo2;
#endif

  next = 1; curr = 0;
  dt_1 = dt1;
  buf[curr] = buf_in = (char *) in;
  buf_out = (char *) out;
  int cnt=0;
  int prev_t = 1;

  int *dims1 = grid1->Ldims;
  size_t size0 = grid1->Ldims[0]*grid1->Ldims[1]*grid1->Ldims[2]*dt1;
  size_t size_out = grid2->Ldims[0]*grid2->Ldims[1]*grid2->Ldims[2]*dt2;
  bool DevAlloc = false;
  bool HostAlloc = false;
  
#ifdef CUDA
  int prevLoc=LocHost;
  cudaStream_t *stream;
  int size1 = grid1->Ldims[0] * grid1->Ldims[1] * grid1->Ldims[2];
  event_t *event_hold=NULL;

  if(!work_dev) {
    checkCudaErrors(cudaMalloc(&work_dev,WorkSpaceDev));
    DevAlloc = true;
  }
  if(!work_host) {
    checkCudaErrors(cudaMallocHost(&work_host,WorkSpaceHost));
    HostAlloc=true;
  }
  //  checkCudaErrors(cudaMalloc((&(DevBuf[0])),size1*sizeof(Type1)));
  /*
  if(size_out >= maxDevSize && OutLoc == LocDevice)
    DevBuf[currdev] = (char *) out;
  else {
    checkCudaErrors(cudaMalloc(&(DevBuf[currdev]),maxDevSize));
    DevAlloc = true;
#ifdef DEBUG
    printf("Allocated DevBuf size %d\n",maxDevSize);
#endif
}*/

  /*
  if(InLoc == LocHost && Stages->kind != MPI_ONLY) {
  // if(!types[0]->isEmpty) ...
    for(i=0;i<nslices;i++) {
      stream = &(streams[i]);
      checkCudaErrors(cudaMemcpyAsync(DevBuf[0]+offset[i]*sizeof(Type1),in+offset[i],mysize[i]*sizeof(Type1),cudaMemcpyHostToDevice,*stream));
      checkCudaErrors(cudaEventRecord(Stages->EVENT_H2D,*stream));
      event_hold = &(Stages->EVENT_H2D);
    }
    buf[curr] = DevBuf[0];
    ndev = 1;
    currdev = 0;
    nextdev = 1;
    DevAlloc = true;
  }
  else { */

  buf[curr] = buf_in;
    
  prevLoc = InLoc;
  char *curr_work_dev = work_dev;
#else
  if(!work_host) {
    work_host = new char[WorkSpaceHost];
    HostAlloc = true;
  }
  int *event_hold=NULL;
  int prevLoc=0;
  char *curr_work_dev=NULL;
#endif

  stage *st;
  char *curr_work_host = work_host;
  char *tmp;

  for(stage *curr_stage=Stages;curr_stage != NULL;curr_stage = curr_stage->next) {

#ifdef DEBUG
    printf("stage %d, kind=%d\n",stage_cnt,curr_stage->kind);
#endif  
    stage_cnt++;
	//	transplan<Type1,Type2> *st = (transplan<Type1,Type2> *) curr_stage;
    st = curr_stage;
    int size1 = st->dims1[0] * st->dims1[1] * st->dims1[2]; 
    int size2 = st->dims2[0] * st->dims2[1] * st->dims2[2]; 
    dt_1 = curr_stage->dt1;
    dt_2 = curr_stage->dt2;

#ifdef CUDA
    if(curr_stage->kind == TRANS_ONLY && prevLoc != ((transplan<Type1,Type2> *) curr_stage)->InLoc)
      printf("Error in transform3D: discrepancy between input and output buffer locations, at stage %d, %d %d\n",stage_cnt,prevLoc,((transplan<Type1,Type2> *) curr_stage)->InLoc);
#endif

    if(!curr_stage->next) /* && OutLoc == LocDevice && (void *) out != (void *) buf[curr]) {
      // Last stage
      buf[next] = buf_out;
      //if(curr_stage->OutLoc == LocHost)
      //  printf("Warning: last stage output location differs from end location\n");
    }
    else if(!curr_stage->next && OutLoc == LocHost) {
    // Last stage*/
      buf[next] = buf_out;
      //      if(curr_stage->OutLoc == LocDevice)
      //  printf("Warning: last stage output location differs from end location\n");
    
#ifdef CUDA
    else if(curr_stage->kind == TRANS_ONLY && ((transplan<Type1,Type2> *) curr_stage)->OutLoc == LocDevice) {
      buf[next] = work_dev;
      curr_work_dev = max(curr_work_dev,work_dev + size2 * dt_2 *st->stage_prec); 
    }
    else
     // Not last stage; output = host
      //      if (curr_stage->InLoc != LocHost)
      if(OW && size2*dt_2 <= size0 && InLoc == LocHost) 
	// Can use input buffer and avoid allocating a temp. buffer
	buf[next] = buf_in;
      else if(size2*dt_2 <= size_out && OutLoc == LocHost)
	// Can use output buffer and avoid allocating a temp. buffer
	buf[next] = buf_out;
      else { 
	buf[next] = work_host;
	curr_work_host = max(curr_work_host,work_host+size2*dt_2*st->stage_prec);
      }
  	  // Need to allocate a temp. buffer on host
	  /*
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
	  checkCudaErrors(cudaMallocHost(&(var[nextvar]),size2*dt_2*st->stage_prec));
	  buf[next] = var[nextvar];  //[nextvar]
	  nvar++;
	  nextvar = 1-nextvar;
	  }*/

	/*
      else // inloc == lochost
	if(size2*dt_2 <= size1*dt_1 && (OW || (void *) buf[curr] != (void *) buf_in))
	  // Reuse the current buffer
	  buf[next] = buf[curr];
	else {
	  // Need to allocate a temp. buffer
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
	  checkCudaErrors(cudaMallocHost(&(var[nextvar]),size2*dt_2*st->stage_prec));
	  buf[next] = var[nextvar];  //[nextvar]
	  nvar++;
	  nextvar = 1-nextvar;
	}
	*/
    
#else //non-CUDA
    
    else 
      if(!OW && buf[curr] == buf_in) { //|| size2*dt_2 > size1 *dt_1) {
	//  out-place
	// check if we can use the destination array as the output buffer
	if(dt_2 * size2 <= dt2 * grid2->Ldims[0] *grid2->Ldims[1] *grid2->Ldims[2]) {
	  buf[next] = buf_out;
#ifdef DEBUG
	  printf("%d: Using output buffer in out-of-place transform\n",taskid);
#endif
	}
	
	else if(size2 * dt_2 > size0) { //  allocate a new buffer
	  buf[next] = work_host;
	  curr_work_host += size2*dt_2*st->stage_prec;
	  /*
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
	  var[nextvar] = new char[size2*dt_2*st->stage_prec];
	  buf[next] = var[nextvar];
	  nvar++;
	  nextvar = 1-nextvar;
	  */
	}
		  
      }
      else { // In-place
	buf[next] = work_host;
	curr_work_host = max(curr_work_host,work_host+size2*dt_2*st->stage_prec);
      }

#endif
      
      if(curr_stage->kind == TRANS_ONLY) {
      // Only transform, no exchange
      
	transplan<Type1,Type2> *tr = (transplan<Type1,Type2> *) st;

#ifdef DEBUG
	  mo2 = tr->mo2;
#endif

#ifdef CUDA    
	  //	  if(tr->OutLoc != LocHost && buf[next] != buf_out)
	  //  buf[next] = DevBuf[currdev];

	  //	  if(!tr->trans_type->is_empty) {

	    if(tr->InLoc == LocHost) {
	      char *devbuf = work_dev;
	      curr_work_dev += size1*dt_1*st->stage_prec;

	//        if(currdev < 0) {
	// currdev = 1 - nextdev;
          //      nextdev = 1;
	  
	//	checkCudaErrors(cudaMalloc(&(DevBuf[currdev],st->stage_prec*max(size1*dt_1,size2*dt_2)));
	      for(i=0;i<nslices;i++) {
		stream = &(streams[i]);
		checkCudaErrors(cudaMemcpyAsync(devbuf+tr->offset1[i]*st->stage_prec*dt_1,buf[curr]+tr->offset1[i] * st->stage_prec*dt_1,tr->mysize1[i] * st->stage_prec*dt_1,cudaMemcpyHostToDevice,*stream));
	  //checkCudaErrors(cudaMemcpy(DevBuf[currdev]+st->offset1[i]*st->stage_prec*dt_1,buf[curr]+st->offset1[i] * st->stage_prec*dt_1,st->mysize1[i] * st->stage_prec*dt_1,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaEventRecord(tr->EVENT_H2D,*stream));
	      }
	      event_hold = &(tr->EVENT_H2D);
	      buf[curr] = devbuf;
	    }

	    tmp = curr_work_dev;
	    if(tr->OutLoc == LocDevice)
	      curr_work_dev = work_dev + size2*dt_2*st->stage_prec;
	    else
	      curr_work_dev = work_dev;
      
        /*      
      nextdev = 1 - currdev;
      checkCudaErrors(cudaMalloc(&(DevBuf[nextdev]),size2 * st->stage_prec*dt_2));
      DevAlloc2 = true;
      if(st->OutLoc != LocHost && buf[next] != buf_out)
        buf[next] = DevBuf[nextdev];
      */
	    /*
	    if(tr->InLoc == LocHost)
	      tr->DevBuf = DevBuf[currdev];
	    if(!arcmp(tr->mo1,tr->mo2,3) && size2*dt_2 > size1*dt_1) {
	      dev_out = devbuf + size1*dt_1*st->stage_prec;

		checkCudaErrors(cudaMalloc(&(DevBuf[nextdev]),size2*dt_2*st->stage_prec));
	      tr->DevBuf2 = DevBuf[nextdev];
	      DevAlloc2 = true;

	      }
	      else
		dev_out = devbuf; //tr->DevBuf2 = DevBuf[currdev];
		}*/
	    //}

	  /*
	  if(tr->grid1->MemOrder[0] == tr->grid2->MemOrder[1] && tr->grid1->MemOrder[1] == tr->grid2->MemOrder[0]) { 
	    checkCudaErrors(cudaMalloc(&(tr->DevBuf3),size2*dt_2*st->stage_prec));
	    DevAlloc3 = true;
	  }
	  */
#else
	    tmp = curr_work_host;
#endif


	  if(dt_1 == dt_2)
	    if(prev_t == 1) {
	      transplan<Type1,Type1> *tr = (transplan<Type1,Type1> *) st;
	      if(!tr->is_empty)
		for(i=0;i<nslices;i++)
		  tr->exec_slice(buf[curr],buf[next],idir,i,nslices,event_hold,OW || buf[curr] != (char *) in,tmp);
	      else
		tr->exec_slice(buf[curr],buf[next],idir,0,1,event_hold,OW || buf[curr] != (char *) in,tmp);

	}
	else {
	  transplan<Type2,Type2> *tr = (transplan<Type2,Type2> *) st;
	  if(!tr->is_empty)
	    for(i=0;i<nslices;i++)
	      tr->exec_slice(buf[curr],buf[next],idir,i,nslices,event_hold,OW || buf[curr] != (char *) in,tmp);
	  else
	    tr->exec_slice(buf[curr],buf[next],idir,0,1,event_hold,OW || buf[curr] != (char *) in,tmp);
	}
	  else {
	    transplan<Type1,Type2> *tr = (transplan<Type1,Type2> *) st;
	    if(!tr->is_empty)
	      for(i=0;i<nslices;i++)
		tr->exec_slice(buf[curr],buf[next],idir,i,nslices,event_hold,OW || buf[curr] != (char *) in,tmp);
	    else
	      tr->exec_slice(buf[curr],buf[next],idir,0,1,event_hold,OW || buf[curr] != (char *) in,tmp);
	    
	    prev_t = 2;
	  }
      
	  /*	  if(DevAlloc2) {
	    checkCudaErrors(cudaFree(DevBuf[nextdev]));
	//	currdev = 1-currdev;
	//nextdev = 1-nextdev;
	    DevAlloc2 = false;
	  }
	  */
#ifdef CUDA
	  if(tr->OutLoc == LocHost) { // && Dev Alloc2
	    cudaEventSynchronize(tr->EVENT_D2H);
	//      cudaFree(devbuf);
	    event_hold = NULL;
	  }
	  else
	    event_hold = &(tr->EVENT_EXEC);
	  prevLoc = tr->OutLoc;
#endif

      }    
      else if(curr_stage->kind == MPI_ONLY) { // Only MPI plan (exchange, no transform)
	if(prev_t == 1) {
	  MPIplan<Type1> *tr = (MPIplan<Type1> *) curr_stage; 
#ifdef DEBUG
	  mo2 = tr->mo2;
#endif
	  tr->exec(buf[curr],buf[next],curr_work_host);
	}
	else {
	  MPIplan<Type2> *tr = (MPIplan<Type2> *) curr_stage; 
#ifdef DEBUG
	  mo2 = tr->mo2;
#endif
	  tr->exec(buf[curr],buf[next],curr_work_host);
	}
	prevLoc = LocHost;
      }
      else { // MPI and transform combined
	
	if(dt_1 == dt_2)
	  if(prev_t == 1) {
	    trans_MPIplan<Type1,Type1> *tr = (trans_MPIplan<Type1,Type1> *) (transplan<Type1,Type1> *) curr_stage;
#ifdef DEBUG
	    mo2 = tr->trplan->mo2;
#endif
	    tr->exec(buf[curr],buf[next],idir,event_hold,OW  || buf[curr] != (char * ) in,curr_work_host,curr_work_dev);
	  }
	  else {
	    trans_MPIplan<Type2,Type2> *tr = (trans_MPIplan<Type2,Type2> *) (transplan<Type2,Type2> *) curr_stage;
#ifdef DEBUG
	    mo2 = tr->trplan->mo2;
#endif
	    tr->exec(buf[curr],buf[next],idir,event_hold,OW  || buf[curr] != (char *) in,curr_work_host,curr_work_dev);
	  }
	else {
	  trans_MPIplan<Type1,Type2> *tr = (trans_MPIplan<Type1,Type2> *) (transplan<Type1,Type2> *) curr_stage;
#ifdef DEBUG
	  mo2 = tr->trplan->mo2;
#endif
	  tr->exec(buf[curr],buf[next],idir,event_hold,OW  || buf[curr] != (char *) in,curr_work_host,curr_work_dev);
	  prev_t = 2;
	}
#ifdef CUDA
	prevLoc = LocHost;
	curr_work_dev = work_dev;
#endif
      }
      
#ifdef DEBUG
      int imo[3];
      char str[80];
      for(i=0;i<3;i++)
	imo[mo2[i]] = i; 
      sprintf(str,"exec-out.%d.%d",stage_cnt,taskid);
      write_buf<Type2>((Type2 *) buf[next],str,curr_stage->dims2,imo);
#endif
      
      /*
      if(nvar > 1) { // Keep track and delete buffers that are no longer used
#ifdef CUDA
	checkCudaErrors(cudaFreeHost(var[nextvar]));
#else
	delete [] var[nextvar];
#endif
	nvar--;
      }
      */

      next = 1-next;
      curr = 1-curr;

      /*     if(currdev >= 0)
	     checkCudaErrors(cudaFree(DevBuf[currdev]));
    if(DevAlloc2) {
      currdev = nextdev;
      nextdev = 1 - currdev;
      DevAlloc2 = false;
    }
    else
      currdev = -1;
  */
  }
  
#ifdef CUDA
  if(prevLoc != OutLoc) {
    size_t size = grid2->Ldims[0]*grid2->Ldims[1]*grid2->Ldims[2]*sizeof(Type2);
    checkCudaErrors(cudaEventSynchronize(*event_hold));
    if(OutLoc == LocHost) {
      checkCudaErrors(cudaMemcpy(out,buf[curr], size, cudaMemcpyDeviceToHost));
    }
    else
      checkCudaErrors(cudaMemcpy(out,buf[curr],  size, cudaMemcpyHostToDevice));

    //    checkCudaErrors(cudaFree(buf[curr]));
  }

#endif

  /*
  if(nvar > 0)
#ifdef CUDA
    checkCudaErrors(cudaFreeHost(var[1-nextvar]));
#else
  delete [] var[1-nextvar];
#endif
  */
#ifdef CUDA
  if(DevAlloc)
    checkCudaErrors(cudaFree(work_dev));
  if(HostAlloc)
    checkCudaErrors(cudaFreeHost(work_host));
#else
  if(HostAlloc)
    delete [] work_host;
#endif

    
  
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
  template <class Type1,class Type2> void transform3D<Type1,Type2>::exec(Type1 *in,Type2 *out, bool OW, char *work_host, char *work_dev)
//void transform3D::exec(char *in,char *out,int OW)
{
  // Call exec_deriv with -1 (void) for derivative dimension
  exec_deriv(in,out,-1, OW, work_host, work_dev);
}

  // Execute 1D transform, combined with local transpose as needed
// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]

// Execute transform (combined with local transpose if needed), on a slice of data asynchronously via streams
// Expect the calling function to set up and initilize device buffer (asynchronously)
  template <class Type1,class Type2> void transplan<Type1,Type2>::exec(char *in_,char *out_, int dim_deriv,bool OW,char *tmpbuf)
{
  int i;
  bool tmpalloc = false;
  

#ifdef CUDA
  size_t size1 = dims1[0]*dims1[1]*dims1[2] * sizeof(Type1);
  size_t size2 = dims2[0]*dims2[1]*dims2[2] * sizeof(Type2);
  cudaStream_t *stream;
  cudaEvent_t *event_hold=NULL;

  if(!tmpbuf) {
    checkCudaErrors(cudaMalloc(&tmpbuf,WorkSpace));
    tmpalloc=true;
  }
      /*
  if(InLoc == LocDevice)
    DevBuf = in_;
  else {
    checkCudaErrors(cudaMalloc((&(DevBuf)),size1));
    DevAlloc = true;
  }
  if(OutLoc == LocDevice)
    DevBuf2 = out_;
  else if(size2 <= size1)
    DevBuf2 = DevBuf;
  else {
    checkCudaErrors(cudaMalloc((&(DevBuf2)),size2));
    DevAlloc2 = true;
  }
*/
  for(i=0;i<nslices;i++) {
    stream = &(streams[i]);
    if(InLoc != LocDevice) {
      checkCudaErrors(cudaMemcpyAsync(tmpbuf+offset1[i]*sizeof(Type1),in_ + offset1[i]*sizeof(Type1),mysize1[i]*sizeof(Type1),cudaMemcpyHostToDevice,*stream));
      //checkCudaErrors(cudaMemcpy(DevBuf+offset1[i]*sizeof(Type1),in_ + offset1[i]*sizeof(Type1),mysize1[i]*sizeof(Type1),cudaMemcpyHostToDevice));
      checkCudaErrors(cudaEventRecord(EVENT_H2D,*stream));
      event_hold = &EVENT_H2D;
      exec_slice(tmpbuf,out_,dim_deriv,i,nslices,event_hold,OW,tmpbuf + size1);
    }
    else
      exec_slice(in_,out_,dim_deriv,i,nslices,event_hold,OW,tmpbuf);
  }

  cudaDeviceSynchronize();
  if(tmpalloc)
    checkCudaErrors(cudaFree(tmpbuf));
#else
  if(!tmpbuf) {
    tmpbuf = new char[WorkSpace];
    tmpalloc=true;
  }
  for(i=0;i<nslices;i++) 
    exec_slice(in_,out_,dim_deriv,i,nslices,NULL,OW,tmpbuf);  
  if(tmpalloc)
    delete [] tmpbuf;
#endif

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
//#ifdef CUDA
  template <class Type1,class Type2> int transplan<Type1,Type2>::reorder_trans_slice(Type1 *in,Type2 *out,int *mo1,int *mo2,int dim_deriv,int slice, int nslices,event_t *event_hold, bool OW, char *tmpbuf)
   //#else
 // template <class Type1,class Type2> int transplan<Type1,Type2>::reorder_trans_slice(Type1 *in,Type2 *out,int *mo1,int *mo2,int dim_deriv,int slice, int nslices,event_t *event_hold, bool OW)
 //#endif
{
int mc[3],i,j,k,ii,jj,kk,i2,j2,k2,cmpl;
  void rel_change(int *,int *,int *);

  int imo1[3],imo2[3],d1[3],d2[3];
  int scheme,nb13,nb31,nb23,nb32;
  
#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(grid1->Pgrid->mpi_comm_glob,&taskid);
  printf("%d: In transplan::reorder_trans, mo1=%d %d %d, mo2=%d %d %d, dims2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2],dims2[0],dims2[1],dims2[2]);
#endif
 
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
  bool deriv;
  if(dim_deriv == trans_dim)
    deriv = true;
  else
    deriv = false;

  if(mo1[trans_dim] == 0)
    scheme = TRANS_IN;
  else if(mo2[trans_dim] == 0)
    scheme = TRANS_OUT;
  else if(is_empty) 
    scheme = TRANS_OUT;
  else {
    printf("Error in reorder_trans: expected dimension %d to be the leading dimension for input or output\n",trans_dim);
    return (-1);
  }

  // Find inverse memory mappings
  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  // Find local storage dimensions
  for(i=0;i<3;i++) {
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }

  rel_change(imo1,imo2,mc);

#ifdef DEBUG
  printf("In reorder_trans, mc=%d %d %d, scheme=%s\n",mc[0],mc[1],mc[2],(scheme == TRANS_IN) ? "IN" : "OUT"); 
  char str[80];
  static int cnt_reorder_trans=0;
#endif

  bool inplace=false;

  if(scheme == TRANS_IN) { // Scheme is transform before reordering
  
  switch(mc[0]) {
 case 1:
     switch(mc[1]) {
      case 0: //1,0,2
  
        if(OW && dt2 <= dt1)
	  inplace=true;
        rot102in_slice(in,out,inplace,d1,d2,trans_type->exec,plan,slice,nslices,deriv,tmpbuf);
        cmpl = SLICE;
  
        break;
      	
      case 2: //1,2,0

        if(slice < nslices-1)
          return(NONE);
        cmpl = FULL;
#ifdef CUDA
        if(event_hold != NULL)
          checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
        rot120in(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);
	
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
        if(slice < nslices-1)
          return(NONE);
        cmpl = FULL;
#ifdef CUDA
        if(event_hold != NULL)
          checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
        rot210in(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);

        break;

	break;

      
      case 0: //2,0,1
        if(slice < nslices-1)
          return(NONE);
        cmpl = FULL;
#ifdef CUDA
        if(event_hold != NULL)
          checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
        rot201in(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);

	break;
      }
      break;
 case 0: //0,2,1
    if(mc[1] == 2) {
      cmpl = FULL;
      if((void *) in == (void *) out) {
        if(slice < nslices-1)
          return(NONE);
#ifdef CUDA
        if(event_hold != NULL)
          checkCudaErrors(cudaEventSynchronize(*event_hold));	
#endif
        rot021_ip(in,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);
      }
      else 
        rot021_op_slice(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,slice,nslices,deriv);
      break;
    }
    
  }
  }
  else { // Scheme is transform after reordering

    if(cmpmo(mc,102)) {
  
      rot102out_slice(in,out,d1,d2,trans_type->exec,plan,dt1,dt2,slice,nslices,deriv,tmpbuf);
      cmpl = SLICE;
      /*
#ifdef DEBUG
      char str[80];
      sprintf(str,"reorder102.out.%d",grid1->Pgrid->taskid);
#ifdef CUDA
      size_t size1 = d2[0]*d2[1]*d2[2];
      Type2 *tmp=new Type2[size1];
      checkCudaErrors(cudaMemcpy(tmp,out,size1*sizeof(Type2),cudaMemcpyDeviceToHost));
      write_buf<Type2>(tmp,str,dims2,imo2);
      delete [] tmp;
#else
      write_buf<Type2>(out,str,dims2,imo2);
#endif
#endif
      */
    }

    else if(cmpmo(mc,120)) {

      if(slice < nslices-1)
	return(NONE);
      cmpl = FULL;
#ifdef CUDA
      if(event_hold != NULL)
	checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
      rot120out(in,out,inplace,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);
    }
    else if(cmpmo(mc,210)) {
      if(slice < nslices-1)
	return(NONE);
      cmpl = FULL;
#ifdef CUDA
      if(event_hold != NULL)
	checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
      rot210out(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);
    }
    else if(cmpmo(mc,201)) {
      if(slice < nslices-1)
	return(NONE);
      cmpl = FULL;
#ifdef CUDA
      if(event_hold != NULL)
	checkCudaErrors(cudaEventSynchronize(*event_hold));
#endif
      rot201out(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,deriv,tmpbuf);
    }
    else if(cmpmo(mc,021)) {
      cmpl = FULL;
      if((void *) in == (void *) out) {
        if(slice < nslices-1)
          return(NONE);
#ifdef CUDA
        if(event_hold != NULL)
          checkCudaErrors(cudaEventSynchronize(*event_hold));	
#endif
        rot021_ip(in,d1,d2,trans_type->exec,plan,CACHE_BL,deriv);
      }
      else 
        rot021_op_slice(in,out,d1,d2,trans_type->exec,plan,CACHE_BL,slice,nslices,deriv);
    }
  }

#ifdef DEBUG
  sprintf(str,"reorder_trans.out%d.$d",cnt_reorder_trans++,Pgrid->taskid);
#ifdef CUDA
  Type2 *tmp=new Type2[dims2[0]*dims2[1]*dims2[2]];
  cudaMemcpy(tmp,out,dims2[0]*dims2[1]*dims2[2]*sizeof(Type2),cudaMemcpyDeviceToHost);
  write_buf<Type2>(tmp,str,dims2,imo2);
#else
  write_buf<Type2>(out,str,dims2,imo2);
#endif
#endif

  return(cmpl);
}
         
 bool cmpmo(int mo[3],int rhs)
 {
   if(mo[0] == rhs/100 && mo[1] == rhs/10-mo[0]*10 && mo[2] == rhs-mo[1]*10-mo[0]*100)
     return true;
   else
     return false;
 }

 /*
 template<class Type>  void ro102(Type *in,Type *out,int d1,int d2) {

   int i,j;
#ifdef MKL_BLAS
   blas_trans<Type2>(d1,d2,1.0,in,d1,out,d2);
#else

   for(j=0;j < d1;j++) {
     Type *pin = in + j;
     for(i=0;i < d2;i++) {
       *out++ = *pin;
       pin += d1;
     }

   }
#endif
 }
 */ 




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
  template <class Type> void MPIplan<Type>::exec(char *in_,char *out_, char *tmpbuf) {
  Type *in,*out;
  in = (Type *) in_;
  out = (Type *) out_;
  //  Type *sendbuf = new Type[dims1[0]*dims1[1]*dims1[2]*4];
  Type *sendbuf = (Type *) tmpbuf;
  tmpbuf += dims1[0]*dims1[1]*dims1[2]*sizeof(Type);
  //  void MPI_Alltoallv4(char *sendbuf,int *sndcnts,int *sndstrt,char *recvbuf,int *rcvcnts,int *rcvstar,MPI_Comm mpicomm,int n);

#ifdef TIMERS
  double t1=MPI_Wtime();
#endif
  pack_sendbuf(sendbuf,(Type *) in);
#ifdef TIMERS
  timers.packsend += MPI_Wtime() -t1;
#endif
      //  Type *recvbuf = new Type[dims2[0]*dims2[1]*dims2[2]*4];
  Type *recvbuf = (Type *) tmpbuf;
  tmpbuf += dims2[0]*dims2[1]*dims2[2]*sizeof(Type);
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  MPI_Alltoallv(sendbuf,SndCnts,SndStrt,MPI_BYTE,recvbuf,RcvCnts,RcvStrt,MPI_BYTE,Pgrid->mpicomm[comm_id]);
  //  MPI_Alltoallv4((float *)sendbuf,SndCnts,SndStrt,(float *) recvbuf,RcvCnts,RcvStrt,grid1.mpicomm[mpicomm_ind],numtasks);
#ifdef TIMERS
  timers.alltoall += MPI_Wtime() -t1;
#endif
  //  delete [] sendbuf;
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  unpack_recvbuf((Type *) out,recvbuf);
#ifdef TIMERS
      timers.unpackrecv += MPI_Wtime() -t1;
#endif

      //  delete [] recvbuf;

#ifdef DEBUG
  char str[80];
  sprintf(str,"mpiplan.out%d.%d",cnt_mpi++,Pgrid->taskid);
  int imo2[3];
  inv_mo(mo2,imo2);
  write_buf<Type>((Type *) out,str,dims2,imo2);
#endif

}

  void MPI_Alltoallv4(char *sendbuf,int *sndcnts,int *sndstrt,char *recvbuf,int *rcvcnts,int *rcvstrt,MPI_Comm mpicomm,int n)
  {
    int i;
    MPI_Status status[n];
    MPI_Request req[n];

    for(i=0;i<n;i++)
      MPI_Irecv(recvbuf+rcvstrt[i],rcvcnts[i],MPI_BYTE,i,MPI_ANY_TAG,mpicomm,req+i);
    for(i=0;i<n;i++)
      MPI_Send(sendbuf+sndstrt[i],sndcnts[i],MPI_BYTE,i,0,mpicomm);

    MPI_Waitall(n,req,status);
  }

template <class Type> void MPIplan<Type>::pack_sendbuf(Type *sendbuf,Type *src)
{
  int i,x,y,z;
  Type *p1,*p0;

  int j,istart[numtasks][3],iend[numtasks][3],d[3];
  int imo1[3],imo2[3];
  //  int ds=sizeof(Type)/4;

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);

  for(i=0;i<numtasks;i++)
    for(j=0;j<3;j++) {
      istart[i][j] = 0;
      iend[i][j] = dims1[j];
      //      istart[i][j] = grid1->st[comm_id][i][j]; iend[i][j] = grid1->en[comm_id][i][j];
    }
  for(j=0;j<numtasks;j++) {
    istart[j][d2] = grid2->st[d2][j]; iend[j][d2] = grid2->en[d2][j];
  }    


  for(i=0;i<3;i++)
    d[i] = dims1[imo1[i]];
  
  for(i=0;i < numtasks;i++) {
    p1 = sendbuf + *(SndStrt+i)/sizeof(Type);
    for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++)
      for(y=istart[i][imo1[1]];y< iend[i][imo1[1]];y++) {
	p0 = src+d[0]*(z*d[1]+y) +istart[i][imo1[0]];
	for(x=istart[i][imo1[0]];x < iend[i][imo1[0]];x++)
	  *p1++ = *p0++;;
      }
  }

}

template <class Type> void MPIplan<Type>::unpack_recvbuf(Type *dest,Type *recvbuf)
{
  int i,ii,x,y,z,k,x2,y2,z2;
  //  int ds=sizeof(Type)/4;
  Type *p1,*p0,*pin,*pin1,*pout,*pout1;

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
    isize[j][d1] = grid1->sz[d1][j];
    istart[j][d1] = grid1->st[d1][j]; iend[j][d1] = grid1->en[d1][j];
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
      p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
      for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++)
	for(y=istart[i][imo1[1]];y < iend[i][imo1[1]];y++) {
	  p0 = dest + d[0]*(z*d[1]+y) + istart[i][imo1[0]];
	  for(x=istart[i][imo1[0]];x < iend[i][imo1[0]];x++)
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
	  p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
	  p0 = dest + istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]]); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];
	  for(z=0;z < zsize;z++) 
	    for(y=0;y < ysize;y++) {
	      pout = p0 + d[0]*d[1]*z +y;
	      for(x=0;x < xsize;x++) {
		*pout = *p1++;
		pout += d[0];
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
	  p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
	  p0 = dest + istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]]); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];

	  for(z=0;z<zsize;z+=nb31) { 
	    z2 = min(zsize,z+nb31);

	    for(x=0;x < xsize;x+=nb13) {
	      x2 = min(xsize,x+nb13);
	      
	      for(k=z;k<z2;k++) {
		pin1 = p1 + xsize*ysize*k+ x;
		pout1 = p0 + d[0]*(k +x*d[1]);
		for(y=0;y < ysize;y++) {
		  pin = pin1;
		  pout = pout1;
		  //		  printf("%d: y,z=%d %d, x from %d to %d\n",taskid,y,k,x,x2);
		  for(ii=x;ii<x2;ii++) {
		    *pout++ = *pin++;
		    pout += d[0]*d[1];
		  }
		  pin1 += xsize;
		  pout1++;
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
	  p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
	  p0 = dest + istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]]); 
	  xsize = isize[i][imo1[0]];
	  ysize = isize[i][imo1[1]];
	  zsize = isize[i][imo1[2]];

	  for(z=0;z<zsize;z+=nb31) { 
	    z2 = min(zsize,z+nb31);
	    for(x=0;x < xsize;x+=nb13) {
	      x2 = min(xsize,x+nb13);
	      
	      for(k=z;k<z2;k++) {
		pin1 = p1 + xsize*ysize*k + x;
		pout1 = p0 + k + x*d[0]*d[1];
		for(y=0;y < ysize;y++) {
		  pin = pin1;
		  pout = pout1 + y*d[0];
		  for(ii=x;ii<x2;ii++) {
		    *pout++ = *pin++;
		    pout += d[0]*d[1];
		  }
		  pin1 += xsize;
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
	  p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
	  p0 = dest + istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]]); 
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
		pin = p1 + xsize*(k*ysize + y);
		pout1 = p0 + k + y*d[0]*d[1];
		for(j=y;j<y2;j++) {
		  // pin = pin1;
		  pout = pout1;
		  for(x=0;x < xsize;x++) {
		    *pout++ = *pin++;
		    pout += d[0];
		  }
		  pout1+=d[0]*d[1];
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
	  p1 = recvbuf + *(RcvStrt+i)/sizeof(Type);
	  p0 = dest + istart[i][imo2[0]] + d[0]*(istart[i][imo2[1]] +d[1]*istart[i][imo2[2]]); 
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
		pin = p1 + xsize*(k*ysize+ y);
		pout1 = p0 + d[0]*(k + y*d[1]);
		for(j=y;j<y2;j++) {
		  //pin = pin1;
		  pout = pout1;
		  for(x=0;x < xsize;x++) {
		    *pout++ = *pin++;
		  }
		  pout1 += d[0]*d[1];
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
    template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::exec(char *in,char *out, int dim_deriv,event_t *event_hold,bool OW,char *tmpbuf,char *devbuf) {
  //  Type1 *in;
  //Type2 *out;
  //in = (Type1 *) in_;
  //out = (Type2 *) out_;

   if(trplan->trans_type->is_empty) {
     mpiplan->exec(in,out);
     return;
   }
   
  int *tmpdims;

  tmpdims = trplan->grid2->Ldims;

  //  Type2 *sendbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
  Type2 *sendbuf = (Type2 *) tmpbuf;
  tmpbuf += tmpdims[0]*tmpdims[1]*tmpdims[2]*sizeof(Type2);

#ifdef TIMERS
  double t1=MPI_Wtime();
#endif
  pack_sendbuf_trans(sendbuf,in,dim_deriv,event_hold,devbuf,OW,tmpbuf);
#ifdef TIMERS
      timers.packsend_trans += MPI_Wtime() -t1;
#endif
  tmpdims = mpiplan->grid2->Ldims;

  //  Type2 *recvbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
  Type2 *recvbuf = (Type2 *) tmpbuf;
  tmpbuf += tmpdims[0]*tmpdims[1]*tmpdims[2]*sizeof(Type2);

#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  //  printf("%d: Calling mpi_alltoallv; tmpdims= %d %d %d, SndCnts=%d %d, RcvCnts=%d %d\n",mpiplan->taskid,tmpdims[0],tmpdims[1],tmpdims[2],mpiplan->SndCnts[0],mpiplan->SndCnts[1],mpiplan->RcvCnts[0],mpiplan->RcvCnts[1]);
  MPI_Alltoallv(sendbuf,mpiplan->SndCnts,mpiplan->SndStrt,MPI_BYTE,recvbuf,mpiplan->RcvCnts,mpiplan->RcvStrt,MPI_BYTE,mpiplan->Pgrid->mpicomm[mpiplan->comm_id]);
#ifdef TIMERS
      timers.alltoall += MPI_Wtime() -t1;
#endif

      //  delete [] sendbuf;
#ifdef TIMERS
  t1=MPI_Wtime();
#endif
  mpiplan->unpack_recvbuf((Type2 *) out,recvbuf);
#ifdef TIMERS
      timers.unpackrecv += MPI_Wtime() -t1;
#endif
      //  delete [] recvbuf;

#ifdef DEBUG
  char str[80];
  sprintf(str,"transmpi.out%d.%d",cnt_trans++,mpiplan->Pgrid->taskid);
  int imo2[3];
  inv_mo(mpiplan->mo2,imo2);
  write_buf<Type2>((Type2 *) out,str,tmpdims,imo2);
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


  template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::pack_sendbuf_trans(Type2 *sendbuf,char *src,  int dim_deriv,event_t *event_hold,char *devbuf,bool OW,char *tmpbuf)
{
  int i,nt,x,y,z,l;

  int d1 = mpiplan->d1;
  int d2 = mpiplan->d2;
  int *dims1 = mpiplan->dims1;
  int *mo1 = trplan->mo1;
  int *mo2 = trplan->mo2;
  int d[3],*tmpdims;
  //int ds=sizeof(Type2)/4;
  Type2 *p1,*p0,*pz,*p2,*buf;
  int imo1[3],imo2[3],xs,xe,ys,ye;
  bool buf_alloc=false;
#ifdef CUDA
  bool DevAlloc = false;
  bool DevAlloc2 = false;
#endif

  nt = mpiplan->numtasks;
  int *SndStrt = mpiplan->SndStrt;
  int j,istart[nt][3],iend[nt][3];
  //  int comm_id = mpiplan-> comm_id;

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
    istart[j][d2] = mpiplan->grid2->st[d2][j]; iend[j][d2] = mpiplan->grid2->en[d2][j];
  }    

  tmpdims = trplan->grid1->Ldims;
  size_t size0=tmpdims[0]*tmpdims[1]*tmpdims[2]*sizeof(Type1);

  tmpdims = trplan->grid2->Ldims;
  size_t size1 = tmpdims[0]*tmpdims[1]*tmpdims[2]*sizeof(Type2);

#ifdef CUDA
  if(trplan->OutLoc != LocHost)
    printf("Error in trans_MPIplan::pack_sendbuf_trans: expected LocHost for outLoc\n");
#endif

  if(OW && trplan->dt2 <= trplan->dt1)  // CC or C2R
    buf = (Type2 *) src;
  else {
    //    checkCudaErrors(cudaMallocHost(&buf,size1*sizeof(Type2)));
    //buf_alloc = true;
    buf = (Type2 *) tmpbuf;
    tmpbuf += size1;
  }

#ifdef DEBUG
  char str[80];
  sprintf(str,"pack_send_trans.in%d.%d",cnt_pack,mpiplan->Pgrid->taskid);
  //  inv_mo(mo2,imo2);
  // Optimize in the future
  write_buf<Type2>((Type2 *)src,str,dims1,imo1);
#endif

#ifdef CUDA
  cudaStream_t *stream;

  //  printf("offset1=%d %d, mysize1=%d %d\n",trplan->offset1[0],trplan->offset1[1],trplan->mysize1[0],trplan->mysize1[1]);
  if(InLoc == LocHost) {
    /*
    if(devbuf != NULL)
      trplan->DevBuf = devbuf;
    else {
      checkCudaErrors(cudaMalloc(&devbuf,size0*sizeof(Type1)));
#ifdef DEBUG
      printf("pack_sendbuf_trans: Allocated DevBuf size %d\n",size0*sizeof(Type1));
#endif
      DevAlloc = true;
    }
    */
    for(i=0;i<nslices;i++) {
      stream = &(streams[i]);
      //#ifdef DEBUG
	// printf("pack_sendbuf_trans: offset = %d,slice=%d\n",trplan->offset1[i],i);
      //#endif
      //checkCudaErrors(cudaMemcpy(trplan->DevBuf+trplan->offset1[i]*sizeof(Type1),src + trplan->offset1[i]*sizeof(Type1),trplan->mysize1[i]*sizeof(Type1),cudaMemcpyHostToDevice));
      checkCudaErrors(cudaMemcpyAsync(devbuf+trplan->offset1[i]*sizeof(Type1),src + trplan->offset1[i]*sizeof(Type1),trplan->mysize1[i]*sizeof(Type1),cudaMemcpyHostToDevice,*stream));
      checkCudaErrors(cudaEventRecord(trplan->EVENT_H2D,*stream));
    }
    event_hold = &(trplan->EVENT_H2D);
  }
  /*
  if(size1*sizeof(Type2) > size0*sizeof(Type1)) {
    checkCudaErrors(cudaMalloc(&(trplan->DevBuf2),size1*sizeof(Type2)));
    DevAlloc2 = true;
  }
  else
    trplan->DevBuf2 = trplan->DevBuf;
  */

  //  checkCudaErrors(cudaMalloc(&(trplan->DevBuf2),size1*sizeof(Type2)));
  //  checkCudaErrors(cudaMemcpyAsync(buf,trplan->DevBuf2, 65536,cudaMemcpyDeviceToHost,streams[0]));
#endif

#ifdef TIMERS
  double t1,t2;
  t1 = MPI_Wtime();
#endif


#ifdef CUDA
  for(i=0;i<nslices;i++)
    trplan->exec_slice(devbuf,(char *) buf,dim_deriv,i,nslices,event_hold,OW,devbuf+size0);
  //  checkCudaErrors(cudaEventRecord(EVENT_D2H));
  checkCudaErrors(cudaEventSynchronize(trplan->EVENT_D2H));
  //  checkCudaErrors(cudaFree(trplan->DevBuf2));
  // if(DevAlloc2) 
  //  checkCudaErrors(cudaFree(trplan->DevBuf2));
  //if(DevAlloc)
  //  checkCudaErrors(cudaFree(trplan->DevBuf));
#else
  for(i=0;i<nslices;i++)
    trplan->exec_slice(src,(char *) buf,dim_deriv,i,nslices,event_hold,OW,tmpbuf);
#endif

#ifdef TIMERS
  t2 = MPI_Wtime();
  timers.packsend_trans -= t2-t1;
  timers.trans_exec += t2-t1;
#endif

#ifdef DEBUG
  sprintf(str,"pack_send_trans.out%d.%d",cnt_pack++,mpiplan->Pgrid->taskid);
  write_buf<Type2>((Type2 *) buf,str,tmpdims,imo2);
#endif
  
  for(i=0;i<3;i++)
    d[i] = dims1[imo2[i]];

  for(i=0;i < nt;i++) {
    xs = istart[i][imo2[0]];
    xe = iend[i][imo2[0]];
    ys = istart[i][imo2[1]];
    ye = iend[i][imo2[1]];
    p1 = sendbuf + *(SndStrt+i)/sizeof(Type2);
    p0 = buf + istart[i][imo2[0]];  // Assume MPIplan doesn't change mem. order
    for(z=istart[i][imo2[2]];z < iend[i][imo2[2]];z++) {
      pz = p0 + d[0]*d[1]*z;
      for(y=ys;y< ye;y++) {
	p2 = pz + d[0]*y;
	for(x=xs;x < xe;x++)
	  *p1++ = *p2++;;
      }
    }
}
  /*
  if(buf_alloc)
#ifdef CUDA
    cudaFreeHost(buf);
#else
  delete [] buf;
#endif
  */
}

  // Execute 1D transform, combined with local transpose as needed
// Input: in[dims1[imo1[0]]][dims1[imo1[1]]][dims1[imo1[2]]]
// Output: out[dims2[imo2[0]]][dims2[imo2[1]]][dims2[imo2[2]]]

// Execute transform (combined with local transpose if needed), on a slice of data asynchronously via streams
// Expect the calling function to set up and initilize device buffer (asynchronously)
  template <class Type1,class Type2> void transplan<Type1,Type2>::exec_slice(char *in_,char *out_, int dim_deriv,int slice,int nslices,event_t *event_hold, bool OW,char *tmpbuf)
{
  int L,N,m,mocurr[3],mc[3],cmpl;
  Type1 *in;
  Type2 *out;
  Type2 *buf=NULL;
  void swap0(int newmo[3],int mo[3],int L, int *next=NULL);

  in = (Type1 *) in_;
  out = (Type2 *) out_;
  size_t size1 = dims1[0]*dims1[1]*dims1[2] * sizeof(Type1);
  size_t size2 = dims2[0]*dims2[1]*dims2[2] * sizeof(Type2);
    
  if(trans_type->is_empty) {
#ifdef CUDA
    if(InLoc == LocDevice && OutLoc == LocDevice) {
      int d1[3];
      for(int i=0;i<3;i++)
	d1[mo1[i]] = dims1[i];
#ifdef TIMERS
      timers.reorder_out -= MPI_Wtime();
#endif
#ifdef CUTENSOR
      ro_cutensor_out((Type2 *) in,out,mc,d1,slice,nslices,streams[slice]);
#endif
#ifdef TIMERS
      timers.reorder_out += MPI_Wtime();
#endif
      return;
    }
    else if(InLoc == LocHost && OutLoc == LocHost) 
#endif
      if((void *) in == (void *) out)
	reorder_in(in,mo1,mo2,dims1,tmpbuf);
      else
	reorder_out((Type2 *) in,out,mo1,mo2,dims1);
    
    return;
  }

#ifdef CUDA
  //  if(useCuda) {
  bool DevAlloc = false;
  bool DevAlloc2 = false;
  cudaStream_t *stream;

  if(OutLoc == LocHost) {
    out = (Type2 *) tmpbuf;
    tmpbuf += size2;
  }

  //  in = (Type1 *) DevBuf;
  //out = (Type2 *) DevBuf2;

  /*
  if(InLoc == LocHost) {  // Need to transfer data from host to device
    if(DevBuf == NULL) { // Highly unlikely: expect the higher level function to allocate this buffer and initiate transfers
      size_t size1 = dims1[0]*dims1[1]*dims1[2] * sizeof(Type1);
      size_t size2 = dims2[0]*dims2[1]*dims2[2] * sizeof(Type2);
      if(size2 <= size1 || OutLoc == LocDevice)
	DevBuf = out_;
      else {
	size_t maxsize = max(size1,size2);
	checkCudaErrors(cudaMalloc((&(DevBuf)),maxsize));
	DevAlloc = true;
      }
      stream = &(streams[slice]);
      checkCudaErrors(cudaMemcpyAsync(DevBuf+offset1[slice]*sizeof(Type1),in_ + offset1[slice]*sizeof(Type1),mysize1[slice]*sizeof(Type1),cudaMemcpyHostToDevice,*stream));
      checkCudaErrors(cudaEventRecord(EVENT_H2D,*stream));
    }
  
    in = (Type1 *) DevBuf;
  }
  if(OutLoc == LocHost)  { // If no device space is given as input, allocate device space
    size_t size1 = dims1[0]*dims1[1]*dims1[2] * sizeof(Type1);
    size_t size2 = dims2[0]*dims2[1]*dims2[2] * sizeof(Type2);
    if(DevBuf == NULL) {
      if(size2 <= size1)
	DevBuf = (char *) in;
      else {
	checkCudaErrors(cudaMalloc((&(DevBuf)),size2));
	DevAlloc = true;
      }
    }
    if(size2 > size1) {
      checkCudaErrors(cudaMalloc((&(DevBuf2)),size2));
      DevAlloc2 = true;
      out = (Type2 *) DevBuf2;
    }
    else
      out = (Type2 *) DevBuf;
      }*/
#endif

#ifdef DEBUG
  int taskid;
  MPI_Comm_rank(grid1->Pgrid->mpi_comm_glob,&taskid);
#ifdef CUDA
    printf("%d: In transplan::exec, mo1=%d %d %d, mo2=%d %d %d, in=%ld, out=%ld, out=%ld\n, inloc=%d, outloc=%d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2],(long int ) in, (long int) out, (long int ) out_,InLoc,OutLoc);
#else
    printf("%d: In transplan::exec, mo1=%d %d %d, mo2=%d %d %d, in=%ld, out=%ld, out_=%ld\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2],(long int ) in, (long int) out, (long int ) out_);
#endif
#endif



  L = trans_dim;

  if(mo1[L] == 0 || mo2[L] == 0) {

    if(!arcmp(mo1,mo2,3)) { // If there is no change in memrory mapping, there is no need to do transpose. Just call 1D transform.
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      if((void *) in != (void *) out)
	(*(trans_type->exec))(plan->libplan_out[slice],in+offset1[slice],out+offset2[slice]);
      else if(dt2 > dt1) {
	(*(trans_type->exec))(plan->libplan_out[slice],in+offset1[slice],tmpbuf);
#ifdef CUDA
	checkCudaErrors(cudaMemcpyAsync(out+offset2[slice],tmpbuf, size2, cudaMemcpyDeviceToDevice,streams[slice]));
#else
	memcpy(out+offset2[slice],tmpbuf, size2);
#endif 	
      }
      else
	(*(trans_type->exec))(plan->libplan_in[slice],in+offset1[slice],out+offset2[slice]);

      if(dim_deriv == L) {
        int sdims[3]; // Find storage dimensions
        for(int i=0;i<3;i++)
          sdims[mo1[i]] = grid2->Ldims[i];
        sdims[2] = mysize2[slice]/(sdims[0]*sdims[1]);
#ifdef CUDA
        //      compute_deriv_loc_cu<<<>>>(out+offset2[slice],out+offset2[slice],sdims);
#else
	compute_deriv_loc(out+offset2[slice],out+offset2[slice],sdims);
#endif
      }
#ifdef CUDA
      checkCudaErrors(cudaEventRecord(EVENT_EXEC,streams[slice]));
#endif
      cmpl = SLICE;

#ifdef TIMERS
      timers.trans_exec += MPI_Wtime() -t1;
#endif
    }
    else { // Otherwise, combine transpose with transform
#ifdef TIMERS
      double t1=MPI_Wtime();
#endif
      cmpl = reorder_trans_slice(in,out,mo1,mo2,dim_deriv,slice,nslices,event_hold ,OW,tmpbuf);
#ifdef TIMERS
      timers.reorder_trans += MPI_Wtime() -t1;
#endif
    }	
  }
  
  else {  //If first dimension is not the one to be transformed, need to reorder first, combined with transform

    swap0(mocurr,mo1,L);
    size_t size = dims2[0]*dims2[1]*dims2[2]*sizeof(Type2);
    //    checkCudaErrors(cudaMalloc(reinterpret_cast<void **>(&buf), size));
    buf = (Type2 *) tmpbuf;
    tmpbuf += size;
    //}
    //    else

#ifdef TIMERS
    double t1=MPI_Wtime();
#endif

    cmpl = reorder_trans_slice(in,buf,mo1,mocurr,dim_deriv,slice,nslices,event_hold,OW,tmpbuf);

#ifdef TIMERS
    timers.reorder_trans += MPI_Wtime() -t1;
    t1=MPI_Wtime();
#endif

    if(cmpl == SLICE) reorder_out_slice(buf,out,mocurr,mo2,dims2,slice,nslices);
    else if(cmpl == FULL) reorder_out(buf,out,mocurr,mo2,dims2);

#ifdef TIMERS
    timers.reorder_out += MPI_Wtime() -t1;
#endif

    /*
#ifdef CUDA
    //    if(useCuda)
    checkCudaErrors(cudaFree(buf));
    //else
#else
    delete [] buf;

#endif
    */

  }
#ifdef CUDA
  //  if(useCuda) {
  if(OutLoc == LocHost) {
#ifdef TIMERS
    timers.gpu_transfer -= MPI_Wtime();
#endif
    if(cmpl == FULL && slice == nslices-1) {
      size_t size2 = dims2[0]*dims2[1]*dims2[2]*sizeof(Type2);
      cudaDeviceSynchronize();
      checkCudaErrors(cudaMemcpy(out_,out, size2, cudaMemcpyDeviceToHost)); 
    }
    else if(cmpl == SLICE) {
      checkCudaErrors(cudaMemcpyAsync(out_+offset2[slice]*sizeof(Type2),out+offset2[slice], mysize2[slice]*sizeof(Type2),cudaMemcpyDeviceToHost,streams[slice]));
      //	 checkCudaErrors(cudaMemcpy(out_+offset2[slice]*sizeof(Type2),out+offset2[slice], mysize2[slice]*sizeof(Type2),cudaMemcpyDeviceToHost));
      checkCudaErrors(cudaEventRecord(EVENT_D2H,streams[slice]));
      //      printf("%d: offset2=%d %d, mysize2=%d %d\n",taskid,offset2[0],offset2[1],mysize2[0],mysize2[1]);
    }
//     cudaDeviceSynchronize();
#ifdef TIMERS
    timers.gpu_transfer += MPI_Wtime();
#endif
    /*
    if(DevAlloc)
      checkCudaErrors(cudaFree(DevBuf));
    if(DevAlloc2)
      checkCudaErrors(cudaFree(DevBuf2));
    */
  }
#endif
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

template class transform3D<float,float>;
template class transform3D<double,double>;
template class transform3D<mycomplex,float>;
template class transform3D<complex_double,double>;
template class transform3D<float,mycomplex>;
template class transform3D<double,complex_double>;
template class transform3D<mycomplex,mycomplex>;
template class transform3D<complex_double,complex_double>;

template class transplan<float,float>;
template class transplan<double,double>;
template class transplan<mycomplex,float>;
template class transplan<complex_double,double>;
template class transplan<float,mycomplex>;
template class transplan<double,complex_double>;
template class transplan<mycomplex,mycomplex>;
template class transplan<complex_double,complex_double>;


}
