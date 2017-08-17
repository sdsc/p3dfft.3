#include "p3dfft.h"

int arcmp(int *A,int *B,int N);
void inv_mo(int mo[3],int imo[3]);
static int cnt_pack=0;
static int cnt_trans=0;

template <class Type1,class Type2> void transform3D<Type1,Type2>::exec(Type1 *in,Type2 *out,int OW)
//void transform3D::exec(char *in,char *out,int OW)
{
  int nvar=0;
  //  int nvar_l = 1;
  //int nup = upper_ind_var(nstages,trans_plans);
  char *buf[2],*var[2],*buf_in,*buf_out;
  int next,curr,dt_1,dt_2,nextvar,stage_cnt=0;

  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);


  next = 1; curr = 0;
  buf[curr] = buf_in = (char *) in;
  buf_out = (char *) out;
  dt_1 = dt1;
  nextvar = 0;

  for(stage *curr_stage=Stages;curr_stage != NULL;curr_stage = curr_stage->next) {
    printf("stage %d, kind=%d\n",stage_cnt++,curr_stage->kind);
  
/*
  for(int nst=0; nst < nstages;nst++) {
    trans_plan *trans_pl = &trans_plans[nst]; 
    MPIplan *mpi_pl = &MPIplans[nst];
    if(nst == nup) {
    }
    else {
    }
*/

/*
    if(!curr_stage->is_set) {
      cout << "Error in transform3D::exec: stage is not set up" << endl;
      return;
    }
*/

    if(curr_stage->kind == TRANS_ONLY) {
      // Only transform, no exchange
	//	transplan<Type1,Type2> *st = (transplan<Type1,Type2> *) curr_stage;
	stage *st = curr_stage;
	int size1 = st->dims1[0] * st->dims1[1] * st->dims1[2]; 
	int size2 = st->dims2[0] * st->dims2[1] * st->dims2[2]; 
	dt_1 = curr_stage->dt1;
	dt_2 = curr_stage->dt2;
	if(!curr_stage->next)
	  buf[next] = buf_out;
	else if(!(st->inplace || (!OW && buf[curr] == buf_in) || size2*dt_2 > size1 *dt_1)) {
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
	  var[nextvar] = new char[size2*dt_2*st->stage_prec];
	  buf[next] = var[nextvar];
	  nvar++;
	  nextvar = 1-nextvar;
	}
	else
	  buf[next] = buf[curr];

	st->myexec(buf[curr],buf[next]);
	dt_1 = dt_2;
      }
    else if(curr_stage->kind == MPI_ONLY) { // Only MPI plan (exchange, no transform)
	int size1 = curr_stage->dims1[0] * curr_stage->dims1[1] * curr_stage->dims1[2]; 
	int size2 = curr_stage->dims2[0] * curr_stage->dims2[1] * curr_stage->dims2[2]; 
	if(!curr_stage->next)
	  buf[next] = buf_out;
	else if((!OW && buf[curr] == buf_in) || size2 > size1) {
	  var[nextvar] = new char[size2*dt_1*curr_stage->stage_prec];
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
	  buf[next] = var[nextvar];
	  nvar++;
	  nextvar = 1-nextvar;
	}
	else 
	  buf[next] = buf[curr];

	curr_stage->myexec(buf[curr],buf[next]);
      }
      else { // MPI and transform combined
	int size1 = curr_stage->dims1[0] * curr_stage->dims1[1] *curr_stage->dims1[2]; 
	int size2 = curr_stage->dims2[0] * curr_stage->dims2[1] *curr_stage->dims2[2]; 
	dt_1 = curr_stage->dt1;
	dt_2 = curr_stage->dt2;
	if(!curr_stage->next)
	  buf[next] = buf_out;
	else if((!OW && buf[curr] == buf_in) || size2*dt_2 > size1*dt_1) {
	  var[nextvar] = new char[size2*dt_2*curr_stage->stage_prec];
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
	  buf[next] = var[nextvar];
	  nvar++;
	  nextvar = 1-nextvar;
	}
	else
	  buf[next] = buf[curr];
	curr_stage->myexec(buf[curr],buf[next]);
	dt_1 = dt_2;
      }
    
    if(nvar > 1) {
      delete [] var[nextvar];
      nvar--;
    }
    next = 1-next;
    curr = 1-curr;
  }
  if(nvar > 0)
    delete [] var[1-nextvar];

}

void stage::myexec(char *in,char *out)
{
  switch(kind) {
  case TRANS_ONLY: 
    if(dt1 == REAL)
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<float,float> *st=(transplan<float,float> *) this;
	  st->exec(in,out);
	}
	else {
	  transplan<double,double> *st=(transplan<double,double> *) this;
	  st->exec(in,out);
	}
      else
	if(stage_prec == 4) {
	  transplan<float,mycomplex> *st=(transplan<float,mycomplex> *) this;
	  st->exec(in,out);
	}
	else {
	  transplan<double,complex_double> *st=(transplan<double,complex_double> *) this;
	  st->exec(in,out);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  transplan<mycomplex,float> *st=(transplan<mycomplex,float> *) this;
	  st->exec(in,out);
	}
	else {
	  transplan<complex_double,double> *st=(transplan<complex_double,double> *) this;
	  st->exec(in,out);
	}
      else
	if(stage_prec == 4) {
	  transplan<mycomplex,mycomplex> *st=(transplan<mycomplex,mycomplex> *) this;
	  st->exec(in,out);
	}
	else {
	  transplan<complex_double,complex_double> *st=(transplan<complex_double,complex_double> *) this;
	  st->exec(in,out);
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
	  st->exec(in,out);
	}
	else {
	  trans_MPIplan<double,double> *st=(trans_MPIplan<double,double> *) (transplan<double,double> *) this;
	  st->exec(in,out);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<float,mycomplex> *st=(trans_MPIplan<float,mycomplex> *) (transplan<float,mycomplex> *) this;
	  st->exec(in,out);
	}
	else {
	  trans_MPIplan<double,complex_double> *st=(trans_MPIplan<double,complex_double> *) (transplan<double,complex_double> *) this;
	  st->exec(in,out);
	}
    else
      if(dt2 == REAL)
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,float> *st=(trans_MPIplan<mycomplex,float> *) (transplan<mycomplex,float> *) this;
	  st->exec(in,out);
	}
	else {
	  trans_MPIplan<complex_double,double> *st=(trans_MPIplan<complex_double,double> *) (transplan<complex_double,double> *) this;
	  st->exec(in,out);
	}
      else
	if(stage_prec == 4) {
	  trans_MPIplan<mycomplex,mycomplex> *st=(trans_MPIplan<mycomplex,mycomplex> *)  (transplan<mycomplex,mycomplex> *) this;
	  st->exec(in,out);
	}
	else {
	  trans_MPIplan<complex_double,complex_double> *st=(trans_MPIplan<complex_double,complex_double> *) (transplan<complex_double,complex_double> *) this;
	  st->exec(in,out);
	}
    break;
  }

}

// Input: in[dims1[mo1[0]]][dims1[mo1[1]]][dims1[mo1[2]]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
//template <class Type1,class Type2> 
template <class Type1,class Type2> void transplan<Type1,Type2>::exec(char *in_,char *out_)
{
  int L,N,m,mocurr[3],mc[3];
  //  lib_plan plan;
  //void (*doexec)(void *,void *);
  Type1 *in;
  Type2 *out;
  Type2 *buf=NULL;
  bool alloc=false;
  void swap0(int newmo[3],int mo[3],int L);

  in = (Type1 *) in_;
  out = (Type2 *) out_;

  L = trans_dim;
  //  prec = sizeof(Type1)/dt1;

  rel_change(mo1,mo2,mc);
  int newL = mc[L];
  /*
  N = dims2[newL];
  m = dims2[0] *dims2[1] *dims2[2] /N;
// ???
  if(trans_type->type_ID == C2RFFT) {
    str1 = N/2+1+pad;
    str2 = N;
  }
  else if(trans_type->type_ID == R2CFFT) {
      N = dims1[L];
      m = dims1[0] *dims1[1] *dims1[2] /N;
      str1 = N+pad;
      str2 = N/2+1;
    }
    else {
      str1 = N+pad;
      str2 = N;
    }
  */

  if(mo2[L] == 0) {

    if(!arcmp(mo1,mo2,3)) {
      memcpy(mocurr,mo1,3*sizeof(int));
      //mocurr = mo1;
      /*
      if(prec == 4) {
	doplan = trans_type->splan;
	doexec = trans_type->sexec;
      }
      else {
	doplan = trans_type->dplan;
	doexec = trans_type->dexec;
      }
      //Either find existing plan with suitable parameters, or create a new one
      plan = find_plan(N,m,inplace,dt1,dt2,prec,1,str1,1,str2,doplan);
      */
      if(!trans_type->is_empty)
	(*(trans_type->exec))(lib_plan,in,out);
    }
    else {
      reorder_trans(in,out,mo1,mo2,dims1);   
      memcpy(mocurr,mo2,3*sizeof(int));
      //      mocurr = mo2;
    }
  }
  else { //If first dimension is not the one to be transformed, need to reorder first, combined with transform
    //mocurr = mo1;
    //    memcpy(mocurr,mo1,3*sizeof(int));
    swap0(mocurr,mo1,L);
	  /*    if(mo1[0] != L) {
      mocurr[0] = mo1[L];
      mocurr[L] = mo1[0];
    }
	  */

  //    Try to organize so out of place transpose is used, if possible
    //    if(mocurr[0] != mo2[0] || mocurr[1] != mo2[1] || mocurr[2] != mo2[2]) {
      buf = new Type2[dims2[0]*dims2[1]*dims2[2]];
      alloc = true;
      //}
      //else
      //buf = out;
    reorder_trans(in,buf,mo1,mocurr,dims1);   
    // }
    //if(mocurr[0] != mo2[0] || mocurr[1] != mo2[1] || mocurr[2] != mo2[2]) {
    int currdims[3],i;
    for(i=0; i < 3;i++)
      currdims[mocurr[i]] = dims2[mo2[i]];
    reorder_out(buf,out,mocurr,mo2,currdims);
    delete [] buf;
  }
}

int arcmp(int *A,int *B,int N)
{
  int i;
  for(i=0; i < N; i++)
    if(A[i] != B[i])
      return(1);

  return(0);
}

// Input: in[dims1[mo1[0]]][dims1[mo1[1]]][dims1[mo1[2]]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);
  Type1 *pin,*pin1,*p1;
  Type2 *pout1,*p2,*ptran2,*pout;
  //  libplan plan;

  rel_change(mo1,mo2,mc);

  printf("In reorder_trans, mc=%d %d %d\n",mc[0],mc[1],mc[2]); 


  if(sizeof(Type1) == sizeof(Type2)) { // Can assume the transform doesn't change dimensions,
    // so we can do in-place transforms

  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2

      //      m /= dims1[2]; // Need finer grain transform, to reuse cache
    //Either find existing plan with suitable parameters, or reate a new one
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
  
      for(k=0;k <dims1[2];k++) {
	p2 = pout = out +k*dims1[0]*dims1[1];
	for(j=0;j < dims1[1];j++)
	  for(i=0;i < dims1[0];i++) {
	    *((Type1 *) pout) = *in++;
	    pout += dims1[0];
	  }	
      if(!trans_type->is_empty)
	(*(trans_type->exec))(plan,p2,p2);
      }

      break;

    case 2: //1,2,0
      int nb31 = CACHE_BL / (sizeof(Type1)*dims1[0]*dims1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type1)*dims2[0]*dims1[1]);
      if(nb13 < 1) nb13 = 1;
      //m = 1;
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
      for(k=0;k <dims1[2];k+=nb31) {
	k2 = min(k+nb31,dims1[2]);
	for(i=0;i < dims1[0];i+=nb13) {
	  i2 = min(i+nb13,dims1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*dims1[0]*dims1[1];

	    ptran2 = pout1 = out + kk*dims2[1];
	    for(j=0;j < dims1[1];j++) {
	      pin = pin1+i;
	      pout = pout1+i*dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*((Type1 *) pout) = *pin++;
		pout += dims2[0]*dims2[1];
	      }
	      pin1 += dims1[0];
	      pout1++;
	    }
	    if(!trans_type->is_empty)
	      for(ii=i; ii < i2; ii++) {
		(*(trans_type->exec))(plan,ptran2,ptran2);
		ptran2+=dims2[0]*dims2[1];
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
      int nb31 = CACHE_BL / (sizeof(Type1)*dims1[0]*dims1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type1)*dims2[0]*dims2[1]);
      if(nb13 < 1) nb13 = 1;
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);

      for(k=0;k <dims1[2];k+=nb31) {
	k2 = min(k+nb31,dims1[2]);
	for(i=0;i < dims1[0];i+=nb13) {
	  i2 = min(i+nb13,dims1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in + kk*dims1[0]*dims1[1];
	    pout1 = out + kk;
	    for(j=0;j < dims1[1];j++) {
	      pin = pin1 + i;
	      pout = pout1 + i * dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*((Type1 *) pout) = *pin++;
		pout += dims1[0]*dims2[1];
	      }
	      pin1 += dims1[0];
	      pout1 += dims2[0];
	    }
	  }
	}
      }
      if(!trans_type->is_empty)
	(*(trans_type->exec))(plan,out,out);
      
      break;
    case 0: //2,0,1
      int nb32 = CACHE_BL / (sizeof(Type1)*dims1[0]*dims1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type1)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);

      for(k=0;k <dims1[2];k+=nb32) {
	k2 = min(k+nb32,dims1[2]);
	for(j=0;j < dims1[1];j+=nb23) {
	  j2 = min(j+nb23,dims1[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = in +kk*dims1[0]*dims1[1] +j*dims1[0];
	    pout1 = out +kk +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < dims1[0];i++) {
		*((Type1 *) pout) =  *pin++;
		pout += dims2[1];
	      }
	      pin1 += dims1[0];
	      pout1 += dims2[0]*dims2[1];
	    }
	  }
	}
      }

      if(!trans_type->is_empty)
	(*(trans_type->exec))(plan,out,out);
      break;
    }
    break;
  case 0: //0,2,1
    if(mc[1] == 2) {
      int nb32 = CACHE_BL / (sizeof(Type1)*dims1[0]*dims1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type1)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      //m = 1;
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);

      for(k=0;k <dims1[2];k+=nb32) {
	k2 = min(k+nb32,dims1[2]);
	for(j=0;j < dims1[1];j+=nb23) {
	  j2 = min(j+nb23,dims1[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*dims1[0]*dims1[1] +j*dims1[0];
	    pout1 = out +kk*dims2[0] +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      ptran2 = pout = pout1;
	      for(i=0;i < dims1[0];i++) 
		*((Type1 *) pout++) =  *pin++;
	      if(trans_type->exec != NULL)
		(*(trans_type->exec))(plan,ptran2,ptran2);
	      pin += dims1[0];
	      pout += dims2[0]*dims2[1];
	    }
	  }
	}
      }      
    }
    else if(!trans_type->is_empty)

      //      plan = find_plan(N,m,(in==out),dt1,dt2,prec,1,str1,1,str2,doplan);
      (*(trans_type->exec))(plan,in,out);
    
    break;
  }
  }
  else { // changing datatype in transform


  }
}


// Input: in[dims1[mo1[0]]][dims1[mo1[1]]][dims1[mo1[2]]]
// dims2[mo2[i]] = dims1[mo1[i]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int *dims_init)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);
  Type2 *pin,*pout,*pin1,*pout1;

  rel_change(mo1,mo2,mc);

  pin = in;
  pout = pin+1;
  *pout = *pin;
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2
      for(k=0;k <dims_init[2];k++) {
	pout = out +k*dims_init[0]*dims_init[1];
	for(j=0;j < dims_init[1];j++)
	  for(i=0;i < dims_init[0];i++) {
	    *(pout)  = *(pin++);
	    pout += dims_init[0];
	  }
      }
      break;

    case 2: //1,2,0
      int nb31 = CACHE_BL / (sizeof(Type2)*dims_init[0]*dims_init[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type2)*dims2[0]*dims_init[1]);
      if(nb13 < 1) nb13 = 1;
      for(k=0;k <dims_init[2];k+=nb31) {
	k2 = min(k+nb31,dims_init[2]);
	for(i=0;i < dims_init[0];i+=nb13) {
	  i2 = min(i+nb13,dims_init[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*dims_init[0]*dims_init[1];
	    pout1 = out + kk*dims2[1];
	    for(j=0;j < dims_init[1];j++) {
	      pin = pin1+i;
	      pout = pout1+i*dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += dims2[0]*dims2[1];
	      }
	      pin1 += dims_init[0];
	      pout1++;
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
      int nb31 = CACHE_BL / (sizeof(Type2)*dims_init[0]*dims_init[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type2)*dims2[0]*dims2[1]);
      if(nb13 < 1) nb13 = 1;
      for(k=0;k <dims_init[2];k+=nb31) {
	k2 = min(k+nb31,dims_init[2]);
	for(i=0;i < dims_init[0];i+=nb13) {
	  i2 = min(i+nb13,dims_init[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in + kk*dims_init[0]*dims_init[1];
	    pout1 = out + kk;
	    for(j=0;j < dims_init[1];j++) {
	      pin = pin1 + i;
	      pout = pout1 + i * dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += dims_init[0]*dims2[1];
	      }
	      pin1 += dims_init[0];
	      pout1 += dims2[0];
	    }
	  }
	}
      }
      
      break;
    case 0: //2,0,1
      int nb32 = CACHE_BL / (sizeof(Type2)*dims_init[0]*dims_init[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type2)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      for(k=0;k <dims_init[2];k+=nb32) {
	k2 = min(k+nb32,dims_init[2]);
	for(j=0;j < dims_init[1];j+=nb23) {
	  j2 = min(j+nb23,dims_init[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = in +kk*dims_init[0]*dims_init[1] +j*dims_init[0];
	    pout1 = out +kk +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < dims_init[0];i++) {
		*pout = *pin++;
		pout += dims2[1];
	      }
	      pin1 += dims_init[0];
	      pout1 += dims2[0]*dims2[1];
	    }
	  }
	}
      }

      break;
    }
    break;
  case 0: //0,2,1
    if(mc[1] == 2) {
      int nb32 = CACHE_BL / (sizeof(Type2)*dims_init[0]*dims_init[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type2)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      for(k=0;k <dims_init[2];k+=nb32) {
	k2 = min(k+nb32,dims_init[2]);
	for(j=0;j < dims_init[1];j+=nb23) {
	  j2 = min(j+nb23,dims_init[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*dims_init[0]*dims_init[1] +j*dims_init[0];
	    pout1 = out +kk*dims2[0] +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < dims_init[0];i++) 
		*pout++ = *pin++;
	      pin += dims_init[0];
	      pout += dims2[0]*dims2[1];
	    }
	  }
	}
      }      
    }
    else
      for(k=0;k<dims_init[2];k++)
	for(j=0;j < dims_init[1];j++)
	  for(i=0;i < dims_init[0];i++)
	    *out++ = *in++;
    break;
  }


}

// Input: in[dims_init[mo1[0]]][dims_init[mo1[1]]][dims_init[mo1[2]]]
// dims2[mo2[i]] = dims_init[mo1[i]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int dims_init[3])
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  Type *pin,*pin1,*pout,*pout1,*tmp;
  void rel_change(int *,int *,int *);

  rel_change(mo1,mo2,mc);

  int pad = CACHEPAD/sizeof(Type);
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2
      tmp = new Type[dims_init[0]+pad];
      for(k=0;k <dims_init[2];k++) {
	pout = tmp;
	for(j=0;j < dims_init[1];j++) {
	  for(i=0;i < dims_init[0];i++) 
	    *pout  = *in++;
	  pout+= pad;//Cache shift
	}
	pin = tmp;
	pout = out +k*dims_init[0]*dims_init[1];
	for(j=0;j < dims_init[1];j++)
	  for(i=0;i < dims_init[0];i++) {	    
	    *pout  = *pin++;
	    pout += dims_init[0]+pad;
	  }
	delete [] tmp;
      }
      break;

    case 2: //1,2,0
      int nb31 = CACHE_BL / (sizeof(Type)*dims_init[0]*dims_init[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type)*dims2[0]*dims_init[1]);
      if(nb13 < 1) nb13 = 1;
      tmp = new Type[(dims_init[0]+1)*dims_init[1]*dims_init[2]];
      pout = tmp;
      for(k=0;k <dims_init[2];k++) 
	for(j=0;j < dims_init[1];j++) {
	  for(i=0;i < dims_init[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}
      
      for(k=0;k <dims_init[2];k+=nb31) {
	k2 = min(k+nb31,dims_init[2]);
	for(i=0;i < dims_init[0];i+=nb13) {
	  i2 = min(i+nb13,dims_init[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*dims_init[0]*dims_init[1];
	    pout1 = out + kk*dims2[1];
	    for(j=0;j < dims_init[1];j++) {
	      pin = pin1+i;
	      pout = pout1+i*dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += dims2[0]*dims2[1];
	      }
	      pin1 += dims_init[0]+1;
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
      int nb31 = CACHE_BL / (sizeof(Type)*dims_init[0]*dims_init[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type)*dims2[0]*dims2[1]);
      if(nb13 < 1) nb13 = 1;
      tmp = new Type[(dims_init[0]+1)*dims_init[1]*dims_init[2]];
      pout = tmp;
      for(k=0;k <dims_init[2];k++) 
	for(j=0;j < dims_init[1];j++) {
	  for(i=0;i < dims_init[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}

      for(k=0;k <dims_init[2];k+=nb31) {
	k2 = min(k+nb31,dims_init[2]);
	for(i=0;i < dims_init[0];i+=nb13) {
	  i2 = min(i+nb13,dims_init[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp + kk*dims_init[0]*dims_init[1];
	    pout1 = out + kk;
	    for(j=0;j < dims_init[1];j++) {
	      pin = pin1 + i;
	      pout = pout1 + i * dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += dims_init[0]*dims2[1];
	      }
	      pin1 += dims_init[0]+1;
	      pout1 += dims2[0];
	    }
	  }
	}
      }

      delete [] tmp;
      
      break;
    case 0: //2,0,1
      int nb32 = CACHE_BL / (sizeof(Type)*dims_init[0]*dims_init[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      tmp = new Type[(dims_init[0]+1)*dims_init[1]*dims_init[2]];
      pout = tmp;
      for(k=0;k <dims_init[2];k++) 
	for(j=0;j < dims_init[1];j++){
	  for(i=0;i < dims_init[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}

      for(k=0;k <dims_init[2];k+=nb32) {
	k2 = min(k+nb32,dims_init[2]);
	for(j=0;j < dims_init[1];j+=nb23) {
	  j2 = min(j+nb23,dims_init[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = tmp +kk*dims_init[0]*dims_init[1] +j*dims_init[0];
	    pout1 = out +kk +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < dims_init[0];i++) {
		*pout = *pin++;
		pout += dims2[1];
	      }
	      pin1 += dims_init[0]+1;
	      pout1 += dims2[0]*dims2[1];
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
      int nb32 = CACHE_BL / (sizeof(Type)*dims_init[0]*dims_init[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      tmp = new Type[(dims_init[0]+1)*dims_init[1]*dims_init[2]];
      pout = tmp;
      for(k=0;k <dims_init[2];k++) 
	for(j=0;j < dims_init[1];j++) {
	  for(i=0;i < dims_init[0];i++) 
	    *pout  = *in++;
	  pout++; //Cache shift
	}

      for(k=0;k <dims_init[2];k+=nb32) {
	k2 = min(k+nb32,dims_init[2]);
	for(j=0;j < dims_init[1];j+=nb23) {
	  j2 = min(j+nb23,dims_init[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*dims_init[0]*dims_init[1] +j*dims_init[0];
	    pout1 = out +kk*dims2[0] +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < dims_init[0];i++) 
		*pout++ = *pin++;
	      pin += dims_init[0]+1;
	      pout += dims2[0]*dims2[1];
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
  pack_sendbuf(sendbuf,(Type *) in);
  Type *recvbuf = new Type[dims2[0]*dims2[1]*dims2[2]];
  MPI_Alltoallv(sendbuf,SndCnts,SndStrt,MPI_REAL,recvbuf,RcvCnts,RcvStrt,MPI_REAL,mpicomm);
  delete [] sendbuf;
  unpack_recvbuf((Type *) out,recvbuf);
  delete [] recvbuf;
}

template <class Type> void MPIplan<Type>::pack_sendbuf(Type *sendbuf,Type *src)
{
  int i,x,y,z;
  Type *p1,*p0;

  int j,istart[numtasks][3],iend[numtasks][3],d[3];
  int imo1[3],imo2[3];

  //  inv_mo(mo1,imo1);
  //inv_mo(mo2,imo2);

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
    p1 = sendbuf + *(SndStrt+i)*4/sizeof(Type);
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
  Type *p1,*p0,*pin,*pin1,*pout,*pout1;

  int j,istart[numtasks][3],iend[numtasks][3],isize[numtasks][3],d[3];
  int xstart,xend,ystart,yend,zstart,zend;
  int xsize,ysize,zsize;
  int imo1[3],imo2[3];

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
      p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
	  p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
	int nb31 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	if(nb13 < 1) nb13 = 1;

	for(i=0;i < numtasks;i++) {
	  p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
		pout1 = p0 + k +x*d[0];
		for(y=0;y < ysize;y++) {
		  pin = pin1;
		  pout = pout1;
		  //		  printf("%d: y,z=%d %d, x from %d to %d\n",taskid,y,k,x,x2);
		  for(ii=x;ii<x2;ii++) {
		    *pout = *pin++;
		    pout += d[0];
		  }
		  pin1 += xsize;
		  pout1+= d[0]*d[1];
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
	int nb31 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	if(nb13 < 1) nb13 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
		    *pout = *pin++;
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
	int nb32 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	if(nb23 < 1) nb23 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
		    *pout = *pin++;
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
	int nb32 = CACHE_BL / (sizeof(Type)*dims1[imo1[0]]*dims1[imo1[1]]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type)*d[0]*d[1]);
	if(nb23 < 1) nb23 = 1;
	
	for(i=0;i < numtasks;i++) {
	  p1 = recvbuf + *(RcvStrt+i)*4/sizeof(Type);
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
template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::exec(char *in,char *out) {
  //  Type1 *in;
  //Type2 *out;
  //in = (Type1 *) in_;
  //out = (Type2 *) out_;
  int *tmpdims;

  tmpdims = trplan->grid2->ldims;

  Type2 *sendbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
  pack_sendbuf_trans(sendbuf,in);
  tmpdims = mpiplan->grid2->ldims;
  Type2 *recvbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
  MPI_Alltoallv(sendbuf,mpiplan->SndCnts,mpiplan->SndStrt,MPI_REAL,recvbuf,mpiplan->RcvCnts,mpiplan->RcvStrt,MPI_REAL,mpiplan->mpicomm);
  delete [] sendbuf;
  mpiplan->unpack_recvbuf((Type2 *) out,recvbuf);
  delete [] recvbuf;
  char str[80];
  sprintf(str,"transmpi.out%d",cnt_trans++);
  int imo2[3];
  inv_mo(mpiplan->mo2,imo2);
  write_buf((Type2 *) out,str,tmpdims,imo2);
}
template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::write_buf(Type2 *buf,char *label,int sz[3],int mo[3]) {
  int i,j,k;
  FILE *fp;
  char str[80],filename[80];
  complex_double *p=(complex_double *) buf;

  strcpy(filename,label);
  sprintf(str,".%d",mpiplan->taskid);
  strcat(filename,str);
  fp=fopen(filename,"w");
  for(k=0;k<sz[mo[2]];k++)
    for(j=0;j<sz[mo[1]];j++)
      for(i=0;i<sz[mo[0]];i++) {
	if(abs(*p) > 1.e-7) {
	  fprintf(fp,"(%d %d %d) %lg %lg\n",i,j,k,p->real(),p->imag());
	}
	p++;
      }
  fclose(fp); 
}


template <class Type1,class Type2> void trans_MPIplan<Type1,Type2>::pack_sendbuf_trans(Type2 *sendbuf,char *src)
{
  int i,nt,x,y,z,l;

  int d1 = mpiplan->d1;
  int d2 = mpiplan->d2;
  int *dims1 = mpiplan->dims1;
  int *mo1 = mpiplan->mo1;
  int *mo2 = mpiplan->mo2;
  int d[3],*tmpdims;
  Type2 *p1,*p0,*buf;
  int imo1[3],imo2[3];

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

  if(trplan->inplace) 
    buf = (Type2 *) src;
  else
    buf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];

  char str[80];
  sprintf(str,"pack_send_trans.in%d",cnt_pack);
  //  inv_mo(mo2,imo2);
  // Optimize in the future
  write_buf((Type2 *) src,str,dims1,imo1);
  trplan->exec(src,(char *) buf);
  sprintf(str,"pack_send_trans.out%d",cnt_pack++);
  write_buf(buf,str,tmpdims,imo1);
  
  for(i=0;i<3;i++)
    d[i] = dims1[imo1[i]];
  
  for(i=0;i < nt;i++) {
    p1 = sendbuf + *(SndStrt+i)*4/sizeof(Type2);
    for(z=istart[i][imo1[2]];z < iend[i][imo1[2]];z++)
      for(y=istart[i][imo1[1]];y< iend[i][imo1[1]];y++) {
	p0 = buf+d[0]*(z*d[1]+y) +istart[i][imo1[0]];
	for(x=istart[i][imo1[0]];x < iend[i][imo1[0]];x++)
	  *p1++ = *p0++;;
      }
  }

  if(!trplan->inplace) 
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

template class transform3D<float,float>;
template class transform3D<double,double>;
template class transform3D<mycomplex,float>;
template class transform3D<complex_double,double>;
template class transform3D<float,mycomplex>;
template class transform3D<double,complex_double>;
template class transform3D<mycomplex,mycomplex>;
template class transform3D<complex_double,complex_double>;

