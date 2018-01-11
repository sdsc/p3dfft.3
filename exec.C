#include "p3dfft.h"

namespace p3dfft {

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
#ifdef DEBUG
    printf("stage %d, kind=%d\n",stage_cnt++,curr_stage->kind);
#endif  

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
	else if(!(st->inplace) || (!OW && buf[curr] == buf_in) || size2*dt_2 > size1 *dt_1) {
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
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
#ifdef DEBUG
      	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
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
#ifdef DEBUG
	  printf("%d: exec: Allocating new var %d, dims1=(%d %d %d), dims2=(%d %d %d)\n",taskid,size2,curr_stage->dims1[0],curr_stage->dims1[1],curr_stage->dims1[2],curr_stage->dims2[0],curr_stage->dims2[1],curr_stage->dims2[2]);
#endif
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

  //  rel_change(mo1,mo2,mc);
  //int newL = mc[L];

  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

#ifdef DEBUG
  printf("%d: In transplan::exec, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  if(mo1[L] == 0 || mo2[L] == 0) {

    if(!arcmp(mo1,mo2,3)) {
      (*(trans_type->exec))(lib_plan,in,out);
    }
    else {
      reorder_trans(in,out,mo1,mo2,dims1);   
    }	
  }

  else { //If first dimension is not the one to be transformed, need to reorder first, combined with transform
    swap0(mocurr,mo1,L);
  //    Try to organize so out of place transpose is used, if possible
    //    if(mocurr[0] != mo2[0] || mocurr[1] != mo2[1] || mocurr[2] != mo2[2]) {
    buf = new Type2[dims2[0]*dims2[1]*dims2[2]];
    reorder_trans(in,buf,mo1,mocurr,dims1);   
    /*
    int currdims[3],i,tmpdims[3];
    for(i=0; i < 3;i++) 
      currdims[mocurr[i]] = dims2[mo2[i]];
    */
    reorder_out(buf,out,mocurr,mo2,dims2);
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

#define TRANS_IN 0
#define TRANS_OUT 1


// Input: in[dims1[mo1[0]]][dims1[mo1[1]]][dims1[mo1[2]]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_trans(Type1 *in,Type2 *out,int *mo1,int *mo2,int *dims1)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);
  //  Type1 *pin,*pin1,*p1;
  //Type2 *p2,*ptran2;
  //  libplan plan;

  int imo1[3],imo2[3],d1[3],d2[3];
  int scheme;
  
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

#ifdef DEBUG
  printf("%d: In transplan::reorder_trans, mo1=%d %d %d, mo2=%d %d %d\n",taskid,mo1[0],mo1[1],mo1[2],mo2[0],mo2[1],mo2[2]);
#endif

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  for(i=0;i<3;i++) {
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }
  
    if(dims1[trans_dim] != *inembed) {
      printf("Error in reorder_trans: leading dimension on input %d doesn't match transform type %d\n",dims1[trans_dim],*inembed);
      return;
    }
    if(dims2[trans_dim] != *onembed) {
      printf("Error in reorder_trans: leading dimension on output %d doesnt match transform type %d\n", dims2[trans_dim],*onembed);
      return;
    }

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
  /*
  sprintf(str,"reorder_trans.in%d",cnt_reorder_trans);
  write_buf<Type1>(in,str,dims1,imo1);
  */
#endif


  if(scheme == TRANS_IN) { // Scheme is transform before reordering

    //    if(sizeof(Type1) == sizeof(Type2)) { // Can assume the transform doesn't change dimensions,
    // so we can do in-place transforms

    Type2 *tmp;
    Type2 *pin,*pin1,*pin2,*pout,*pout1,*ptran2,*pout2; 

   switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
	
	  //      m /= dims1[2]; // Need finer grain transform, to reuse cache
	  //Either find existing plan with suitable parameters, or reate a new one
	  //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	tmp = new Type2[d2[0]*d2[1]];
	
	for(k=0;k <d1[2];k++) {
	  (*(trans_type->exec))(lib_plan,in+k*d1[0]*d1[1],tmp);
	  pout = out+k*d2[0]*d2[1];
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

	  pout = out +k*d2[0]*d2[1];
	  for(j=0;j < d2[1];j++) {
	    pin = tmp + j;
	    for(i=0;i < d2[0];i++) {
	      *pout++ = *pin;
	      pin += d2[1];
	    }	
	 
	  }
#endif
	  //      if(!trans_type->is_empty)
	}
	
	delete [] tmp;
	break;
	
      case 2: //1,2,0
	int nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb13 < 1) nb13 = 1;
	//m = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);	
	tmp = new Type2[d2[2]*d2[1]*nb31];
	//double *tmp1=(double *) tmp; 

	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  (*(trans_type->exec))(lib_plan,in+k*d1[0]*d1[1],tmp);
	  for(i=0;i < d2[1];i+=nb13) {
	    i2 = min(i+nb13,d2[1]);
	    pout2 =  out + i*d2[0];
	    pin2 = tmp +i;
	    for(kk=k; kk < k2; kk++) {
	      pin1 = pin2;
	      pout1 = pout2+kk;
	      for(j=0;j < d1[1];j++) {
		pin = pin1  ;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  *pout = *pin++;
		  pout += d2[0];
		}
		pin1 += d2[1];
		pout1+= d2[0]*d2[1];
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
	int nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb13 < 1) nb13 = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	tmp = new Type2[d2[0]*d2[1]*d2[2]];
	(*(trans_type->exec))(lib_plan,in,tmp);
	
	for(k=0;k <d2[0];k+=nb31) {
	  k2 = min(k+nb31,d2[0]);
	  for(i=0;i < d2[2];i+=nb13) {
	    i2 = min(i+nb13,d2[2]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = tmp + kk*d2[2]*d2[1] +i;
	      pout1 = out + kk + i *d2[1]*d2[0];
	      for(j=0;j < d2[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  *pout = *pin++;
		  pout += d2[0]*d2[1];
		}
		pin1 += d2[2];
		pout1 += d2[0];
	      }
	    }
	  }
	}
	
	delete [] tmp;
	break;
      case 0: //2,0,1
	int nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb23 < 1) nb23 = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	
	tmp = new Type2[d2[0]*d2[1]*d2[2]];
	(*(trans_type->exec))(lib_plan,in,tmp);
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = tmp +kk*d2[1]*d2[2] +j*d2[1];
	      pout1 = out +kk +j*d2[0]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d2[1];i++) {
		  *pout =  *pin++;
		  pout += d2[0];
		}
		pin1 += d2[1];
		pout1 += d2[0]*d2[1];
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
	int nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb23 < 1) nb23 = 1;
	//m = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	Type1 *pin;
	Type2 *pout;
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      for(jj=j; jj < j2; jj++) {
		pin = in +kk*d1[0]*d1[1] +jj*d1[0];
		pout = out + kk * d2[0] + jj * d2[0]*d2[1];
		(*(trans_type->exec))(lib_plan,pin,pout);
	      }
	    }
	  }
	}    
      }
      
      break;
    }
    
  }
  else { // Scheme is transform after reordering

    Type1 *tmp;
    Type1 *pin,*pout,*pin1,*pout1,*ptran;
    Type2 *pout2,*ptran2;

    switch(mc[0]) {
    case 1:
      switch(mc[1]) {
      case 0: //1,0,2
	
      //      m /= dims1[2]; // Need finer grain transform, to reuse cache
    //Either find existing plan with suitable parameters, or reate a new one
      //      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
  
	tmp = new Type1[d1[0]*d1[1]];
	for(k=0;k <d1[2];k++) {
	  pin = in + k*d1[0]*d1[1];
#ifdef MKL_BLAS
	  blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmp,d1[1]);
#else
	  for(j=0;j < d1[1];j++) {
	    pout = tmp +j;
	    for(i=0;i < d1[0];i++) {
	      *pout = *pin++;
	      pout += d1[1];
	    }	
	  }
#endif
	  pout2 = out + k*d2[0]*d2[1];
	  (*(trans_type->exec))(lib_plan,tmp,pout2);
	}
	delete [] tmp;

	break;

      case 2: //1,2,0
	int nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb13 < 1) nb13 = 1;
	//m = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);

	tmp = new Type1[d1[0]*d1[2]]; // tmp[k][i] = in[i][j][k]

	for(j=0;j < d1[1];j++) {
	  pin1 = in +j*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + i + kk*d1[0]*d1[1];
		pout = tmp + i*d1[2] +kk;
		
		for(ii=i; ii < i2; ii++) {
		  *pout = *pin++;
		  pout += d1[2];
		}
	      }
	    }
	  }

	  pout2 = out + j*d2[0]*d2[1];
	  (*(trans_type->exec))(lib_plan,tmp,pout2);
	}
	
	delete [] tmp;
	break;
      }
      
      break;
    case 2:
      switch(mc[1]) {
      case 1: //2,1,0
	int nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb31 < 1) nb31 = 1;
	int nb13 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb13 < 1) nb13 = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	
	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  for(i=0;i < d1[0];i+=nb13) {
	    i2 = min(i+nb13,d1[0]);
	    for(kk=k; kk < k2; kk++) {
	      pin1 = in + kk*d1[0]*d1[1] +i;
	      pout1 = tmp + kk +i*d1[2]*d2[1];
	      for(j=0;j < d1[1];j++) {
		pin = pin1;
		pout = pout1;
		for(ii=i; ii < i2; ii++) {
		  *pout = *pin++;
		  pout += d1[2]*d2[1];
		}
		pin1 += d1[0];
		pout1 += d1[2];
	      }
	    }
	  }
	}
	//      if(!trans_type->is_empty)
	(*(trans_type->exec))(lib_plan,tmp,out);
	delete [] tmp;
	
	break;
      case 0: //2,0,1
	int nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb23 < 1) nb23 = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	
	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++){
	      pin1 = in +kk*d1[0]*d1[1] +j*d1[0];
	      pout1 = tmp +kk +j*d1[2]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		pin = pin1;
		pout = pout1;
		for(i=0;i < d1[0];i++) {
		  *pout =  *pin++;
		  pout += d1[2];
		}
		pin1 += d1[0];
		pout1 += d1[2]*d2[1];
	      }
	    }
	  }
	}
	
	(*(trans_type->exec))(lib_plan,tmp,out);
	delete [] tmp;
	
	break;
      }
      break;
    case 0: //0,2,1
      if(mc[1] == 2) {
	int nb32 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb32 < 1) nb32 = 1;
	int nb23 = CACHE_BL / (sizeof(Type1)*d2[0]*d2[1]);
	if(nb23 < 1) nb23 = 1;
	//m = 1;
	//      plan = find_plan(N,m,.true.,dt1,dt2,prec,1,str1,1,str2,doplan);
	
	for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	    j2 = min(j+nb23,d1[1]);
	    for(kk=k; kk < k2; kk++) {
	      pin = in +kk*d1[0]*d1[1] +j*d1[0];
	      ptran2 = out +kk*d2[0] +j*d2[0]*d2[1];
	      for(jj=j; jj < j2; jj++) {
		if(trans_type->exec != NULL)
		  (*(trans_type->exec))(lib_plan,pin,ptran2);
		pin += d1[0];
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
  sprintf(str,"reorder_trans.out%d",cnt_reorder_trans++);
  write_buf<Type2>(out,str,dims2,imo2);
#endif  
}


// Input: in[dims1[mo1[0]]][dims1[mo1[1]]][dims1[mo1[2]]]
// dims2[mo2[i]] = dims1[mo1[i]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type1,class Type2> void transplan<Type1,Type2>::reorder_out(Type2 *in,Type2 *out,int mo1[3],int mo2[3],int *dims_init)
{
  int mc[3],i,j,k,ii,jj,kk,i2,j2,k2;
  void rel_change(int *,int *,int *);
  Type2 *pin,*pout,*pin1,*pout1;
  int imo1[3],imo2[3],d1[3],d2[3];


  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  for(i=0;i<3;i++) {
    d1[i] = dims_init[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }

  rel_change(imo1,imo2,mc);
  pin = in;
  //  pout = pin+1;
  //*pout = *pin;
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2
      for(k=0;k <d1[2];k++) {
	pout1 = out +k*d1[0]*d1[1];
	for(j=0;j < d1[1];j++)
	  pout = pout1 + j;
	  for(i=0;i < d1[0];i++) {
	    *(pout)  = *(pin++);
	    pout += d2[0];
	  }
      }
      break;

    case 2: //1,2,0
      int nb31 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb13 < 1) nb13 = 1;
      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*d1[0]*d1[1] +i;
	    pout1 = out + kk*d2[1] +i*d2[0]*d2[1];
	    for(j=0;j < d1[1];j++) {
	      pin = pin1;
	      pout = pout1;
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += d2[0]*d2[1];
	      }
	      pin1 += d1[0];
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
      int nb31 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb13 < 1) nb13 = 1;
      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in + kk*d1[0]*d1[1] +i;
	    pout1 = out + kk + i * d2[0]*d2[1];
	    for(j=0;j < d1[1];j++) {
	      pin = pin1;
	      pout = pout1;
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += d2[0]*d2[1];
	      }
	      pin1 += d1[0];
	      pout1 += d2[0];
	    }
	  }
	}
      }
      
      break;
    case 0: //2,0,1
      int nb32 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type2)*d2[0]*d2[1]);
      if(nb23 < 1) nb23 = 1;
      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = in +kk*d1[0]*d1[1];
	    pout1 = out +kk;
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1 +jj*d1[0];
	      pout = pout1 +jj*d2[0]*d2[1];
	      for(i=0;i < d1[0];i++) {
		*pout = *pin++;
		pout += d2[0];
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
      int nb32 = CACHE_BL / (sizeof(Type2)*d1[0]*d1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type2)*d2[0]*d2[1]);
      if(nb23 < 1) nb23 = 1;
      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = in +kk*d1[0]*d1[1] +j*d1[0];
	    pout1 = out +kk*d2[0] +j*d2[0]*d2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) 
		*pout++ = *pin++;
	      pin += d1[0];
	      pout += d2[0]*d2[1];
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

// Input: in[d1[mo1[0]]][d1[mo1[1]]][d1[mo1[2]]]
// dims2[mo2[i]] = d1[mo1[i]]
// Output: out[dims2[mo2[0]]][dims2[mo2[1]]][dims2[mo2[2]]]
template <class Type> void reorder_in(Type *in,int mo1[3],int mo2[3],int d1[3])
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
      tmp = new Type[d1[0]+pad];
      for(k=0;k <d1[2];k++) {
	pout = tmp;
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout  = *in++;
	  pout+= pad;//Cache shift
	}
	pin = tmp;
	pout = out +k*d1[0]*d1[1];
	for(j=0;j < d1[1];j++)
	  for(i=0;i < d1[0];i++) {	    
	    *pout  = *pin++;
	    pout += d1[0]+pad;
	  }
	delete [] tmp;
      }
      break;

    case 2: //1,2,0
      int nb31 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type)*dims2[0]*d1[1]);
      if(nb13 < 1) nb13 = 1;
      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}
      
      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*d1[0]*d1[1];
	    pout1 = out + kk*dims2[1];
	    for(j=0;j < d1[1];j++) {
	      pin = pin1+i;
	      pout = pout1+i*dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += dims2[0]*dims2[1];
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
      int nb31 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      if(nb31 < 1) nb31 = 1;
      int nb13 = CACHE_BL / (sizeof(Type)*d2[0]*d2[1]);
      if(nb13 < 1) nb13 = 1;
      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}

      for(k=0;k <d1[2];k+=nb31) {
	k2 = min(k+nb31,d1[2]);
	for(i=0;i < d1[0];i+=nb13) {
	  i2 = min(i+nb13,d1[0]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp + kk*d1[0]*d1[1];
	    pout1 = out + kk;
	    for(j=0;j < d1[1];j++) {
	      pin = pin1 + i;
	      pout = pout1 + i * dims2[0]*dims2[1];
	      for(ii=i; ii < i2; ii++) {
		*pout = *pin++;
		pout += d1[0]*dims2[1];
	      }
	      pin1 += d1[0]+1;
	      pout1 += dims2[0];
	    }
	  }
	}
      }

      delete [] tmp;
      
      break;
    case 0: //2,0,1
      int nb32 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++){
	  for(i=0;i < d1[0];i++) 
	    *pout  = *in++;
	  pout++;//Cache shift
	}

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++){
	    pin1 = tmp +kk*d1[0]*d1[1] +j*d1[0];
	    pout1 = out +kk +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) {
		*pout = *pin++;
		pout += dims2[1];
	      }
	      pin1 += d1[0]+1;
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
      int nb32 = CACHE_BL / (sizeof(Type)*d1[0]*d1[1]);
      if(nb32 < 1) nb32 = 1;
      int nb23 = CACHE_BL / (sizeof(Type)*dims2[0]*dims2[1]);
      if(nb23 < 1) nb23 = 1;
      tmp = new Type[(d1[0]+1)*d1[1]*d1[2]];
      pout = tmp;
      for(k=0;k <d1[2];k++) 
	for(j=0;j < d1[1];j++) {
	  for(i=0;i < d1[0];i++) 
	    *pout  = *in++;
	  pout++; //Cache shift
	}

      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    pin1 = tmp +kk*d1[0]*d1[1] +j*d1[0];
	    pout1 = out +kk*dims2[0] +j*dims2[0]*dims2[1];
	    for(jj=j; jj < j2; jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d1[0];i++) 
		*pout++ = *pin++;
	      pin += d1[0]+1;
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
  double *tmp = (double *) sendbuf;
  Type2 *recvbuf = new Type2[tmpdims[0]*tmpdims[1]*tmpdims[2]];
  MPI_Alltoallv(sendbuf,mpiplan->SndCnts,mpiplan->SndStrt,MPI_REAL,recvbuf,mpiplan->RcvCnts,mpiplan->RcvStrt,MPI_REAL,mpiplan->mpicomm);
  tmp = (double *) recvbuf;
  delete [] sendbuf;
  mpiplan->unpack_recvbuf((Type2 *) out,recvbuf);
  delete [] recvbuf;

#ifdef DEBUG
  char str[80];
  sprintf(str,"transmpi.out%d",cnt_trans++);
  int imo2[3];
  inv_mo(mpiplan->mo2,imo2);
  write_buf<Type2>((Type2 *) out,str,tmpdims,imo2);
#endif
}

template <class Type> void write_buf(Type *buf,char *label,int sz[3],int mo[3]) {
  int i,j,k;
  FILE *fp;
  char str[80],filename[80];
  Type *p=buf;
  complex_double *p1;
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  strcpy(filename,label);
  sprintf(str,".%d",taskid);
  strcat(filename,str);
  fp=fopen(filename,"w");
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

#ifdef DEBUG
  char str[80];
  sprintf(str,"pack_send_trans.in%d",cnt_pack);
  //  inv_mo(mo2,imo2);
  // Optimize in the future
  write_buf<Type2>((Type2 *)src,str,dims1,imo1);
#endif
  trplan->exec(src,(char *) buf);
#ifdef DEBUG
  sprintf(str,"pack_send_trans.out%d",cnt_pack++);
  write_buf<Type2>(buf,str,tmpdims,imo1);
#endif
  
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


}
