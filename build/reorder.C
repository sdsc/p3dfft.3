
/*
#ifdef CUDA 
#include "reorder.cu"
#endif
*/

#include "p3dfft.h"

#ifdef CUDA
#ifndef CUTENSOR 
#include "reorder_kernel.cu"
#endif
#endif

//namespace p3dfft {

#ifdef CUBLAS
cublasHandle_t cublas_handle;
#endif 

#ifdef CUDA
#ifdef CUTENSOR
//template <class Type1,class Type2> void transplan<Type1,Type2>::

template <class T> void ro_cutensor_in(T *in,T *out,int imc[3],int dout[3],int slice,int nslices,cudaStream_t stream)
{
  int imo1[] = {0,1,2};
  //int imo2[] = {1,2,0};
int i;

cutensorStatus_t err;
cutensorHandle_t handle;
cutensorInit(&handle);
cudaDataType_t mytype;

if(typeid(T) == type_float)
  mytype = CUDA_R_32F;
else if(typeid(T) == type_double)
  mytype = CUDA_R_64F;
else if(typeid(T) == type_complex)
  mytype = CUDA_C_32F;
else if(typeid(T) == type_complex_double)
  mytype = CUDA_C_64F;

int64_t extA[3],extB[3];
for(i=0;i<3;i++) {
  extB[i] = dout[i];
  extA[imc[i]] = dout[i];
}

// printf("extA = %ld %ld %ld\n",extA[0],extA[1],extA[2]);
// printf("extB = %ld %ld %ld\n",extB[0],extB[1],extB[2]);

 // printf("imo1=%d %d %d, imo2=%d %d %d\n",0,1,2,imo1[2],imo2[0],imo2[1],imo2[2]);

 T one=1.;
cutensorTensorDescriptor_t descA,descB;
 int64_t strideA[3],strideB[3];
 strideA[0] = 1;
 strideA[1] = extA[0];
 strideA[2] = extA[0]*extA[1];
 strideB[0] = 1;
 strideB[1] = extB[0];
 strideB[2] = extB[0]*extB[1];

 // extB[imc[2]] /= nslices;
 //extA[2] /= nslices;

 err = cutensorInitTensorDescriptor( &handle, &descA,3,extA,strideA,mytype,CUTENSOR_OP_IDENTITY);
err = cutensorInitTensorDescriptor( &handle, &descB,3,extB,strideB,mytype,CUTENSOR_OP_IDENTITY);

if(err != CUTENSOR_STATUS_SUCCESS)
  printf("Error while creating tensor descriptor\n");
// cudaDeviceSynchronize();

/*
 long int dispA;
 switch(imc[2]) {
 case 0: dispA=extA[0]*slice;
   break;
 case 1: dispA=extA[0]*extA[1]*slice;
   break;
 case 2:  dispA=extA[0]*extA[1]*extA[2]*slice;
   break;
 }
*/
 // long int dispB = extB[0]*extB[1]*extB[2]*slice;
 err = cutensorPermutation( &handle, &one, in, &descA, imo1, out, &descB, imc, mytype,stream);
 // cudaDeviceSynchronize();

if(err != CUTENSOR_STATUS_SUCCESS)
  printf("Error while permuting tensor: %s\n",cutensorGetErrorString(err));

}

template <class T> void ro_cutensor_out(T *in,T *out,int imc[3],int din[3],int slice,int nslices,cudaStream_t stream)
{
  int imo1[] = {0,1,2};
  //int imo2[] = {1,2,0};
int i;

cutensorStatus_t err;
cutensorHandle_t handle;
cutensorInit(&handle);
cudaDataType_t mytype;

if(typeid(T) == type_float)
  mytype = CUDA_R_32F;
else if(typeid(T) == type_double)
  mytype = CUDA_R_64F;
else if(typeid(T) == type_complex)
  mytype = CUDA_C_32F;
else if(typeid(T) == type_complex_double)
  mytype = CUDA_C_64F;

int64_t extA[3],extB[3];
for(i=0;i<3;i++) {
  extA[i] = din[i];
  extB[i] = din[imc[i]];
}
 int64_t strideA[3],strideB[3];
 strideA[0] = 1;
 strideA[1] = extA[0];
 strideA[2] = extA[0]*extA[1];
 strideB[0] = 1;
 strideB[1] = extB[0];
 strideB[2] = extB[0]*extB[1];

 // extB[2] /= nslices;
 //extA[imc[2]] /= nslices;

// printf("extA = %ld %ld %ld\n",extA[0],extA[1],extA[2]);
// printf("extB = %ld %ld %ld\n",extB[0],extB[1],extB[2]);

 // printf("imo1=%d %d %d, imo2=%d %d %d\n",imo1[0],imo1[1],imo1[2],imo2[0],imo2[1],imo2[2]);

 T one=1.;
cutensorTensorDescriptor_t descA,descB;
 err = cutensorInitTensorDescriptor( &handle, &descA,3,extA,strideA,mytype,CUTENSOR_OP_IDENTITY);
err = cutensorInitTensorDescriptor( &handle, &descB,3,extB,strideB,mytype,CUTENSOR_OP_IDENTITY);

if(err != CUTENSOR_STATUS_SUCCESS)
  printf("Error while creating tensor descriptor\n");

#ifdef DEBUG
 printf("cutensorPermutation: extA=%d %d %d, extB=%d %d %d, imo1=%d %d %d, imc=%d %d %d,mytype=%d\n",extA[0],extA[1],extA[2],extB[0],extB[1],extB[2],imo1[0],imo1[1],imo1[2],imc[0],imc[1],imc[2],mytype);
#endif

 /*
 long int dispA;
 switch(imc[2]) {
 case 0: dispA=extA[0]*slice;
   break;
 case 1: dispA=extA[0]*extA[1]*slice;
   break;
 case 2:  dispA=extA[0]*extA[1]*extA[2]*slice;
   break;
 }
 */
 long int dispB = extB[0]*extB[1]*extB[2]*slice;
 err = cutensorPermutation( &handle, &one, in, &descA, imo1, out, &descB, imc, mytype,stream);
 // cudaDeviceSynchronize();

#ifdef DEBUG
 printf("cutensorPermutation done\n");
#endif

if(err != CUTENSOR_STATUS_SUCCESS)
  printf("Error while permuting tensor: %s\n",cutensorGetErrorString(err));


}

#elif defined CUTT

#define cuttCheck(stmt) do {                                 \
  cuttResult err = stmt;                            \
  if (err != CUTT_SUCCESS) {                          \
    fprintf(stderr, "%s in file %s, function %s\n", #stmt,__FILE__,__FUNCTION__); \
    exit(1); \
  }                                                  \
} while(0)

template <class T> void ro_cutt_in(T *in,T *out,int imc[3],int dout[3])
{
  int imo1[] = {0,1,2};
  int i, extA[3];

  for(i=0;i<3;i++) {
    extA[imc[i]] = dout[i];
  }

  cuttHandle plan;
  //cuttCheck(
  cuttPlan(&plan, 3, extA, imc, sizeof(T), 0);
  cuttExecute(plan, in, out);
  cuttDestroy(plan);
}
template <class T> void ro_cutt_out(T *in,T *out,int imc[3],int din[3])
{
  cuttHandle plan;
  cuttPlan(&plan, 3, din, imc, sizeof(T), 0);
  cuttExecute(plan, in, out);
  cuttDestroy(plan);
}
#endif

#ifdef CUBLAS
template <class Type1,class Type2> void transplan<Type1,Type2>::rot102in_cublas(int lda,int ldb,Type2 *A,Type2 *C)
{
    cublasOperation_t transa = CUBLAS_OP_T;
    cublasOperation_t transb = CUBLAS_OP_N;
    cublasStatus_t	 stat;
    if(typeid(Type2) == type_float) {
       float alpha=1.0;
       float beta=0.0;	
       float *B = new float[lda*ldb];
       stat = cublasSgeam(cublas_handle,transa,transb,ldb,lda,&alpha,(const float *) A,lda,&beta,B,ldb,(float *) C,ldb);
       delete [] B;	
   }
   else if(typeid(Type2) == type_double) {
       double alpha=1.0;
       double beta=0.0;	
       double *B = new double[lda*ldb];

       stat = cublasDgeam(cublas_handle,transa,transb,ldb,lda,&alpha,(const double *) A,lda,&beta,B,ldb,(double *) C,ldb);
       delete [] B;	
   }
   else if(typeid(Type2) == type_complex) {
       mycomplex alpha(1.0,0.0);
       mycomplex beta(0.0,0.0);	
       cuComplex *B = new cuComplex[lda*ldb];

       stat = cublasCgeam(cublas_handle,transa,transb,ldb,lda,(const cuComplex *) &alpha,(const cuComplex *) A,lda,(const cuComplex *) &beta,B,ldb,(cuComplex *) C,ldb);
       delete [] B;	
   }
   else if(typeid(Type2) == type_complex_double) {
       complex_double alpha(1.0,0.0);
       complex_double beta(0.0,0.0);	
       cuDoubleComplex *B = new cuDoubleComplex[lda*ldb];
       delete [] B;	

       stat = cublasZgeam(cublas_handle,transa,transb,ldb,lda,(const cuDoubleComplex *) &alpha,(const cuDoubleComplex *) A,lda,(const cuDoubleComplex *) &beta,B,ldb,(cuDoubleComplex *) C,ldb);
   }

   if(stat != CUBLAS_STATUS_SUCCESS)
      printf("Error in CUBLAS call from rot102in_cublas\n");

}
#endif

#endif


template <class Type1,class Type2> void transplan<Type1,Type2>::rot102in_slice(Type1 *in,Type2 *out,bool inplace,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int slice,int nslices,bool deriv, char *tmpbuf,int pack_dim,int pack_procs)
{

  //    Type2 *tmp;

#if defined CUDA && !defined CUBLAS  

  int sdims[3];
  if(exec)
    //tmp = (Type2 *) in + offset1[slice];
      //else {
      //cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*mysize2[slice]);
    (*(exec))(plan->libplan_out[slice],in+offset1[slice],(Type2 *) tmpbuf+offset2[slice]);

  d2[2] = mysize2[slice]/(d2[0]*d2[1]);

  if(deriv) {
    sdims[0] = d2[1];
    sdims[1] = d2[0];
    sdims[2] = d2[2];
    //      compute_deriv_loc_cu<<<>>>(tmp,tmp,sdims);
  }

  int imc[3] = {1,0,2};

#ifdef CUTENSOR
    ro_cutensor_in<Type2>((Type2 *) tmpbuf+offset2[slice],out+offset2[slice],imc,d2,slice,nslices,streams[slice]);
#elif defined CUTT
	ro_cutt_in<Type2>((Type2 *) tmpbuf+offset2[slice],out+offset2[slice],imc,d2,slice,nslices,streams[slice]);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro102in_cu<Type2>(gridDim,blockSize,(Type2 *) tmpbuf+offset1[slice],out+offset2[slice],inplace,d1,d2,slice,nslices,streams[slice]);
//    rot102in_cu<<<gridDim,blockSize>>>(in,out,inplace,d1,d2,exec,plan,slice,nslice,streams[slice]);
#endif
    //    if(!is_empty)
    //  cudaFree(tmp);
#else
   
    Type2 *pout,*pin,*out1;
     int i,j,k;
     int ipack,mypacksize,mystart,start,myen,l,nl,sz,dmult;
     /*
     if(OW && dt2 <= dt1) 
       inplace=true;
       else {  */
	 /* for(k=0;k <d1[2];k++) 
	   (*(trans_type->exec))(plan->libplan_out,in+k*d1[0]*d1[1],tmp+k*d2[0]*d2[1]);
       }
       else
       */
     //     if(!is_empty) 
       /*
       if((void *) in == (void *) out && dt2 > dt1) 
	 tmp = new Type2[d2[0]*d2[1]*d2[2]];
       else
	 tmp = new Type2[d2[0]*d2[1]];
       */

     int kst = offset1[slice]/(d1[0]*d1[1]);
     int ken = (offset1[slice]+mysize1[slice])/(d1[0]*d1[1]);

     if(pack_dim == 0 || pack_dim == 1) {
       sz = d2[pack_dim]/pack_procs;
       l = d2[pack_dim] % pack_procs;
       nl = pack_procs - l;
       dmult = d2[1-pack_dim];
     }
     for(k=kst ;k < ken;k++) {
       /*
       if(inplace) {
	 tmp = (Type2 *) in + k*d2[0]*d2[1];
	 (*(trans_type->exec))(plan->libplan_in,in+k*d1[0]*d1[1],tmp);
       }
       else */
	 //	 if((void *) in != (void *) out || dt2 <= dt1)
         if(exec == NULL) 
	   tmpbuf = (char *) (in + offset1[slice] + k*d2[0]*d2[1]); 
	 else 
	   /*
	   if((void *) in == (void *) out && dt2 > dt1)  {
	     (*exec)(plan->libplan_out[nslices],in+k*d1[0]*d1[1],(Type2 *) tmpbuf+k*d2[0]*d2[1]);
	     continue;
	     }
	     else*/ 
	     (*exec)(plan->libplan_out[slice],in+k*d1[0]*d1[1],(Type2 *) tmpbuf);
	 
#ifdef MKL_BLAS
	 pout = out+k*d2[0]*d2[1];
	 blas_trans<Type2>(d2[1],d2[0],1.0,(Type2 *) tmpbuf,d2[1],pout,d2[0]);
#else

	 if(pack_dim == 0 || pack_dim == 1) {
	   mystart = start = 0;

	   for(ipack=0;ipack<pack_procs;ipack++) {
	     if(ipack >= nl)
	       mypacksize = sz+1;
	     else
	       mypacksize = sz;
	     start = mystart * dmult*(ken-kst);
	     pout = out+start+(k-kst)*dmult*mypacksize;

	     myen = mystart+mypacksize;
	     if(pack_dim == 1) 
	       for(j=mystart;j < myen;j++) {
		 pin = (Type2 *) tmpbuf + j;
		 for(i=0;i < d2[0];i++) {
		   *pout++ = *pin;
		   pin += d2[1];
		 }	
	       }
	     else if(pack_dim == 0) 
	       for(j=0;j < d2[1];j++) {
		 pin = (Type2 *) tmpbuf + j;
		 for(i=mystart;i < myen;i++) {
		   *pout++ = *pin;
		   pin += d2[1];
		 }	
	       }
	     mystart = myen;
	   }

	 }
	 else { // Paccking Z dim. or no packing
	   if(pack_dim == 2) 
	     out1 = out-kst*d2[0]*d2[1];
	   else 
	     out1 = out;
	   
	   pout = out1 + k*d2[0]*d2[1];
	   for(j=0;j < d2[1];j++) {
	     pin = (Type2 *) tmpbuf + j;
	     for(i=0;i < d2[0];i++) {
	       *pout++ = *pin;
	       pin += d2[1];
	     }	
	   }
	 }
     }
#endif

     /*
     if((void *) in == (void *) out && dt2 > dt1)  
       for(k=offset2[slice]/(d2[0]*d2[1]);k < (offset2[slice]+mysize2[slice])/(d2[0]*d2[1]);k++) {
	 
#ifdef MKL_BLAS
	 pout = out+k*d2[0]*d2[1];
	 blas_trans<Type2>(d2[1],d2[0],1.0,(Type2 *) tmpbuf,d2[1],pout,d2[0]);
#else
	 if(pack_dim == 0 || pack_dim == 1) {
	   mystart = start = 0;
	   sz = d2[pack_dim]/pack_procs;
	   l = d2[pack_dim] % pack_procs;
	   nl = pack_procs - l;

	   for(ipack=0;ipack<pack_procs;ipack++) {
	     if(ipack >= nl)
	       mypacksize = sz+1;
	     else
	       mypacksize = sz;
	     start = mystart * (MULT3(d2)/d2[pack_dim]);
	     pout = out+start+k*d2[0]*mypacksize;

	     myen = mystart+mypacksize;
	     if(pack_dim == 1) 
	       for(j=mystart;j < myen;j++) {
		 pin = (Type2 *) tmpbuf + j;
		 for(i=0;i < d2[0];i++) {
		   *pout++ = *pin;
		   pin += d2[1];
		 }	
	       }
	     else if(pack_dim == 0) 
	       for(j=0;j < d2[1];j++) {
		 pin = (Type2 *) tmpbuf + j;
		 for(i=mystart;i < myen;i++) {
		   *pout++ = *pin;
		   pin += d2[1];
		 }	
	       }
	     mystart = myen;
	   }

	 }
	 else { // Paccking Z dim. or no packing
	   pout = out +k*d2[0]*d2[1];
	   for(j=0;j < d2[1];j++) {
	     pin = (Type2 *) tmpbuf + j;
	     for(i=0;i < d2[0];i++) {
	       *pout++ = *pin;
	       pin += d2[1];
	     }	
	     
	   }
	 }
#endif
       }
     */
     // if(!is_empty)
     //	delete [] tmp;
#endif
}	 
	 

template <class Type1,class Type2> void transplan<Type1,Type2>::rot120in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice,int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{

  
#ifdef CUDA  

  //  Type2 *tmp;
  int imc[3] = {2,0,1};

  if(exec) {
    //     cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
     (*(exec))(plan->libplan_out[nslices],in,(Type2 *) tmpbuf);
  }
  else 
    tmpbuf = (char *) in;
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,d2);
  

#ifdef CUTENSOR
  ro_cutensor_in<Type2>((Type2 *) tmpbuf,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>((Type2 *) tmpbuf,out+k,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro120in_cu(gridDim,blockSize,(Type2 *) tmpbuf,out,d2);
#endif

    //  if(!is_empty) 
    //cudaFree(tmp);

#ifdef DEBUG
  int size=d2[0]*d2[1]*d2[2];
  Type2 *tmp1=new Type2[size];;  
  cudaMemcpy(tmp1,out, size, cudaMemcpyDeviceToHost);  
  char str[80];
  int taskid;
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  sprintf(str,"rot120in.%d",taskid);
  int imo1[] = {0,1,2};
  int din[3] = {1,1,size};
  // for(int i=0;i<3;i++)
  //  din[imc[i]] = d2[i];
  write_buf<Type2>(tmp1,str,din,imo1);
  delete [] tmp1;
#endif

#else
    {

 // no Cuda
  
      Type2 *tmp,*pin2,*pout2,*pout0,*out1;
    Type2 *pin,*pin1,*pout,*pout1,*ptran2; 
    int nb23, nb32,d0;
    int i,j,k,ii,jj,kk,i2,j2,k2,kinit,k0,lim,dmult;
    int sz,l,nl,ipack,isize,jst,jen,jsize,mystart,myen,start,mypacksize;

    tmp = (Type2 *) tmpbuf;

    if(d1[0]*d1[1] >0)
      nb32 = nb23 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
    else
      nb32 = nb23 = 1;
    if(nb32 < 1) nb32 = 1;
    if(nb23 < 1) nb23 = 1;
    
    //	tmp = new Type2[d2[2]*d2[1]*nb31];
    int sdims[3];
    
    if(deriv) {
      sdims[0] = d2[1];
      sdims[1] = d2[2];
      sdims[2] = 1;
    }
    
    int kst = offset1[slice]/(d1[0]*d1[1]);
    int ken = kst + mysize1[slice]/(d1[0]*d1[1]);
    if(pack_dim == 0) {
      sz = (ken-kst)/pack_procs;
      l = (ken-kst) % pack_procs;
      nl = pack_procs - l;
      dmult = (MULT3(d2)/d2[pack_dim]); 
      lim = sz * nl;
    }
    if(pack_dim == 1) {
      sz = d2[pack_dim]/pack_procs;
      l = d2[pack_dim] % pack_procs;
      nl = pack_procs - l;
      dmult = (ken-kst) * d2[2];
      lim = sz * nl;
    }

    if(pack_dim == 0) 
      kinit = ar3d_cnt(kst,pack_procs,sz,l,dmult);
    //    for(k=0;k <d1[2];k+=nb31) {
    for(k=kst;k < ken;k+=nb32) {
      k2 = min(k+nb32,ken);
      // if(inplace) {
      //  tmp = (Type2 *) (in+k*d1[0]*d1[1]);
      //  (*(trans_type->exec))(plan->libplan_in,in+k*d1[0]*d1[1],tmp);
      // }
      //else
      if(exec) {
	for(kk=k;kk < k2;kk++) {
	  pout = (Type2 *) tmpbuf + d2[1]*d2[2]*kk;
	  (*exec)(plan->libplan_out[slice],in+kk*d1[0]*d1[1],pout);
	  if(deriv)
	    compute_deriv_loc(pout,pout,sdims);
	}
      }
      else if((void *) in == (void *) out)
	memcpy(tmpbuf+k*d1[0]*d1[1],in+k*d1[0]*d1[1],d2[2]*d2[1]*(k2-k)*sizeof(Type2));
      else
	tmpbuf = (char *) (in + k*d1[0]*d1[1]);
      
      if(pack_dim == 1) {

	for(j=0;j < d1[1];j+= nb23) {
	  j2 = min(j+nb23,d1[1]);

	  mystart = 0;
	  for(ipack=0;ipack<pack_procs;ipack++) {
	    if(ipack < nl)
	      jsize = sz;
	    else
	      jsize = sz+1;
	    myen = mystart + jsize;
	    start = mystart * dmult;
	    
	    for(kk=k; kk < k2; kk++) {
	      pin1 =  (Type2 *) tmpbuf + j * d2[1] + kk *d2[1]*d2[2] + mystart;
	      pout1 =  out + start + kk-kst + (ken-kst)* j * jsize;
	      for(jj=j;jj < j2;jj++) {
		pin = pin1;
		pout = pout1;
		for(i=mystart;i < myen;i++) {
		  *pout = *pin++;
		  pout += (ken-kst);
		}
		pin1 += d2[1];
		pout1+= (ken-kst)*jsize;
	      }
	    }
	    mystart += jsize;
	  }
	}
      }
      else if(pack_dim == 0) {
	
	for(j=0;j < d1[1];j+= nb23) {
	  j2 = min(j+nb23,d1[1]);

	  k0 = kst;//kinit;
	  for(kk=k; kk < k2; kk++) {
	    if(k < lim)
	      isize = sz;
	    else
	      isize = sz+1;
	    pin1 =  (Type2 *) tmpbuf + j * d2[1] + kk *d2[1]*d2[2];
	    pout1 =  out + (kk-k0) + isize* j * d2[1] + k0 * dmult;
	    for(jj=j;jj < j2;jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d2[1];i++) {
		*pout = *pin++;
		pout += isize;
	      }
	      pin1 += d2[1];
	      pout1+=isize*d2[1];
	    }
	    if(kk - k0 == isize-1)
	      k0 += isize;
	  }
	}
	kinit = k0;
      }
      else { // pack_dim = 2 or no pack
	if(pack_dim == 2) {
	  out1 = out-kst;
	  d0 = ken-kst;
	}
	else {
	  out1 = out;
	  d0 = d2[0];
	}
	
	for(j=0;j < d1[1];j+= nb23) {
	  j2 = min(j+nb23,d1[1]);
	      
	  for(kk=k; kk < k2; kk++) {
	    pin1 =  (Type2 *) tmpbuf + j * d2[1] + kk *d2[1]*d2[2];
	    pout1 =  out1 + kk + d0 * j * d2[1] ;
	    for(jj=j;jj < j2;jj++) {
	      pin = pin1;
	      pout = pout1;
	      for(i=0;i < d2[1];i++) {
		*pout = *pin++;
		pout += d0;
	      }
	      pin1 += d2[1];
	      pout1+=d0*d2[1];
	    }
	  }
	}
      }    
    }

    }
//    if(!is_empty)
    // delete [] tmp;
#endif
}


inline int ar3d_cnt(int init,int pack_procs,int sz,int l,int od)
{
  int lim = sz * (pack_procs - l);
  
  if(init < lim)
    return(od*init);
  else
    return(od*(lim + ((init-lim)/(sz+1))*(sz+1) ));
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot210in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice,int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{

#ifdef CUDA  

  Type2 *tmp;
  if(exec) {
    //    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
    (*(exec))(plan->libplan_out[nslices],in,tmpbuf);
  }
  else
    tmpbuf = (char *) in;
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,d2);

  int imc[3] = {2,1,0};

#ifdef CUTENSOR
  ro_cutensor_in<Type2>((Type2 *) tmpbuf,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>((Type2 *) tmpbuf,out,imc,d2);
#else // do our own CUDA transpose
  dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
  dim3 blockSize=(TILE_DIM,TILE_DIM);
    //    ro120in_cu<Type2>(gridDim,blockSize,tmpbuf,out,d2);
  ro210in_cu(gridDim,blockSize,(Type2 *) tmpbuf,out,d2);
#endif
    //   if(!is_empty)
    //  cudaFree(tmpbuf);
 
#else 

  Type2 *pout,*pin,*pin1,*pin2,*pout1,*pout0,*out1;
  int i,ii,i2,j,k,kk,k2,nb31,nb13,j0,jinit,lim,kinit,k0,dmult,d0;
  int sz,l,nl,ipack,isize,jst,jen,jsize,mystart,myen,start,mypacksize;


  if(d1[0]*d1[1] >0)
    nb31 = cache_bl / (sizeof(Type2)*d1[0]*d1[1]);
  else
    nb31 = 1;
  if(nb31 < 1) nb31 = 1;
  if(d2[0]*d2[1] >0)
    nb13 = cache_bl / (sizeof(Type2)*d2[0]*d2[1]);
  else
    nb13 = 1;
  if(nb13 < 1) nb13 = 1;
 
  
  //  tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
	
      //	tmpbuf = new Type2[d2[2]*d2[1]*d2[0]];
  int kst = offset1[slice]/(d1[0]*d1[1]);
  int ken = kst + mysize1[slice]/(d1[0]*d1[1]);
  if(pack_dim == 0) {
    sz = (ken-kst)/pack_procs;
    l = (ken-kst) % pack_procs;
    nl = pack_procs - l;
    dmult = (MULT3(d2)/d2[pack_dim]); 
    lim = sz * nl;
    //    kinit = ar3d_cnt(kst,pack_procs,sz,l,dmult);
  }
  else   if(pack_dim == 1) {
    sz = d2[pack_dim]/pack_procs;
    l = d2[pack_dim] % pack_procs;
    nl = pack_procs - l;
    dmult = (ken-kst) * d2[2];
    lim = sz * nl;
  }
  
  int sdims[3];
  if(deriv) {
    sdims[0] = d2[2];
    sdims[1] = d2[1];
    sdims[2] = 1;
  }

  for(k=kst;k < ken;k+=nb31) {
    k2 = min(k+nb31,ken);
    
    if(exec) {
      for(kk=k;kk < k2;kk++) {
	pout = (Type2 *) tmpbuf + d2[1]*d2[2]*kk;
	(*(exec))(plan->libplan_out[slice],in+kk*d1[0]*d1[1],pout);
	if(deriv) 
	  compute_deriv_loc(pout,pout,sdims);
      }
    }
    else if((void *) in == (void *) out)
      memcpy(tmpbuf,in+k*d1[0]*d1[1],d2[2]*d2[1]*(k2-k)*sizeof(Type2));
    else
      tmpbuf = (char *) (in + k*d1[0]*d1[1]);
      
    
    if(pack_dim == 1) {
            
      for(i=0;i < d2[2];i+=nb13) {
	i2 = min(i+nb13,d2[2]);
	for(kk=k; kk < k2; kk++) {
	  pin1 =  (Type2 *) tmpbuf + kk*d2[2]*d2[1] +i;
	  pout1 =  out + (kk-kst) ;
	  j0 = 0;
	  for(j=0;j < d2[1];j++) {
	    if(j < lim)
	      jsize = sz;
	    else
	      jsize = sz+1;
	    pin = pin1;
	    pout = pout1+(ken-kst)*(j0 + i *jsize);
	    for(ii=i; ii < i2; ii++) {
	      *pout = *pin++;
	      pout += (ken-kst)*jsize;
	    }
	    pin1 += d2[2];
	    pout1 += (ken-kst);
	  }
	  if(j - j0 == jsize-1)
	    j0 += jsize;
	}
      }
    }
    else if(pack_dim == 0) {
    
      for(i=0;i < d2[2];i+=nb13) {
	i2 = min(i+nb13,d2[2]);
	k0 = kst;//kinit;
	for(kk=k; kk < k2; kk++) {
	  pin1 =  (Type2 *) tmpbuf + kk*d2[2]*d2[1] +i;
	  if(kk < lim)
	    isize = sz;
	  else
	    isize = sz+1;
	  pout1 =  out + (kk -k0 + i *d2[1]*isize) + (k0-kst) * dmult; //ar3d_cnt(kk,pack_procs,sz,l,dmult) + 
	  for(j=0;j < d2[1];j++) {
	    pin = pin1;
	    pout = pout1;
	    for(ii=i; ii < i2; ii++) {
	      *pout = *pin++;
	      pout += isize*d2[1];
	    }
	    pin1 += d2[2];
	    pout1 += isize;
	  }
	  
	  if(kk-k0 == isize-1)
	    k0 += isize;
	}
      }
      kinit = k0;
    }
    
    else { // pack_dim = 2 or no pack
	if(pack_dim == 2) {
	  out1 = out-kst;
	  d0 = ken-kst;
	}
	else {
	  out1 = out;
	  d0 = d2[0];
	}
      for(i=0;i < d2[2];i+=nb13) {
	i2 = min(i+nb13,d2[2]);
	for(kk=k; kk < k2; kk++) {
	  pin1 =  (Type2 *) tmpbuf + kk*d2[2]*d2[1] +i;
	  pout1 =  out1 + kk + i *d2[1]*d0;
	  for(j=0;j < d2[1];j++) {
	    pin = pin1;
	    pout = pout1;
	    for(ii=i; ii < i2; ii++) {
	      *pout = *pin++;
	      pout += d0*d2[1];
	    }
	    pin1 += d2[2];
	    pout1 += d0;
	  }
	}
      }
    }
  }
  //  delete [] tmpbuf;
#endif
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot201in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice,int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{

#ifdef CUDA  
  //  Type2 *tmpbuf;
  if(exec) {
    //    cudaMalloc(reinterpret_cast<void **> (&tmpbuf),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
    (*(exec))(plan->libplan_out[slice],in,tmpbuf);
  }
  else
    tmpbuf = (char *) in;
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmpbuf,tmpbuf,d2);

  int imc[3] = {1,2,0};

#ifdef CUTENSOR
  ro_cutensor_in<Type2>((Type2 *) tmpbuf,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>((Type2 *) tmpbuf,out,imc,d2);
#else // do our own CUDA transpose
  dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
  dim3 blockSize=(TILE_DIM,TILE_DIM);
  //    ro201in_cu<Type2>(gridDim,blockSize,tmpbuf,out,d2);
    //ro201in_cu(gridDim,blockSize,(Type2 *) tmpbuf,out,d2);
#endif
    //    if(!is_empty)
    //  cudaFree(tmpbuf);

#else
  Type2 *pout,*pin,*pin1,*pin2,*pout1,*out1;
  int i,j,jj,j2,k,kk,k2,nb31,nb13,j0,jinit,lim,kinit,k0,i2,ii;
  int sz,l,nl,ipack,jsize,isize,dother,dmult,dk;


  if(d2[0]*d2[1] > 0)
    nb13 = cache_bl / (sizeof(Type2)*d2[0]*d2[1]);
  else
    nb13 = 1;
  if(nb13 < 1) nb13 = 1;
  if(d2[0]>0)
    nb31 = cache_bl / (sizeof(Type2)*d2[0]*nb13);
  else
    nb31 = 1;
  if(nb31 < 1) nb31 = 1;
  
  //  tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
  
	//	tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
  int kst = offset1[slice]/(d1[0]*d1[1]);
  int ken = kst + mysize1[slice]/(d1[0]*d1[1]);
  if(pack_dim == 0)  {
    sz = d2[pack_dim]/pack_procs;
    l = d2[pack_dim] % pack_procs;
    nl = pack_procs - l;
    //    dother = d2[0]*d2[1]/d2[pack_dim];
    dmult = (ken-kst) * d2[2];
    lim = sz * nl;
  }
  else if(pack_dim == 1)  {
    sz = (ken-kst)/pack_procs;
    l = (ken-kst) % pack_procs;
    nl = pack_procs - l;
    //    dother = d2[0]*d2[1]/d2[pack_dim];
    dmult = MULT3(d2)/d2[pack_dim];
    lim = sz * nl;
  }
    
  int sdims[3];
  if(deriv) {
    sdims[0] = d2[2];
    sdims[1] = d2[0];
    sdims[2] = 1;
  }

  for(k=kst;k < ken;k+=nb31) {
    k2 = min(k+nb31,ken);
    
    if(exec) {
      for(kk=k;kk < k2;kk++) {
	pout = (Type2 *) tmpbuf + d2[0]*d2[2]*kk;
	(*(exec))(plan->libplan_out[slice],in+kk*d1[0]*d1[1],pout);
	if(deriv) 
	  compute_deriv_loc(pout,pout,sdims);
      }
    }
    else if((void *) in == (void *) out)
      memcpy(tmpbuf,in+k*d1[0]*d1[1],d2[2]*d2[0]*(k2-k)*sizeof(Type2));
    else
      tmpbuf = (char *) (in + k*d1[0]*d1[1]);
    
    
    if(pack_dim == 1) {
      kinit = ar3d_cnt(kst,pack_procs,sz,l,dmult);
      for(i=0;i < d2[2];i+= nb13) {
	i2 = min(i+nb13,d2[2]);
	for(j=0;j <d1[1];j++) {
	  k0 = kst; //kinit;
	  pin1 = (Type2 *) tmpbuf + i + j*d2[2]+k*d2[0]*d2[2];
	  pout1 = out + j;
	  for(kk=k; kk < k2; kk++) {
	    if(kk < lim)
	      jsize = sz;
	    else
	      jsize = sz+1;
	    pin = pin1;
	    pout =  pout1 + d2[0]*(kk -k0) + d2[0]*jsize*i + k0 * dmult;
	    for(ii=i;ii < i2;ii++) {
	      *pout = *pin++;
	      pout += d2[0]*jsize;
	    }
	    pin1 += d2[2]*d2[0];
	    pout1+= d2[0];
	    if(kk - k0 == jsize-1)
	      k0 += jsize;
	  }
	}
      }
    }
    else if(pack_dim == 0) {
      
      for(i=0;i < d2[2];i+= nb13) {
	i2 = min(i+nb13,d2[2]);
	j0 = 0;
	for(j=0;j <d1[1];j++) {
	  if(j < lim) 
	    isize = sz;
	  else
	    isize = sz+1;
	  pin1 = (Type2 *) tmpbuf + i + j*d2[2]+k*d2[0]*d2[2];
	  pout1 =  out + j -j0 +isize*(k-kst + i *(ken-kst)) + j0 * dmult;//ar3d_cnt(jj,pack_procs,sz,l,dmult);;
	  for(kk=k; kk < k2; kk++) {
	    pin = pin1;
	    pout = pout1;
	    for(ii=i;ii < i2;ii++) {
	      *pout = *pin++;
	      pout += isize*(ken-kst);
	    }
	    pin1 += d2[2]*d2[0];
	    pout1+= isize;
	  }
	  if(j-j0 == isize-1)
	    j0 += isize;
	}
      }
    }
    
    else {
      	if(pack_dim == 2) {
	  out1 = out-kst*d2[0];
	  dk = ken-kst;
	}
	else {
	  out1 = out;
	  dk = d2[1];
	}

      for(i=0;i < d2[2];i+= nb13) {
	i2 = min(i+nb13,d2[2]);
	for(j=0;j <d1[1];j++) {
	  pin1 = (Type2 *) tmpbuf + i + j*d2[2]+k*d2[0]*d2[2];
	  pout1 =  out1 + j + d2[0]*(k + i * dk);
	  for(kk=k; kk < k2; kk++) {
	    pin = pin1;
	    pout = pout1;
	    for(ii=i;ii < i2;ii++) {
	      *pout = *pin++;
	      pout += d2[0]*dk;
	    }
	    pin1 += d2[2]*d2[0];
	    pout1+= d2[0];
	  }
	}
      }
    }
  }
  
  //  delete [] tmpbuf;
#endif
}

int find_disp(int dim,int slice,int nslices);

int find_disp(int dim,int slice,int nslices)
{
  int d=0;
  int i,l,s;
  
  s = dim/nslices;
  l = dim % nslices;
  for(i=0;i< min(l,slice);i++)
    d += s+1;
  for(;i<slice;i++)
    d += s;
  return(d);
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot021_op_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice,int nslices,bool deriv, int pack_dim, int pack_procs)
{

  if((void *) in == (void *) out) {
    printf("Error in rot021_op_slice: expected non-overlapping input and output\n");
    return;
  }

  Type1 *pin_t1;
  Type2 *pout,*pin,*pin1,*pin2,*pout1,*pout2,*out1;
  int i,j,k,nb32,nb23,kk,k2,jj,j2,dk;
  int sz,l,nl,ipack,jsize,isize,mysz,start,myst;

  if(d1[0]*d1[1] > 0)
    nb32 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else
    nb32 = 1;
  if(nb32 < 1) nb32 = 1;
  if(d1[0]*d1[1] > 0)
    nb23 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else nb23 = 1;
  if(nb23 < 1) nb23 = 1;
  

  int kst = offset1[slice]/(d1[0]*d1[1]);
  int ken = kst + mysize1[slice]/(d1[0]*d1[1]);
  int sdims[3] = {d2[2],1,1};

  if(pack_dim == 0) {
    sz = d2[0]/pack_procs;
    l = d2[0] % pack_procs;
    nl = pack_procs - l;
    Type2 *tmpbuf;
    if(exec)
      tmpbuf=new Type2[d2[0]];

    for(k=kst;k < ken;k+=nb32) {
      k2 = min(k+nb32,ken);
      for(j=0;j < d1[1];j+=nb23) {
	j2 = min(j+nb23,d1[1]);
	for(kk=k; kk < k2; kk++) {
	  for(jj=j; jj < j2; jj++) {
	    pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
	    if(exec) {
	      pout = tmpbuf;
	      (*(exec))(plan->libplan_out[slice],pin_t1,pout);
	      if(deriv) 
		compute_deriv_loc(pout,pout,sdims);
	    }
	    else
	      pout = (Type2 *) pin_t1;

	    start=0;myst=0;
	    for(ipack = 0;ipack < pack_procs;ipack++) {
	      if(ipack < nl)
		mysz = sz;
	      else
		mysz = sz+1;
	      pout2 = out + start + mysz * (kk + jj * (ken-kst));
	      start += mysz*(ken-kst)*d2[2];
	      memcpy(pout2,pout+myst*sizeof(Type2),mysz*sizeof(Type2));
	      myst += mysz;
	    }
	  }
	}
      }
    }
    if(exec)
      delete [] tmpbuf;
  }
  else if(pack_dim == 1) {

    sz = (ken-kst)/pack_procs;
    l = (ken-kst) % pack_procs;
    nl = pack_procs - l;
    int dmult = (MULT3(d2)/d2[pack_dim]); 

    for(k=kst;k < ken;k+=nb32) {
      k2 = min(k+nb32,ken);
      for(j=0;j < d1[1];j+=nb23) {
	j2 = min(j+nb23,d1[1]);
	for(kk=k; kk < k2; kk++) {
	  ipack = kk/pack_procs;
	  if(ipack < nl)
	    jsize = sz;
	  else
	    jsize = sz+1;
	  pout = out + ar3d_cnt(kk,pack_procs,sz,l,dmult);
	  for(jj=j; jj < j2; jj++) {
	    pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
	    pout2 = pout + kk * d2[0] + jj * d2[0]*jsize;
	    if(exec) {
	      (*(exec))(plan->libplan_out[slice],pin_t1,pout2);
	      if(deriv) 
		compute_deriv_loc(pout2,pout2,sdims);
	    }
	    else
	      memcpy(pout2,pin_t1,d2[0]*sizeof(Type2));
	  }
	}
      }
    }
  }
  else {
    if(pack_dim == 2) {
      out1 = out-kst*d2[0];
      dk = ken-kst;
    }
    else {
      out1 = out;
      dk = d2[1];
    }
    
    for(k=kst;k < ken;k+=nb32) {
      k2 = min(k+nb32,ken);
      for(j=0;j < d1[1];j+=nb23) {
	j2 = min(j+nb23,d1[1]);
	for(kk=k; kk < k2; kk++) {
	  for(jj=j; jj < j2; jj++) {
	    pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
	    pout2 = out1 + kk * d2[0] + jj * d2[0]*dk;
	    if(exec) {
	      (*(exec))(plan->libplan_out[slice],pin_t1,pout2);
	      if(deriv) 
		compute_deriv_loc(pout2,pout2,sdims);
	    }
	    else
	      memcpy(pout2,pin_t1,d2[0]*sizeof(Type2));
	  }
	}
      }
    }
  }
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot021_ip(Type1 *in,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,bool deriv, char *tmpbuf)
{
    /*
    if(dt2 <= dt1) {
#ifdef CUDA
      cudaCheckErrors(cudaMalloc(&tmpbuf,d2[0]*sizeof(Type2)));
#else
      tmpbuf = new Type2[d2[0]];
    //	  tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
      int sdims[3] = {d2[0],1,1};
#endif    

      float d0 = d2[0]*dt1/dt2;
      int ind;
      for(k=kst;k < ken;k+=nb32) {
	k2 = min(k+nb32,ken);
	for(j=0;j < d2[2];j+=nb23) {
	  j2 = min(j+nb23,d2[2]);
	  for(kk=k; kk < k2; kk++) {
	    for(jj=j; jj < j2; jj++) {
	      (*(exec)(plan->libplan_out,in+d1[0]*(jj+kk*d1[1]),tmpbuf);
	      ind = (int) d0*(kk+jj*d2[1]);
	      pout = (Type2 *) in+d1[0]*(jj+kk*d1[1]);
	      (*(exec))(plan->libplan_out,in + ind,pout);
#ifdef CUDA
	      checkCudaErrors(cudaMemcpy(in + ind,tmpbuf,d2[0]*sizeof(Type2),cudaMemcpyDeviceToDevice));
#else
	      if(deriv) { 
		compute_deriv_loc(tmpbuf,tmpbuf,sdims);
		compute_deriv_loc(pout,pout,sdims);
	      }
	      memcpy(in + ind,tmpbuf,d2[0]*sizeof(Type2));
#endif
	    }
	  }
	}
      }
#ifdef CUDA
      cudaCheckErrors(cudaFree(tmpbuf));
#else
      delete [] tmpbuf;
#endif
    }
    else {
    */
  Type2 *pout;
  int i,j,k,nb32,nb23,kk,k2,jj,j2;
  if(d1[0]*d1[1] > 0)
    nb32 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else
    nb32 = 1;
  if(nb32 < 1) nb32 = 1;
  if(d1[0]*d1[1] > 0)
    nb23 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else nb23 = 1;
  if(nb23 < 1) nb23 = 1;

  /*
#ifdef CUDA
  checkCudaErrors(cudaMalloc(&tmpbuf,d2[0]*d2[1]*d2[2]*sizeof(Type2)));
#else
  tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
  //	  tmpbuf = new Type2[d2[0]*d2[1]*d2[2]];
#endif    
  */
  int sdims[3] = {d2[0],1,1};
  
  for(k=0;k < d2[1];k+=nb32) {
    k2 = min(k+nb32,d2[1]);
    for(j=0;j < d2[2];j+=nb23) {
      j2 = min(j+nb23,d2[2]);
      for(kk=k; kk < k2; kk++) {
	for(jj=j; jj < j2; jj++) {
	  pout = (Type2 *) tmpbuf+d2[0]*(kk+jj*d2[1]);
	  if(exec) {
	    (*(exec))(plan->libplan_out[0],in+d1[0]*(jj+kk*d1[1]),pout);
	    if(deriv) 
	      compute_deriv_loc(pout,pout,sdims);
	  }
	  else
	    memcpy(pout,in+d1[0]*(jj+kk*d1[1]),d2[0]*sizeof(Type2));
	}
      }
    }
  }
#ifdef CUDA
  checkCudaErrors(cudaMemcpy(in,tmpbuf,d2[0]*d2[1]*d2[2]*sizeof(Type2),cudaMemcpyDeviceToDevice));
  //  checkCudaErrors(cudaFree(tmpbuf));
#else
  memcpy(in,tmpbuf,d2[0]*d2[1]*d2[2]*sizeof(Type2));
  //  delete [] tmpbuf;
#endif
}
 

template <class Type1,class Type2> void transplan<Type1,Type2>::rot102out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int slice,int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{



#ifdef CUDA  
  /*
      Type1 *tmpbuf;
    if(!is_empty)
      cudaMalloc(reinterpret_cast<void **> (&tmpbuf),sizeof(Type1)*mysize1[slice]);
    else
      tmpbuf = (Type1 *) out;
  */
  int imc[3] = {1,0,2};
  d1[2] = mysize1[slice]/(d1[0]*d1[1]);
  int sdims[3];

#ifdef CUTENSOR
  //    ro_cutensor_out<Type1>(in+offset1[slice],(Type1 *) (out+offset2[slice]),imc,d1,slice,nslices,streams[slice]);
  ro_cutensor_out<Type1>(in+offset1[slice],((Type1 *) tmpbuf)+offset1[slice],imc,d1,slice,nslices,streams[slice]);
#elif defined CUTT
    ro_cutt_out<Type1>(in+offset1[slice],out+offset2[slice],imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro102out_cu(gridDim,blockSize,in+offset1[slice],out+offset2[slice],d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    if(exec)
      (*(exec))(plan->libplan_out[slice],((Type1 *) tmpbuf) + offset1[slice],out+offset2[slice]);
    if(deriv) {
      /*
      int l=d2[2] % nslices;
      int mysz;
      if(slice < l)
	mysz = d2[2]/nslices+1;
      else
	mysz = d2[2]/nslices;
      int sdims[3]={d2[0],d2[1],mysz};
      */
      sdims[0] = d2[0];
      sdims[1] = d2[1];
      sdims[2] = d1[2];
      //      compute_deriv_loc_cu<<<>>>(out+offset2[slice],out+offset2[slice],sdims);
    }
    //    if(!is_empty)
    // cudaFree(tmpbuf);
#else

    //    if((void *) in == (void *) out)
    //  printf("Error in rot102out: expected different, nonoverlapping input and output arrays\n");


  Type1 *pin,*pout,*pout1;
  Type2 *pout2,*pin2,*pout3;
  int i,j,k;
  int ipack,mypacksize,mystart,start,myen,l,nl,sz,dmult,dother;

  int sdims[3];
  if(deriv) {
    sdims[0] = d2[0];
    sdims[1] = d2[1];
    sdims[2] = 1;
  }
  int kst = offset1[slice]/(d1[0]*d1[1]);
  int ken = kst + mysize1[slice]/(d1[0]*d1[1]);

  if(pack_dim == 0 || pack_dim == 1)  {
    sz = d2[pack_dim]/pack_procs;
    l = d2[pack_dim] % pack_procs;
    nl = pack_procs - l;
    dother = d2[1-pack_dim];
    dmult = (ken-kst) * dother;
  }

  pout1 = (Type1 *) tmpbuf;

  for(k=kst;k < ken;k++) {
    pin = in + k*d1[0]*d1[1];
    if(is_empty && pack_dim != 0 && pack_dim != 1) 
      //      tmpbuf = new Type1[d1[0]*d1[1]];
      //else
      tmpbuf = (char *) (out + d2[0]*d2[1]*k);
#ifdef MKL_BLAS
    blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmpbuf,d1[1]);
#else
    for(j=0;j < d1[1];j++) {
      pout = pout1 +j;
      for(i=0;i < d1[0];i++) {
	*pout = *pin++;
	pout += d1[1];
      }	
    }
#endif
      
    if(pack_dim == 0 || pack_dim == 1) 
      pout2 = (Type2*) tmpbuf + d2[0]*d2[1];
    else if(pack_dim == 2)
      pout2 = out + d2[0]*d2[1]*(k-kst);
    else
      pout2 = out + d2[0]*d2[1]*k;
    
    if(exec)
      (*(exec))(plan->libplan_out[slice],pout1,pout2);
    else
      pout2 = (Type2 *) tmpbuf;
    if(deriv)
      compute_deriv_loc(pout2,pout2,sdims);
    
    if(pack_dim == 0 || pack_dim == 1) {
      mystart  = 0;
      
      for(ipack=0;ipack<pack_procs;ipack++) {
	if(ipack >= nl)
	  mypacksize = sz+1;
	else
	  mypacksize = sz;
	start = mystart * dmult;
	pout3 = out+start+k*dother*mypacksize;
	myen = mystart+mypacksize;
	if(pack_dim == 1) {
	  pin2 = pout2+d2[0]*mystart; 
	  memcpy(pout3,pin2,d2[0]*mypacksize*sizeof(Type2));
	}
	else 
	  for(j=0;j < d2[1];j++)  {
	    pin2 = pout2+mystart+j*d2[0]; 
	    memcpy(pout3,pin2,mypacksize*sizeof(Type2));
	    pout3 += mypacksize;
	  }
	mystart = myen;
      }
    }
  }
  
  //  if(!is_empty)
  //  delete [] tmpbuf;
  
#endif
}
 

// Assume IN != OUT
template <class Type1,class Type2> void transplan<Type1,Type2>::rot120out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice, int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{

#ifdef CUDA  
  //    Type1 *tmpbuf;
    if(is_empty)
      //  cudaMalloc(reinterpret_cast<void **> (&tmpbuf),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
      //else
      tmpbuf = (char *) out;
    int imc[3] = {2,0,1}; //{1,2,0};


#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,(Type1 *) tmpbuf,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,(Type1 *) tmpbuf,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
  //  ro120in_cu<Type2>(gridDim,blockSize,tmpbuf,out,d2);
    ro120out_cu(gridDim,blockSize,in,(Type1 *) tmpbuf,d1);
#endif

    if(exec) {    
      (*(exec))(plan->libplan_out[nslices],(Type1 *) tmpbuf,out);
    //     if(deriv)
    //  compute_deriv_loc_cu<<<>>>(out,out,d2);
   //    cudaDeviceSynchronize();
      //      cudaFree(tmpbuf);
    }
    
#else
 // no Cuda

      //    if((void *) in == (void *) out)
      //printf("Error in rot120out: expected different, nonoverlapping input and output arrays\n");

    Type1 *pin,*pin1,*pout,*pout1;
    Type2 *pout2,*pout3,*pin2,*pout4;
    int i,j,k,jj,j2,kk,k2,nb23,nb32,dmult,dother;
    int sz,l,nl,ipack,isize,jsize,mystart,myen,start,mypacksize;
    
    if(d1[0]*d1[1] > 0)
      nb32 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
    else
      nb32 = 1;
    if(nb32 < 1)  nb32= 1;
    if(d1[0]*d1[2] > 0)
      nb23 = cache_bl / (sizeof(Type1)*d1[0]*d1[2]);
    else
      nb23 = 1;
    if(nb23 < 1)  nb23= 1;
    
    int sdims[3] = {d2[0],d2[1],1};
    bool alg_inplace;
    if((void *) in == (void *) out)
    alg_inplace = true;
    else
      alg_inplace = false;
    
    //	tmpbuf = new Type1[d1[0]*d1[1]*d1[2]]; // tmpbuf[k][i][j] = in[i][j][k]
    
    int jst = offset2[slice]/(d2[0]*d2[1]);
    int jen = jst + mysize2[slice]/(d2[0]*d2[1]);
    if(pack_dim == 0 || pack_dim == 1) {
      mystart  = 0;
      sz = d2[pack_dim]/pack_procs;
      l = d2[pack_dim] % pack_procs;
      nl = pack_procs - l;
      dother = d2[1-pack_dim];
      dmult = (jen-jst)*dother;
    }
    
    Type1 *pout0 = (Type1 *) tmpbuf;
    Type2 *pen = (Type2 *) (tmpbuf + d1[0]*d1[0]*d1[2]*sizeof(Type1)); 
    
    for(j=jst;j < jen;j+=nb23) {
      j2 = min(j+nb23,d2[2]);
      
      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	
	for(jj=j;jj<j2;jj++) {
	  
	  pin1 = in + d1[0]*(jj + k*d1[1]);
	  pout1 = pout0 + jj*d1[2]*d1[0] +k;
	  for(kk=k; kk < k2; kk++) {	
	    pin = pin1;
	    pout = pout1;
	    for(i=0;i < d1[0];i++) {
	      *pout = *pin++;
	      pout += d1[2];
	    }
	    pin1 += d1[0]*d1[1];
	    pout1++;
	  }
	}
      }
      
      if(!alg_inplace) {
	
	for(jj=j;jj<j2;jj++) {
	  
	  pin = pout0 + jj*d1[2]*d1[0];
	  if(pack_dim == 0 || pack_dim == 1) 
	    pout2 = pen;
	  else if(pack_dim == 2)
	    pout2 = out + (jj-jst)*d2[0]*d2[1];
	  else
	    pout2 = out + jj*d2[0]*d2[1];
	  
	  if(exec)	    
	    (*(exec))(plan->libplan_out[slice],pin,pout2);
	  else
	    memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
	  if(deriv) 
	    compute_deriv_loc(pout2,pout2,sdims);
	  
	  if(pack_dim == 1) {
	    
	    mystart = 0;
	    for(ipack=0;ipack<pack_procs;ipack++) {
	      if(ipack >= nl)
		mypacksize = sz+1;
	      else
		mypacksize = sz;
	      start = mystart * dmult;
	      pout3 = out+start+(jj-jst)*dother*mypacksize;
	      myen = mystart+mypacksize;
	      pin2 = pout2+d2[0]*mystart;
	      memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	    }
	    mystart = myen;
	  }
	  else if(pack_dim == 0) {
	    for(i=0;i < d2[1];i++)  {
	      mystart = 0;
	      for(ipack=0;ipack<pack_procs;ipack++) {
		if(ipack >= nl)
		  mypacksize = sz+1;
		else
		  mypacksize = sz;
		start = mystart * dmult;
		pout3 = out +(jj-jst)*dother*mypacksize +start+i*mypacksize;
		myen = mystart+mypacksize;
		//     pin2 = pout2+mystart+i*d2[0]; 
		pin2 = pout2+mystart+i*d2[0]; 
		memcpy(pout3,pin2,mypacksize*sizeof(Type2));
		mystart = myen;
	      }
	    }
	  }
	}      
      }
    }
    
    if(alg_inplace) 
      for(j=jst;j<jen;j++) {
	pin = pout0 + j*d1[2]*d1[0];
	if(pack_dim == 0 || pack_dim == 1) 
	  pout2 = pen;
	else if(pack_dim == 2)
	  pout2 = out + (j-jst)*d2[0]*d2[1];
	else
	  pout2 = out + j*d2[0]*d2[1];
	if(exec)	    
	  (*(exec))(plan->libplan_out[slice],pin,pout2);
	else
	  memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
	if(deriv) 
	  compute_deriv_loc(pout2,pout2,sdims);
	
	if(pack_dim == 1) {
	  
	  mystart = 0;
	  for(ipack=0;ipack<pack_procs;ipack++) {
	    if(ipack >= nl)
	      mypacksize = sz+1;
	    else
	      mypacksize = sz;
	    start = mystart * dmult;
	    pout3 = out+start+(j-jst)*dother*mypacksize;
	    myen = mystart+mypacksize;
	    pin2 = pout2+d2[0]*mystart; 
	    memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	  }
	  mystart = myen;
	}
	else if(pack_dim == 0) 
	  for(i=0;i < d2[1];i++)  {
	    mystart = 0;
	    for(ipack=0;ipack<pack_procs;ipack++) {
	      if(ipack >= nl)
		mypacksize = sz+1;
	      else
		mypacksize = sz;
	      start = mystart * dmult;
	      pout3 = out+start+(j-jst)*dother*mypacksize+i*mypacksize;
	      myen = mystart+mypacksize;
	      pin2 = pout2+mystart+i*d2[0]; 
	      memcpy(pout3,pin2,mypacksize*sizeof(Type2));
	      mystart = myen;
	    }
	  }
      }
    
#endif
}




template <class Type1,class Type2> void transplan<Type1,Type2>::rot210out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice, int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{


#ifdef CUDA  
  //    Type1 *tmpbuf;
    if(is_empty)
      //  cudaMalloc(reinterpret_cast<void **> (&tmpbuf),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
      //else
      tmpbuf = (char *) out;
    int imc[3] = {2,1,0};
#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,(Type1 *) tmpbuf,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,(Type1 *) tmpbuf,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro210out_cu(gridDim,blockSize,in,(Type1 *) tmpbuf,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    if(exec) {
      (*(exec))(plan->libplan_out[nslices],(Type1 *) tmpbuf,out);
      //      cudaFree(tmpbuf);
    }
    
#else
 // no Cuda

    //    if((void *) in == (void *) out)
    //  printf("Error in rot210out: expected different, nonoverlapping input and output arrays\n");

    Type1 *pin,*pin1,*pout,*pout1,*tmp;
  Type2 *pout2,*pout3,*pin2;
  int i,j,k,ii,i2,kk,k2,nb13,nb31,dmult,dother;
  int sz,l,nl,ipack,isize,jst,jen,jsize,mystart,myen,start,mypacksize;
  
  if(d1[0]*d1[1] >0)
    nb31 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else nb31 = 1;
  if(nb31 < 1) nb31 = 1;
  if(d2[0]*d2[1] >0)
    nb13 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else nb13 = 1;
  if(nb13 < 1) nb13 = 1;
	
  //  tmpbuf = new Type1[d1[0]*d1[1]*d1[2]];
	
  int sdims[3] = {d2[0],d2[1],1};

  int kst = offset2[slice]/(d2[0]*d2[1]);
  int ken = kst + mysize2[slice]/(d2[0]*d2[1]);
  if(pack_dim == 0 || pack_dim == 1)  {
    sz = d2[pack_dim]/pack_procs;
    l = d2[pack_dim] % pack_procs;
    nl = pack_procs - l;
    dother = d2[1-pack_dim];
    dmult = (ken-kst) * dother;
    tmp = (Type1 *) tmpbuf + d1[2]*d1[1]*nb13;
  }

  Type1 *pout0 = (Type1 *) tmpbuf;
  Type2 *pen = (Type2 *) (tmpbuf + d1[0]*d1[1]*d1[2]*sizeof(Type1)); 
  bool alg_inplace;
  
  if((void *) in == (void *) out)
    alg_inplace = true;
  else
    alg_inplace = false;
  
  for(k=kst;k < ken;k+=nb13) {
    k2 = min(k+nb13,ken);
    for(i=0;i < d1[2];i+=nb31) {
      i2 = min(i+nb31,d1[2]);
      for(kk=k; kk < k2; kk++) {
	pin1 = in + (i*d1[0]*d1[1] +kk);
	pout1 = pout0 + i + kk * d1[2]*d1[1];
	for(j=0;j < d1[1];j++) {
	  pin = pin1;
	  pout = pout1;
	  for(ii=i; ii < i2; ii++) {
	    *pout++ = *pin;
	    pin += d1[0]*d1[1];
	  }
	  pin1 += d1[0];
	  pout1 += d1[2];
	}
      }
    }

    if(!alg_inplace)
    
    for(kk=k; kk < k2; kk++) {
      //    for(kk=0; kk < d1[0]; kk++) {
      if(pack_dim == 0 || pack_dim == 1) 
	pout2 = pen;
      else if(pack_dim == 2)
	pout2 = out + (kk-kst)*d2[0]*d2[1];
      else
	pout2 = out + kk*d2[0]*d2[1];
      pin = pout0 + kk*d1[2]*d1[1];
      
      if(exec)
	(*(exec))(plan->libplan_out[slice],pin,pout2);
      else
	memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
      if(deriv)
	compute_deriv_loc(pout2,pout2,sdims);
          
      if(pack_dim == 0 || pack_dim == 1) {
	
	mystart = 0;
	for(ipack=0;ipack<pack_procs;ipack++) {
	  if(ipack >= nl)
	    mypacksize = sz+1;
	  else
	    mypacksize = sz;
	  start = mystart * dmult;
	  pout3 = out+start+(kk-kst)*dother*mypacksize;
	  myen = mystart+mypacksize;
	  if(pack_dim == 1) {
	    pin2 = pout2+d2[0]*mystart;
	    memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	  }
	  else 
	    for(j=0;j < d2[1];j++)  {
	      pin2 = pout2+mystart+j*d2[0]; 
	      memcpy(pout3,pin2,mypacksize*sizeof(Type2));
	      pout3 += mypacksize;
	    }
	  mystart = myen;
	}
	
      }
    }
  }


  if(alg_inplace)
    
    for(kk=kst; kk < ken; kk++) {
      //    for(kk=0; kk < d1[0]; kk++) {
      if(pack_dim == 0 || pack_dim == 1) 
	pout2 = pen;
      else if(pack_dim == 2)
	pout2 = out + (kk-kst)*d2[0]*d2[1];
      else
	pout2 = out + kk*d2[0]*d2[1];

      pin = pout0 + kk*d1[2]*d1[1];
      
      if(exec)
	(*(exec))(plan->libplan_out[slice],pin,pout2);
      else
	memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
      if(deriv)
	compute_deriv_loc(pout2,pout2,sdims);
          
      if(pack_dim == 0 || pack_dim == 1) {
	
	mystart = 0;
	for(ipack=0;ipack<pack_procs;ipack++) {
	  if(ipack >= nl)
	    mypacksize = sz+1;
	  else
	    mypacksize = sz;
	  start = mystart * dmult;
	  pout3 = out+start+(kk-kst)*dother*mypacksize;
	  myen = mystart+mypacksize;
	  if(pack_dim == 1) {
	    pin2 = pout2+d2[0]*mystart;
	    memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	  }
	  else 
	    for(j=0;j < d2[1];j++)  {
	      pin2 = pout2+mystart+j*d2[0]; 
	      memcpy(pout3,pin2,mypacksize*sizeof(Type2));
	      pout3 += mypacksize;
	    }
	  mystart = myen;
	}
	
      }
    }


#endif
}
	

template <class Type1,class Type2> void transplan<Type1,Type2>::rot201out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice, int nslices,bool deriv,char *tmpbuf,int pack_dim,int pack_procs)
{

#ifdef CUDA  
  //    Type1 *tmpbuf;
    if(is_empty)
      //  cudaMalloc(reinterpret_cast<void **> (&tmpbuf),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
      // else
      tmpbuf = (char *) out;

    int imc[3] = {1,2,0};
#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,(Type1 *) tmpbuf,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,(Type1 *) tmpbuf,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro201out_cu(gridDim,blockSize,in,(Type1 *) tmpbuf,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    if(exec) {
      (*(exec))(plan->libplan_out[nslices],(Type1 *) tmpbuf,out);
      //      cudaFree(tmpbuf);
    }
   
#else

// no Cuda
    bool alg_inplace;
  if((void *) in == (void *) out)
    alg_inplace = true;
  else
    alg_inplace = false;

  Type1 *pin,*pin1,*pout,*pout1,*tmp;
  Type2 *pout2,*pout3,*pin2;
  int i,j,k,jj,j2,kk,k2,nb13,nb31,i2,ii,dmult,dother;
  int sz,l,nl,ipack,isize,jst,jen,jsize,mystart,myen,start,mypacksize;
  //std::unique_ptr<Type1> myptr1;

  
  if(d1[0]*d1[1] >0)
    nb31 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  if(nb31 < 1) nb31 = 1;
  if(d1[2]*d1[1] >0)
    nb13 = cache_bl / (sizeof(Type1)*d1[2]*d1[1]);
  if(nb13 < 1) nb13 = 1;
  
  int sdims[3] = {d2[0],d2[1],1};
  int ist = offset2[slice]/(d2[0]*d2[1]);
  int ien = ist + mysize2[slice]/(d2[0]*d2[1]);
  if(pack_dim == 0 || pack_dim == 1)  {
    sz = d2[pack_dim]/pack_procs;
    l = d2[pack_dim] % pack_procs;
    nl = pack_procs - l;
    dother = d2[1-pack_dim];
    dmult = (ien-ist) * dother;
  }

  //  tmpbuf = new Type1[d1[0]*d1[1]*d1[2]];
  Type1 *pout0 = (Type1 *) tmpbuf;
  Type2 *pen = (Type2 *) (tmpbuf + d1[0]*d1[1]*d1[2]*sizeof(Type1)); 
  
  for(i=ist;i < ien;i+=nb13) {
    i2 = min(i+nb13,ien);

    for(k=0;k < d1[2];k+=nb31) {
      k2 = min(k+nb31,d1[2]);

      for(j=0;j <d1[1];j++) {
	pin1 = in + (j*d1[0] +k*d1[0]*d1[1])+i;
	pout1 = pout0 + (j +k*d1[1]) + i*d1[2]*d1[1];
	for(kk=k; kk < k2; kk++) {
	  pin = pin1;
	  pout = pout1;
	  for(ii=i;ii<i2;ii++) {
	    *pout = *pin++;
	    pout += d1[2]*d1[1];
	  }
	  pin1 += d1[0]*d1[1];
	  pout1 += d1[1];
	}
      }
    }
    if(!alg_inplace)
      for(ii=i;ii<i2;ii++) {
	
	if(pack_dim == 0 || pack_dim == 1) 
	  pout2 = pen;
	else if(pack_dim == 2)
	  pout2 = out + (ii-ist)*d2[0]*d2[1];
	else
	  pout2 = out + ii*d2[0]*d2[1];
	pin = pout0 + d1[1]*d1[2]*ii;
	if(exec)
	  (*(exec))(plan->libplan_out[slice],pin,pout2);
	else
	  memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
	if(deriv)
	  compute_deriv_loc(pout2,pout2,sdims);
	
	mystart = 0;
	if(pack_dim == 0 || pack_dim == 1) {
	  
	  for(ipack=0;ipack<pack_procs;ipack++) {
	    if(ipack >= nl)
	      mypacksize = sz+1;
	    else
	      mypacksize = sz;
	    start = mystart * dmult;
	    pout3 = out+start+(ii-ist)*dother*mypacksize;
	    myen = mystart+mypacksize;
	    if(pack_dim == 1) {
	      pin2 = pout2+d2[0]*mystart; 
	      memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	    }
	    else 
	      for(k=0;k < d2[1];k++)  {
		pin2 = pout2+mystart+k*d2[0]; 
		memcpy(pout3,pin2,mypacksize*sizeof(Type2));
		pout3 += mypacksize;
	      }
	    mystart = myen;
	  }
	}
      }
    
  }


  if(alg_inplace)
    for(ii=ist;ii<ien;ii++) {
      
      if(pack_dim == 0 || pack_dim == 1) 
	pout2 = pen;
      else if(pack_dim == 2)
	pout2 = out + (ii-ist)*d2[0]*d2[1];
      else
	pout2 = out + ii*d2[0]*d2[1];
      pin = pout0+d1[1]*d1[2]*ii;
      if(exec)
	(*(exec))(plan->libplan_out[slice],pin,pout2);
      else
	memcpy(pout2,pin,d2[0]*d2[1]*sizeof(Type2));
      if(deriv)
	compute_deriv_loc(pout2,pout2,sdims);
      
      mystart = 0;
      if(pack_dim == 0 || pack_dim == 1) {
	
	for(ipack=0;ipack<pack_procs;ipack++) {
	  if(ipack >= nl)
	    mypacksize = sz+1;
	  else
	    mypacksize = sz;
	  start = mystart * dmult;
	  pout3 = out+start+(ii-ist)*dother*mypacksize;
	  myen = mystart+mypacksize;
	  if(pack_dim == 1) {
	    pin2 = pout2+d2[0]*mystart; 
	    memcpy(pout3,pin2,mypacksize*d2[0]*sizeof(Type2));
	  }
	  else 
	    for(k=0;k < d2[1];k++)  {
	      pin2 = pout2+mystart+k*d2[0]; 
	      memcpy(pout3,pin2,mypacksize*sizeof(Type2));
	      pout3 += mypacksize;
	    }
	  mystart = myen;
	}
      }
    }
  
#endif
}






