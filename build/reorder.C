
/*
#ifdef CUDA 
#include "reorder.cu"
#endif
*/

//#include "p3dfft.h"

#ifdef CUDA
#include "reorder_kernel.cu"
#endif

//namespace p3dfft {

#ifdef CUBLAS
cublasHandle_t cublas_handle;
#endif 

#ifdef CUDA
#ifdef CUTENSOR
//template <class Type1,class Type2> void transplan<Type1,Type2>::

template <class T> void ro_cutensor_in(T *in,T *out,int imc[3],int dout[3],int slice,int nslice,cudaStream_t stream)
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

 extB[imc[2]] /= nslice;
 extA[2] /= nslice;

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

template <class T> void ro_cutensor_out(T *in,T *out,int imc[3],int din[3],int slice,int nslice,cudaStream_t stream)
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

 extB[2] /= nslice;
 extA[imc[2]] /= nslice;

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


template <class Type1,class Type2> void transplan<Type1,Type2>::rot102in_slice(Type1 *in,Type2 *out,bool inplace,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int slice,int nslice,bool deriv, bool OW)
{

#if defined CUDA && !defined CUBLAS  

    Type2 *tmp;

    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*mysize2[slice]);
    (*(exec))(plan->libplan_out[slice],in+offset1[slice],tmp);
    if(deriv) {
      int sdims[3] = {d2[1],d2[0],mysize2[slice]/(d2[0]*d2[1])};
      //      compute_deriv_loc_cu<<<>>>(tmp,tmp,sdims);
    }

    int imc[3] = {1,0,2};

#ifdef CUTENSOR
    ro_cutensor_in<Type2>(tmp,out+offset2[slice],imc,d2,slice,nslice,streams[slice]);
#elif defined CUTT
	ro_cutt_in<Type2>(tmp,out,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro102in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    rot102in_cu<<<gridDim,blockSize>>>(in,out,inplace,d1,d2,exec,plan);
#endif

    cudaFree(tmp);
#else
   {
     Type2 *tmp,*pout,*pin;

     /*
     if(OW && dt2 <= dt1) 
       inplace=true;
     else {
              if((void *) in == (void *) out && dt2 > dt1) {
	 tmp = new Type2[d2[0]*d2[1]*d2[2]];
	 for(k=0;k <d1[2];k++) 
	   (*(trans_type->exec))(plan->libplan_out,in+k*d1[0]*d1[1],tmp+k*d2[0]*d2[1]);
       }
       else
       */
     tmp = new Type2[d2[0]*d2[1]];
     
     for(k=0;k <d1[2];k++) {
       /*
       if(inplace) {
	 tmp = (Type2 *) in + k*d2[0]*d2[1];
	 (*(trans_type->exec))(plan->libplan_in,in+k*d1[0]*d1[1],tmp);
       }
       else */
	 //	 if((void *) in != (void *) out || dt2 <= dt1)
       (*(trans_type->exec))(plan->libplan_out,in+k*d1[0]*d1[1],tmp);
       
       pout = out+k*d2[0]*d2[1];
#ifdef MKL_BLAS
       blas_trans<Type2>(d2[1],d2[0],1.0,tmp,d2[1],pout,d2[0]);
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
     
     //     if(!inplace)
     delete [] tmp;

     break;     
   }
#endif
   
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot120in(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,bool deriv)
{

  
#ifdef CUDA  

  Type2 *tmp;
  cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
  int imc[3] = {2,0,1};

  (*(exec))(plan->libplan_out[nslices],in,tmp);
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,d2);
  

#ifdef CUTENSOR
  ro_cutensor_in<Type2>(tmp,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>(tmp,out+k,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro120in_cu(gridDim,blockSize,tmp,out,d2);
#endif

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
  



    Type1 *pin_t1,*pIN;
    Type2 *tmp,*pin2,*pout2;
    Type2 *pin,*pin1,*pout,*pout1,*ptran2; 
    int nb13, nb31;
    int i,j,k,ii,jj,kk,i2,j2,k2;

	if(d1[0]*d1[1] >0)
	  nb31 = nb13 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = nb13 = 1;
	if(nb31 < 1) nb31 = 1;
	if(nb13 < 1) nb13 = 1;

	tmp = new Type2[d2[2]*d2[1]*nb31];
	int sdims[3];

	if(deriv) {
	  sdims[0] = d2[1];
	  sdims[1] = d2[2];
	  sdims[2] = nb31;
	}

	for(k=0;k <d1[2];k+=nb31) {
	  k2 = min(k+nb31,d1[2]);
	  // if(inplace) {
	  //  tmp = (Type2 *) (in+k*d1[0]*d1[1]);
	  //  (*(trans_type->exec))(plan->libplan_in,in+k*d1[0]*d1[1],tmp);
	  // }
	  //else
	  (*(trans_type->exec))(plan->libplan_out[nslices],pIN+k*d1[0]*d1[1],tmp);
	  if(deriv)
	    compute_deriv_loc(tmp,tmp,sdims);

	  for(i=0;i < d2[1];i+=nb13) {
	    i2 = min(i+nb13,d2[1]);
	    pout2 =  out + i*d2[0];
	    pin2 =  tmp + i;
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
//	if(!inplace)
	  delete [] tmp;

    }
#endif
}


template <class Type1,class Type2> void transplan<Type1,Type2>::rot210in(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,bool deriv)
{

#ifdef CUDA  

  Type2 *tmp;
  cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
  (*(exec))(plan->libplan_out[nslices],in,tmp);
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,d2);

  int imc[3] = {2,1,0};

#ifdef CUTENSOR
  ro_cutensor_in<Type2>(tmp,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>(tmp,out,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro210in_cu(gridDim,blockSize,tmp,out,d2);
#endif

    cudaFree(tmp);
  
#else 
    {

  Type2 *tmp,*pout,*pin,*pin1,*pin2,*pout1;
  std::unique_ptr<Type2> myptr2;
  int i,ii,i2,j,k,kk,k2,nb31,nb13;


  if(d1[0]*d1[1] >0)
    nb31 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else
    nb31 = 1;
  if(nb31 < 1) nb31 = 1;
  if(d2[0]*d2[1] >0)
    nb13 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else
    nb13 = 1;
  if(nb13 < 1) nb13 = 1;
  
  myptr2 =  std::unique_ptr<Type2> (new Type2[d2[0]*d2[1]*d2[2]]);
  tmp = myptr2.get();
	
      //	tmp = new Type2[d2[2]*d2[1]*d2[0]];
  (*(exec))(plan->libplan_out[nslices],in,tmp);
  if(deriv) {
    int sdims[3] = {d2[2],d2[1],d2[0]};
    compute_deriv_loc(tmp,tmp,sdims);
  }

  for(k=0;k <d2[0];k+=nb31) {
    k2 = min(k+nb31,d2[0]);
    for(i=0;i < d2[2];i+=nb13) {
      i2 = min(i+nb13,d2[2]);
      for(kk=k; kk < k2; kk++) {
	pin1 =  tmp + (kk*d2[2]*d2[1] +i);
	pout1 =  out + (kk + i *d2[1]*d2[0]);
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

    }
#endif
}

template <class Type1,class Type2> void transplan<Type1,Type2>::rot201in(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,bool deriv)
{

#ifdef CUDA  
  Type2 *tmp;
  cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*d2[0]*d2[1]*d2[2]);
  (*(exec))(plan->libplan_out[nslices],in,tmp);
  //  if(deriv)
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,d2);

  int imc[3] = {1,2,0};

#ifdef CUTENSOR
  ro_cutensor_in<Type2>(tmp,out,imc,d2,0,1,0);
#elif defined CUTT
  ro_cutt_in<Type2>(tmp,out,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro201in_cu(gridDim,blockSize,tmp,out,d2);
#endif

    cudaFree(tmp);

#else
      {
      
      Type2 *tmp,*pout,*pin,*pin1,*pin2,*pout1;
  std::unique_ptr<Type2> myptr2;
  int i,j,jj,j2,k,kk,k2,nb32,nb23;


  if(d2[0]*d2[1] > 0)
    nb23 = cache_bl / (sizeof(Type2)*d2[0]*d2[1]);
  else
    nb23 = 1;
  if(nb23 < 1) nb23 = 1;
  if(d2[0]>0)
    nb32 = cache_bl / (sizeof(Type2)*d2[0]*nb23);
  else
    nb32 = 1;
  if(nb32 < 1) nb32 = 1;
  
  myptr2 =  std::unique_ptr<Type2> (new Type2[d2[0]*d2[1]*d2[2]]);
  tmp = myptr2.get();
  
	//	tmp = new Type2[d2[0]*d2[1]*d2[2]];
  (*(exec))(plan->libplan_out[nslices],in,tmp);
  if(deriv) {
    int sdims[3] = {d2[2],d2[0],d2[1]};
    compute_deriv_loc(tmp,tmp,sdims);
  }

  for(k=0;k <d1[1];k+=nb32) {
    k2 = min(k+nb32,d1[1]);
    for(j=0;j < d1[2];j+=nb23) {
      j2 = min(j+nb23,d1[2]);
      for(kk=k; kk < k2; kk++){
	pin1 = tmp + (kk*d2[2] +j*d2[0]*d2[2]);
	pout1 =  out + (kk +j*d2[0]);
	for(jj=j; jj < j2; jj++) {
	  pin = pin1;
	  pout = pout1;
	  for(i=0;i < d2[2];i++) {
	    *pout = *pin++;
	    pout += d2[0]*d2[1];
	  }
	  pin1 += d2[2]*d2[0];
	  pout1+= d2[0];
	}
      }
    }
  }

    }
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
template <class Type1,class Type2> void transplan<Type1,Type2>::rot021in_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,int slice,int nslice,bool deriv)
{

#ifdef CUDA  
  //  if(useCuda) {

  Type2 *tmp;
  cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type2)*mysize1[slice]);
  /*
  if (cudaDeviceSynchronize() != cudaSuccess){
    fprintf(stderr, "Cuda error: Failed to synchronize\n");
    return;
  }
  */
  int imc[3] = {0,2,1};

  (*(exec))(plan->libplan_out[slice],in+offset1[slice],tmp);
  //  if(deriv) {
  //  int sdims[3] = {d2[0],d1[1],mysize1[slice]/(d1[1]*d2[0])};
  //  compute_deriv_loc_cu<<<>>>(tmp,tmp,sdims);
  //}

  /*
  `for(int j=0;j<d1[1]/nslice;j++) {
    int j1 = j + d1[1]*slice/nslice;
    for(int k=0;k<d1[2];k++) 
      (*(exec))(plan->libplan_out[slice],in+d1[0]*(j1+d1[1]*k),tmp+d2[0]*(j+d1[1]*k));
  }
*/

  int myoffset = d2[0]*find_disp(d1[2],slice,nslice);

#ifdef CUTENSOR
  ro_cutensor_in<Type2>(tmp,out+myoffset,imc,d2,slice,nslice,streams[slice]);
#elif defined CUTT
  ro_cutt_in<Type2>(tmp,out,imc,d2);
#else // do our own CUDA transpose
    dim3 gridDim=((d2[1]+1)/TILE_DIM,(d2[0]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
  //  ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro021in_cu(gridDim,blockSize,tmp,out,d2);
    //  rot021in_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif

    cudaFree(tmp);
    cudaDeviceSynchronize();
    //  }
    //else
#else
  {
  Type1 *pin_t1;
  Type2 *tmp,*pout,*pin,*pin1,*pin2,*pout1,*pout2;
  std::unique_ptr<Type2> myptr2;
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
  
  if((void *) in != (void *) out) {
    if(deriv) {
      int sdims[3] = {d2[2],1,1};
    
      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    for(jj=j; jj < j2; jj++) {
	      pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
	      pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
	      (*(exec))(plan->libplan_out,pin_t1,pout2);
	    }
	  }
	}
      }
    }

    else
      for(k=0;k <d1[2];k+=nb32) {
	k2 = min(k+nb32,d1[2]);
	for(j=0;j < d1[1];j+=nb23) {
	  j2 = min(j+nb23,d1[1]);
	  for(kk=k; kk < k2; kk++) {
	    for(jj=j; jj < j2; jj++) {
	      pin_t1 = in +kk*d1[0]*d1[1] +jj*d1[0];
	    pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
	    (*(exec))(plan->libplan_out,pin_t1,pout2);
	    }
	  }
	}
      }
  }
  else { // In-place
    
    myptr2 =  std::unique_ptr<Type2> (new Type2[d2[0]*d2[1]*d2[2]]);
    tmp = myptr2.get();
    //	  tmp = new Type2[d2[0]*d2[1]*d2[2]];
    (*(exec))(plan->libplan_out,in,tmp);
    if(deriv) {
      int sdims[3] = {d2[0],d2[2],d2[1]};
      compute_deriv_loc(tmp,tmp,sdims);
    }

    for(k=0;k <d2[1];k+=nb32) {
      k2 = min(k+nb32,d2[1]);
      for(j=0;j < d2[2];j+=nb23) {
	j2 = min(j+nb23,d2[2]);
	for(kk=k; kk < k2; kk++) {
	  for(jj=j; jj < j2; jj++) {
	    pin2 = tmp +kk*d2[0]*d2[2] +jj*d2[0];
	    pout2 = out + kk * d2[0] + jj * d2[0]*d2[1];
	    memcpy(pout2,pin2,d2[0]);
	  }
	}
      }
    }
    
  }

  }
#endif
}


 template <class Type1,class Type2> void transplan<Type1,Type2>::rot102out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int dt1,int dt2,int slice,int nslice,bool deriv)
{

#ifdef CUDA  

    Type1 *tmp;
    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type1)*mysize1[slice]);
    int imc[3] = {1,0,2};
#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in+offset1[slice],tmp,imc,d1,slice,nslice,streams[slice]);
#elif defined CUTT
    ro_cutt_out<Type1>(in,tmp,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro102out_cu(gridDim,blockSize,in,tmp,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    
    (*(exec))(plan->libplan_out[slice],tmp,out+offset2[slice]);
    if(deriv) {
      int l=d2[2] % nslices;
      int mysz;
      if(slice < l)
	mysz = d2[2]/nslices+1;
      else
	mysz = d2[2]/nslices;
      int sdims[3]={d2[0],d2[1],mysz};
      //      compute_deriv_loc_cu<<<>>>(out+offset2[slice],out+offset2[slice],sdims);
    }
    
    cudaFree(tmp);
#else
    {

  Type1 *tmp,*pin,*pout,*pout1;
  Type2 *pout2;
  int i,j,k;
  std::unique_ptr<Type1> myptr1;

  int sdims[3];
  if(deriv) {
    sdims[0] = d2[0];
    sdims[1] = d2[1];
    sdims[2] = 1;
  }
  if((void *) in == (void *) out && dt2 > dt1) {

    myptr1 =  std::unique_ptr<Type1> (new Type1[d1[0]*d1[1]*d1[2]]);
    tmp = myptr1.get();
    //	tmp = new Type1[d1[0]*d1[1]*d1[2]];
    for(k=0;k <d1[2];k++) {
      pin = in + k*d1[0]*d1[1];
#ifdef MKL_BLAS
      blas_trans<Type1>(d1[0],d1[1],1.0,pin,d1[0],tmp,d1[1]);
#else
      pout1 = tmp + k*d1[0]*d1[1];
      for(j=0;j < d1[1];j++) {
	pout = pout1 + j;
	for(i=0;i < d1[0];i++) {
	  *pout = *pin++;
	  pout += d1[1];
	}	
      }
#endif
      pout2 = out + d2[0]*d2[1]*k;
      (*(exec))(plan->libplan_out,pout1, pout2);
      if(deriv)
	compute_deriv_loc(pout2,pout2,sdims);
    }
    
  }
  
  else {

    myptr1 =  std::unique_ptr<Type1> (new Type1[d1[0]*d1[1]]);
    tmp = myptr1.get();
    //	tmp = new Type1[d1[0]*d1[1]];
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
      (*(exec))(plan->libplan_out,tmp,pout2);
      if(deriv)
	compute_deriv_loc(pout2,pout2,sdims);
    }
    
  }

    }
#endif
}
 

template <class Type1,class Type2> void transplan<Type1,Type2>::rot120out(Type1 *in,Type2 *out,bool inplace,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan, int cache_bl,bool deriv)
{

#ifdef CUDA  
    Type1 *tmp;
    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
    int imc[3] = {2,0,1}; //{1,2,0};


#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,tmp,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,tmp,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
  //  ro120in_cu<Type2>(gridDim,blockSize,tmp,out,d2);
    ro120out_cu(gridDim,blockSize,in,tmp,d1);
#endif
    (*(exec))(plan->libplan_out[nslices],tmp,out);
    //     if(deriv)
    //  compute_deriv_loc_cu<<<>>>(out,out,d2);
   //    cudaDeviceSynchronize();
    cudaFree(tmp);
    
#else
    {
 // no Cuda

  Type1 *tmp,*pin,*pin1,*pout,*pout1;
  Type2 *pout2;
  int i,j,k,ii,i2,kk,k2,nb13,nb31;

	if(d1[0]*d1[1] > 0)
	  nb31 = nb13 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
	else
	  nb31 = nb13 = 1;
	if(nb31 < 1){  nb31 = 1; nb13 = 1; }

	if((void *) in == (void *) out) {

	tmp = new Type1[d1[0]*d1[1]*d1[2]]; // tmp[k][i][j] = in[i][j][k]

	for(j=0;j < d1[1];j++) {
	  pin1 = in + j*d1[0];
	  pout1 = tmp + j*d1[2]*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + i + kk*d1[0]*d1[1];
		pout = pout1 + i*d1[2] +kk;
		
		for(ii=i; ii < i2; ii++) {
                  *pout = *pin++;
		  pout += d1[2];
		}
	      }
	    }
	  }
	}
	if(deriv) {
	  int sdims[3] = {d2[0],d2[1],1};
	  for(j=0;j < d1[1];j++) {
	    pout1 = tmp + j*d1[2]*d1[0];
	    pout2 = out + j*d2[0]*d2[1];
	    (*(trans_type->exec))(plan->libplan_out[nslices],(Type1 *) pout1,pout2);
	    compute_deriv_loc(pout2,pout2,sdims);
	  }
	  else
	    for(j=0;j < d1[1];j++) {
	      pout1 = tmp + j*d1[2]*d1[0];
	      pout2 = out + j*d2[0]*d2[1];
	      (*(trans_type->exec))(plan->libplan_out[nslices],(Type1 *) pout1,pout2);
	    }
	}
	else {
	tmp = new Type1[d1[0]*d1[2]]; // tmp[k][i] = in[i][j][k]

	int sdims[3];
	if(deriv) {
	  sdims[0] = d2[0];
	  sdims[1] = d2[1];
	  sdims[2] = 1;
	}

	for(j=0;j < d1[1];j++) {
	  pin1 = in + j*d1[0];
	  for(k=0;k <d1[2];k+=nb31) {
	    k2 = min(k+nb31,d1[2]);
	    
	    for(i=0;i < d1[0];i+=nb13) {
	      i2 = min(i+nb13,d1[0]);
	      
	      for(kk=k; kk < k2; kk++) {
		
		pin = pin1 + i + kk*d1[0]*d1[1];
		pout =  tmp + i*d1[2] +kk;
		
		for(ii=i; ii < i2; ii++) {
                  *pout = *pin++;
		  pout += d1[2];
		}
	      }
	    }
	  }

	  pout2 = out + j*d2[0]*d2[1];
	  (*(trans_type->exec))(plan->libplan_out[nslices],tmp,pout2);
	  if(deriv)
	    compute_deriv_loc(pout2,pout2,sdims);
	}
	
	}	
	delete [] tmp;

    }
#endif
}




template <class Type1,class Type2> void transplan<Type1,Type2>::rot210out(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,bool deriv)
{


#ifdef CUDA  
    Type1 *tmp;
    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
    int imc[3] = {2,1,0};
#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,tmp,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,tmp,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro210out_cu(gridDim,blockSize,in,tmp,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    (*(exec))(plan->libplan_out[nslices],tmp,out);
    cudaFree(tmp);
    
#else
  {
 // no Cuda

  Type1 *tmp,*pin,*pin1,*pout,*pout1,*pout2;
  int i,j,k,ii,i2,kk,k2,nb13,nb31;
  std::unique_ptr<Type1> myptr1;


  if(d1[0]*d1[1] >0)
    nb31 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else nb31 = 1;
  if(nb31 < 1) nb31 = 1;
  if(d2[0]*d2[1] >0)
    nb13 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else nb13 = 1;
  if(nb13 < 1) nb13 = 1;
	
  myptr1 =  std::unique_ptr<Type1> (new Type1[d1[0]*d1[1]*d1[2]]);
  tmp = myptr1.get();
	//	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
  for(k=0;k <d1[2];k+=nb31) {
    k2 = min(k+nb31,d1[2]);
    for(i=0;i < d1[0];i+=nb13) {
      i2 = min(i+nb13,d1[0]);
      for(kk=k; kk < k2; kk++) {
	pin1 = in + (kk*d1[0]*d1[1] +i);
	pout1 = tmp + (kk +i*d1[2]*d1[1]);
	for(j=0;j < d1[1];j++) {
	  pin = pin1;
	  pout = pout1;
	  for(ii=i; ii < i2; ii++) {
	    *pout = *pin++;
	    pout += d1[2]*d1[1];
	  }
	  pin1 += d1[0];
	  pout1 += d1[2];
	}
      }
    }
  }

#ifdef DEBUG
  printf("(2,1,0): Transform of length %d, batch size %d, in array sized %d\n",N,m,d1[2]*d1[0]*d1[1]);
#endif
  (*(exec))(plan->libplan_out[nslices],tmp,out);
  if(deriv)
    compute_deriv_loc(out,out,d2);
  }
#endif
}
	

template <class Type1,class Type2> void transplan<Type1,Type2>::rot201out(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,bool deriv)
{

#ifdef CUDA  
    Type1 *tmp;
    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type1)*d1[0]*d1[1]*d1[2]);
    int imc[3] = {1,2,0};
#ifdef CUTENSOR
    ro_cutensor_out<Type1>(in,tmp,imc,d1,0,1,0);
#elif defined CUTT
    ro_cutt_out<Type1>(in,tmp,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro210out_cu(gridDim,blockSize,in,tmp,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif
    (*(exec))(plan->libplan_out[nslices],tmp,out);
    cudaFree(tmp);
   
#else
    {
// no Cuda

  Type1 *tmp,*pin,*pin1,*pout,*pout1,*pout2;
  int i,j,k,jj,j2,kk,k2,nb23,nb32;
  std::unique_ptr<Type1> myptr1;

  
  if(d1[0]*d1[1] >0)
    nb32 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  if(nb32 < 1) nb32 = 1;
  nb23 = nb32;
  
  myptr1 =  std::unique_ptr<Type1> (new Type1[d1[0]*d1[1]*d1[2]]);
  tmp = myptr1.get();
  //	tmp = new Type1[d1[1]*d1[2]*d1[0]];
	
  for(k=0;k <d1[1];k+=nb32) {
    k2 = min(k+nb32,d1[1]);
    for(j=0;j < d1[2];j+=nb23) {
      j2 = min(j+nb23,d1[2]);
      for(kk=k; kk < k2; kk++){
	pin1 = in + (kk*d1[0] +j*d1[0]*d1[1]);
	pout1 = tmp + (kk +j*d1[1]);
	for(jj=j; jj < j2; jj++) {
	  pin = pin1;
	  pout = pout1;
	  for(i=0;i < d1[0];i++) {
	    *pout = *pin++;
	    pout += d1[2]*d1[1];
	  }
	  pin1 += d1[0]*d1[1];
	  pout1 += d1[1];
	}
      }
    }
  }

  (*(exec))(plan->libplan_out[nslices],tmp,out);
  if(deriv)
    compute_deriv_loc(out,out,d2);
	//	delete [] tmp;
    }
#endif
}


template <class Type1,class Type2> void transplan<Type1,Type2>::rot021out_slice(Type1 *in,Type2 *out,int d1[3],int d2[3],void (*exec)(...),  Plantype<Type1,Type2> *plan,int cache_bl,int slice,int nslice,bool deriv)
{

#ifdef CUDA  

    Type1 *tmp;
    cudaMalloc(reinterpret_cast<void **> (&tmp),sizeof(Type1)*mysize1[slice]);
    int imc[3] = {0,2,1};

int myoffset = d1[0] * find_disp(d1[2],slice,nslice);

#ifdef CUTENSOR
ro_cutensor_out<Type1>(in+offset1[slice],tmp+myoffset,imc,d1,slice,nslice,streams[slice]);
#elif defined CUTT
    ro_cutt_out<Type1>(in,tmp,imc,d1);
#else // do our own CUDA transpose
    dim3 gridDim=((d1[0]+1)/TILE_DIM,(d1[2]+1)/TILE_DIM);
    dim3 blockSize=(TILE_DIM,TILE_DIM);
    ro021out_cu(gridDim,blockSize,in,tmp,d1);
    //  rot102out_cu<<<gridDim,blockSize>>>(in,out,d1,d2,exec,plan);
#endif

    int d[3];
    d[0] = 1;
    d[1] = 1;
    d[2] = d2[2];
    int offset[nslices],mysize[nslices];
    divide_work(d,offset,mysize,nslices);
    if(deriv) {
      int sdims[3];
      sdims[2] = 1;
      sdims[1] = d2[1];
      sdims[0] = d2[0];

      for(int k=0;k<mysize[slice];k++) {
	int k1 = k + offset[slice];
	for(int j=0;j>d2[1];j++) {
	  (*(exec))(plan->libplan_out[slice],tmp+d1[0]*(j+d2[1]*k),out+d2[0]*(j+d2[1]*k1));
	  //	  compute_deriv_loc_cu<<<>>>(out+d2[0]*(j+d2[1]*k1),out+d2[0]*(j+d2[1]*k1),sdims);
	}
      }
    }
    else
      for(int k=0;k<mysize[slice];k++) {
	int k1 = k + offset[slice];
	for(int j=0;j>d2[1];j++) 
	  (*(exec))(plan->libplan_out[slice],tmp+d1[0]*(j+d2[1]*k),out+d2[0]*(j+d2[1]*k1));
      }
    
    cudaFree(tmp);
    cudaDeviceSynchronize();
#else
    {
 // no Cuda

  Type1 *tmp,*pin,*pin1,*pout,*pout1,*pin2,*pin3,*pout3,*pout4;
  Type2 *pout2,*ptran2;
  int i,j,k,jj,j2,kk,k2,nb23,nb32;
  std::unique_ptr<Type1> myptr1;

  if(d1[0]*d1[1] >0)
    nb32 = cache_bl / (sizeof(Type1)*d1[0]*d1[1]);
  else nb32 = 1;
  if(nb32 < 1) nb32 = 1;
  if(d2[0]*d2[1] >0)
    nb23 = cache_bl / (sizeof(Type1)*d2[0]*d2[1]);
  else nb23 = 1;
  if(nb23 < 1) nb23 = 1;
	
#ifdef DEBUG
  printf("(0,2,1): Transform of length %d, batch size %d, array leading dim %d\n",N,m,d1[0]);
#endif

  if((void *) in == (void *) out) {
    
    myptr1 =  std::unique_ptr<Type1> (new Type1[d1[0]*d1[1]*d1[2]]);
    tmp = myptr1.get();
	//	tmp = new Type1[d1[0]*d1[1]*d1[2]]; // tmp[k][i] = in[i][j][k]
	
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
	
    if(deriv) {
       int sdims[3] = {d2[0],1,1};
       for(k=0;k<d2[2];k++)
	 for(j=0;j<d2[1];j++) {
	    pin2 = tmp+d1[0]*(j+k*d2[1]);
	    pout2 = out + d2[0]*(j+k*d2[1]);
	    (*(exec))(plan->libplan_out,pin2,pout2);
	    compute_deriv_loc(pout2,pout2,sdims);
	}
    }
    else
      for(k=0;k<d2[2];k++)
	for(j=0;j<d2[1];j++) {
	  pin2 = tmp+d1[0]*(j+k*d2[1]);
	  pout2 = out + d2[0]*(j+k*d2[1]);
	  (*(exec))(plan->libplan_out,pin2,pout2);
	}
    
  }
  else {
	
    if(deriv) {
       int sdims[3] = {d2[0],1,1};

       for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	     j2 = min(j+nb23,d1[1]);
	     for(kk=k; kk < k2; kk++) {
	        pin2 = in + kk*d1[0]*d1[1] +j*d1[0];
		ptran2 = out +kk*d2[0] +j*d2[0]*d2[1];
		for(jj=j; jj < j2; jj++) {
	          (*(exec))(plan->libplan_out,pin2,ptran2);
		  compute_deriv_loc(ptran2,ptran2,sdims);
		  pin2 += d1[0];
		  ptran2 += d2[0]*d2[1];
	        }
	      }
           }
	}
    }
    else
      for(k=0;k <d1[2];k+=nb32) {
	  k2 = min(k+nb32,d1[2]);
	  for(j=0;j < d1[1];j+=nb23) {
	     j2 = min(j+nb23,d1[1]);
	     for(kk=k; kk < k2; kk++) {
	        pin2 = in + kk*d1[0]*d1[1] +j*d1[0];
		ptran2 = out +kk*d2[0] +j*d2[0]*d2[1];
		for(jj=j; jj < j2; jj++) {
	           (*(exec))(plan->libplan_out,pin2,ptran2);
		   pin2 += d1[0];
		   ptran2 += d2[0]*d2[1];
	        }
	     }
           }
       }
    }	

  }
#endif
}

