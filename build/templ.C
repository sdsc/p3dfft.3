/*
!
!    P3DFFT++
!
!    Software Framework for Scalable Fourier Transforms in Three Dimensions
!
!    Copyright (C) 2017 Dmitry Pekurovsky
!    Copyright (C) 2017 University of California
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!
!----------------------------------------------------------------------------
*/

#include "p3dfft.h"

namespace p3dfft {

void inv_mo(int mo[3],int imo[3]);

template<class Type1,class Type2> transform3D<Type1,Type2>::transform3D(const grid& grid1_, const grid& grid2_,const trans_type3D *type,bool inplace_)
{

#ifdef DEBUG
  cout << "In transform3D" << endl;
  <print_type3D(type);
#endif

  int prec2,dt;
  Stages = NULL;is_set = false;
  if(typeid(Type1) == type_float) {
    prec = 4;
    dt = 1;
  }
  else if(typeid(Type1) == type_double) {
    prec = 8;
    dt = 1;
  }
  else if(typeid(Type1) == type_complex) {
    prec = 4;
    dt = 2;
  }
  else if(typeid(Type1) == type_complex_double) {
    prec = 8;
    dt = 2;
  }

  if(typeid(Type2) == type_float) {
    prec2 = 4;
  }
  else if(typeid(Type2) == type_double) {
    prec2 = 8;
  }
  else if(typeid(Type2) == type_complex) {
    prec2 = 4;
  }
  else if(typeid(Type2) == type_complex_double) {
    prec2 = 8;
  }
  if(prec != prec2)
    cout << "Error in transform3D: precisions don't match!" << endl;
  int dt_init = dt;

  int pgrid1[3],pgrid2[3],nstages,pgrid[3],proc_order[3],gdims[3],L[3],Lfin,st;
  int ns,d1,d2,nd,dt1,dt2,i,df,L1,L2,L3,excl(int,int),dist(int);
  void swap0(int new_mo[3],int mo[3],int L);
  int monext[3];
  int dims1[3],dims2[3];
  gen_trans_type *tmptype;
  MPI_Comm mpicomm,splitcomm;
  bool reverse_steps;

  /*
  if(!grid1.is_set || !grid2.is_set) {
    printf("Error in tran3D_plan: grid is not set up\n");
    return;
  }
  */

  if((nd = grid1_.nd) != grid2_.nd) {
    printf("ERror in tran3D_plan: dimensions of grids don't match %d %d\n",nd,grid2_.nd);
    MPI_Abort(MPI_COMM_WORLD,0);
  }
  for(i=0; i < nd; i++)
    if(grid1_.P[i] != grid2_.P[i] || grid1_.proc_order[i] != grid2_.proc_order[i]) {
      printf("Error in transform3D: processor grids dont match: %d %d %d\n",i,grid1_.P[i],grid2_.P[i]); 
      MPI_Abort(MPI_COMM_WORLD,0);
    }

  /*
    int tmp;
    MPI_Comm_size(grid1_.mpicomm[0],&tmp);
    printf("size of grid1_ mpicomm=%d\n",tmp);
    MPI_Comm_size(grid1->mpicomm[0],&tmp);
    printf("size of grid1 mpicomm=%d\n",tmp);
  */
  grid1 = new grid(grid1_);
  grid2 = new grid(grid2_);
  inplace = inplace_;
  dt = dt_init;
  dt1 = dt;

  memcpy(pgrid1,grid1_.pgrid,sizeof(int)*3);
  memcpy(pgrid2,grid2_.pgrid,sizeof(int)*3);

  //  nstages = 0;
  stage *prev_stage,*curr_stage;
  int dt_prev = dt;
  bool inpl = false;
  int mocurr[3],mo1[3],mo2[3];
  for(i=0; i < 3; i++) {
    mocurr[i] = mo1[i] = grid1_.mem_order[i];
    mo2[i] = grid2_.mem_order[i];
  }

  mpicomm = grid1_.mpi_comm_glob;

  grid *tmpgrid0 = new grid(grid1_);
  grid *tmpgrid1;
  reverse_steps = false;

  /* Find the order of the three transforms, attempting to minimize reordering and transposes */ 

  L[0] = grid1_.L[0];
  if(nd == 1 && mocurr[L[0]] != 0)
    if(mocurr[grid1_.L[1]] == 0)
      L[0] = grid1_.L[1];
      
  Lfin = L[2] = grid2_.L[0];
  if(L[2] == L[0])
    if(nd == 1)
      Lfin = L[2] = grid2_.L[1];
    else {
      reverse_steps=true;
      L[2] = dist(L[0]);
    }
  L[1] = excl(L[0],L[2]);

  for(i=0;i<3;i++) 
    monext[i] = mocurr[i];

#ifdef DEBUG
  printf("%d: Planning stages: %d %d %d\n",grid1_.taskid,L[0],L[1],L[2]);
#endif

  // Plan the stages
  for(int st=0;st < 3;st++) {
    
    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      gdims[i] = tmpgrid0->gdims[i];
      proc_order[i] = tmpgrid0->proc_order[i];
      pgrid[i] = tmpgrid0->pgrid[i];
    }

    // Determine if an MPI transpose is involved in this stage
    d1 = -1;
    if(st < 2) { 
      swap0(monext,mocurr,L[st+1]);
      for(i=0;i<nd;i++)
	if(L[st+1] == tmpgrid0->D[i]) {
	  d1 = L[st+1];
	  d2 = L[st];
	  splitcomm = tmpgrid0->mpicomm[i];
	  break;
	}
   }
    else   {
      for(i=0;i<3;i++) 
	monext[i] = grid2_.mem_order[i];

      for(i=0;i<nd;i++)
	if(L[st] == tmpgrid0->D[i]) {
	  d1 = L[st];
	  d2 = Lfin;
	  splitcomm = tmpgrid0->mpicomm[i];
	  break;
	}
    }

    //    tmpgrid0->set_mo(monext);
      
    tmptype = types1D[type->types[L[st]]];
    if(tmptype->dt1 < tmptype->dt2) { // Real-to-complex
      gdims[L[st]] = gdims[L[st]]/2+1;
    }
    else if(tmptype->dt2 < tmptype->dt1) { // Complex-to-real
      gdims[L[st]] = (gdims[L[st]]-1)*2;
    }

    inpl = false;

    if(d1 >= 0) { // If a transpose is involved

      // Set up the new grid
      pgrid[d1] = 1;
      pgrid[d2] = tmpgrid0->pgrid[d1];  ////P[0];  // proc_order[0] ???

      tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

#ifdef DEBUG
      printf("Calling init_tran_MPIsplan, trans_dim=%d, d1=%d, d2=%d, gdims2=(%d %d %d), ldims2=(%d %d %d), mem_order=(%d %d %d)\n",L[st],d1,d2,gdims[0],gdims[1],gdims[2],tmpgrid1->ldims[0],tmpgrid1->ldims[1],tmpgrid1->ldims[2],monext[0],monext[1],monext[2]);
#endif

      curr_stage = init_trans_MPIplan(*tmpgrid0,*tmpgrid1,splitcomm,d1,d2,tmptype,L[st],inpl,prec);
      curr_stage->kind = TRANSMPI;
    }
    else { // Only transform

      tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);
#ifdef DEBUG
      printf("Calling init_transplan, trans_dim=%d, gdims2=(%d %d %d), ldims2=(%d %d %d), mem_order=(%d %d %d)\n",L[st],gdims[0],gdims[1],gdims[2],tmpgrid1->ldims[0],tmpgrid1->ldims[1],tmpgrid1->ldims[2],monext[0],monext[1],monext[2]);
#endif

      curr_stage = init_transplan(*tmpgrid0,*tmpgrid1,tmptype,L[st],inpl,prec);

      curr_stage->kind = TRANS_ONLY;
      
    }
    dt_prev = tmptype->dt2;
    delete tmpgrid0;
    tmpgrid0 = tmpgrid1;
    if(st == 0) 
      Stages = prev_stage = curr_stage;
    else {
      prev_stage->next = curr_stage;
      prev_stage = curr_stage;
    }

  }


    //If needed, transpose back to the desired layout, specified by grid2

    if(reverse_steps) {

    d1 = tmpgrid1->D[1];
    df = grid2_.D[1];

#ifdef DEBUG
    cout << "Return steps" << endl;
#endif
    if(d1 != df) {
      // Exchange D1 with L0
      d2 = L1 = tmpgrid1->L[0];
      if(d1 == grid1_.D[0])
	splitcomm = grid1_.mpicomm[0];
    else
      splitcomm = grid1_.mpicomm[1];
      pgrid[d1] = 1;
      pgrid[d2] = tmpgrid1->pgrid[d1];

      tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

      prev_stage = curr_stage;
      curr_stage = init_MPIplan(*tmpgrid1,*tmpgrid0,splitcomm,d1,d2,dt_prev,prec);
      prev_stage->next = curr_stage;
      curr_stage->kind = MPI_ONLY;
      delete tmpgrid1;
      L1 = d1;
      d1 = d2;
    }
    else
      tmpgrid0 = tmpgrid1;
    
    L1 = tmpgrid0->L[0];
    L2 = grid2_.L[0];
    if(d1 != df || L1 != L2) {
	// Exchange L0 with D0
	d2 = tmpgrid0->L[0];
	d1 = tmpgrid0->D[0];
	if(d1 == grid1_.D[0])
	  splitcomm = grid1_.mpicomm[0];
        else
          splitcomm = grid1_.mpicomm[1];
        pgrid[d1] = 1;
	pgrid[d2] = tmpgrid0->pgrid[d1];
	
	tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);
	
	prev_stage = curr_stage;
	curr_stage = init_MPIplan(*tmpgrid0,*tmpgrid1,splitcomm,d1,d2,dt_prev,prec);
	prev_stage->next = curr_stage;
	curr_stage->kind = MPI_ONLY;
	delete tmpgrid0;	
    }
      
    d1 = tmpgrid1->D[1];
    if(d1 != df) {
      // Again exchange D1 with L0
      d2 = L1 = tmpgrid1->L[0];
      if(d1 == grid1_.D[0])
	splitcomm = grid1_.mpicomm[0];
      else
	splitcomm = grid1_.mpicomm[1];
      pgrid[d1] = 1;
      pgrid[d2] = tmpgrid1->pgrid[d1];

      tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

      prev_stage = curr_stage;
      curr_stage = init_MPIplan(*tmpgrid1,*tmpgrid0,splitcomm,d1,d2,dt_prev,prec);
      prev_stage->next = curr_stage;
      curr_stage->kind = MPI_ONLY;
      delete tmpgrid1;
      L1 = d1;
      if(L1 != L2) {
	// Exchange L0 with D0
	d2 = tmpgrid0->L[0];
	d1 = tmpgrid0->D[0];
	pgrid[d1] = 1;
	pgrid[d2] = tmpgrid0->pgrid[d1];
	
	tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);
	
	prev_stage = curr_stage;
	curr_stage = init_MPIplan(*tmpgrid0,*tmpgrid1,tmpgrid0->mpicomm[0],d1,d2,dt_prev,prec);
	prev_stage->next = curr_stage;
	curr_stage->kind = MPI_ONLY;
	delete tmpgrid0;	
      }
    }
    }

      for(i=0;i<3;i++)
      if(pgrid[i] != grid2_.pgrid[i]) {
	cout << "Error in transform3D: processor grids dont match" <<endl;
	return;
      }

    bool iseq = true;
    for(i=0; i < 3; i++) 
      if(monext[i] != grid2_.mem_order[i]) {
	iseq = false;
	break;
      }

    if(!iseq) { //If not in the final memory ordering
      tmpgrid0 = new grid(grid2_);
      prev_stage = curr_stage;
#ifdef DEBUG
      printf("Calling init_transplan, trans_dim=%d\n",L2);
#endif
      curr_stage = init_transplan(*tmpgrid1,*tmpgrid0,types1D[EMPTY_TYPE],L2,inpl,prec);
      curr_stage->kind = TRANS_ONLY;
      prev_stage->next = curr_stage;
      delete tmpgrid0;
    }
    delete tmpgrid1;

  if(nd == 3) {
    cout << "Three-dimensional decomposition is presently not supported" << endl;
    return;
  }

  //  cout << "Done transform3D planning" << endl;
  is_set = true;
  return;

}

template<class Type1,class Type2> transform3D<Type1,Type2>::~transform3D()
{
  stage *pnext,*p=Stages;
  if(p)
    while(pnext = p->next) {
      delete p;
      p = pnext;
    }
  delete grid1,grid2;
}

int excl(int a,int b)
{
  if(a*b == 0) {
    if(max(a,b) == 1)
      return(2);
    else
      return(1);
  }
  else
    return(0);
}

int dist(int a)
{
  switch(a) {
  case 0:
  case 1:
    return(2);
  case 2:
    return(0);
  }
}


template <class Type1,class Type2> transplan<Type1,Type2>::transplan(const grid &gr1,const grid &gr2, const gen_trans_type *type,int d, bool inplace_) 
{
  lib_plan = 0;plan = NULL;fft_flag = DEF_FFT_FLAGS;
  if(!type->is_set) {
    cout << "Error in trans_plan: 1D transform type not set" << endl;
    return;
  }
  if(gr1.ldims[d] != gr1.gdims[d] || gr2.ldims[d] != gr2.gdims[d] ) {
    printf("Error in transplan: dimensions dont match %d, %d, %d\n",gr1.ldims[d],gr2.ldims[d],d);
    return;
  }
  stage_prec = prec = type->prec;
  trans_dim = d; 
  trans_type = (trans_type1D<Type1,Type2> *) type;
  dt1 = type->dt1;
  dt2 = type->dt2;
  kind = TRANS_ONLY;
  inplace = inplace_;
  grid1 = new grid(gr1);
  grid2 = new grid(gr2);

  istride = 1;ostride = 1; 
  idist,odist;
  isign = type->isign;

  for(int i=0;i<3;i++) {
    dims1[i] = gr1.ldims[i];
    dims2[i] = gr2.ldims[i];
    mo1[i] = gr1.mem_order[i];
    mo2[i] = gr2.mem_order[i];
  }
  //  inembed = onembed = (int*) &grid1.gdims[d];
  if(type->dt1 < type->dt2) { //Real to complex
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = N;
    odist = N/2+1;
    /*
    dims2[0] = dims2[0]/2+1;
    mygrid->ldims[d] = (mygrid->ldims[d]+2)/2;
    mygrid->gdims[d] = (mygrid->gdims[d]+2)/2;
    mygrid->sz[0][d] = (mygrid->sz[0][d]+2)/2;  
    mygrid->en[0][d] = (mygrid->en[0][d]+2)/2;  
    mygrid->st[0][d] = (mygrid->st[0][d]+2)/2; 
    */
  }
  else if(type->dt1 > type->dt2) { //Complex to real
    N=gr2.gdims[d];
    m=dims2[0]*dims2[1]*dims2[2]/N;
    odist = N;
    idist = N/2+1;
    /*
    dims2[0] = dims2[0]*2-2;
    mygrid->ldims[d] = mygrid->ldims[d]*2-2;
    mygrid->gdims[d] = mygrid->gdims[d]*2-2;
    mygrid->sz[0][d] = mygrid->sz[0][d]*2-2;
    mygrid->st[0][d] = mygrid->st[0][d]*2-2;
    mygrid->en[0][d] = mygrid->en[0][d]*2-2;
    */
  }
  else { // No change in datatype
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = odist = N;
  }

  m = find_m(mo1,mo2,dims1,dims2,trans_dim);

  if(idist <= gr1.gdims[d]) 
    inembed = (int*) &(grid1->gdims[d]);
  else  {
    printf("Error in transplan: dimension too small %d, N=%d\n",gr1.gdims[d],N);
    return;
  }
  if(odist <= gr2.gdims[d]) 
    onembed = (int*) &(grid2->gdims[d]);
  else {
    printf("Error in transplan: dimension too small %d, N=%d\n",dims2[d],N);
    return;
  }

  lib_plan = find_plan(trans_type); 

}  

template <class Type1,class Type2> transplan<Type1,Type2>::transplan(const grid &gr1,const grid &gr2,int type_ID,int d, bool inplace_) 
{

  lib_plan = 0;plan = NULL;fft_flag = DEF_FFT_FLAGS;
  if(gr1.ldims[d] != gr1.gdims[d] || gr2.ldims[d] != gr2.gdims[d] ) {
    printf("Error in transplan: dimensions dont match %d %d %d\n",gr1.ldims[d],gr2.ldims[d],d);
    return;
  }

  if(gr1.pgrid[d] != 1 || gr1.pgrid[d] != 1) {
    printf("Error in transplan: transform dimension %d must be local.\n",d);
    return;
  }

  trans_type = (trans_type1D<Type1,Type2> *) types1D[type_ID];

  if(!trans_type || !trans_type->is_set) {
    cout << "Error in trans_plan: 1D transform type no set" << endl;
    return;
  }

  dt1 = trans_type->dt1;
  dt2 = trans_type->dt2;
  stage_prec = prec = trans_type->prec;
  kind = TRANS_ONLY;
  trans_dim = d;
  inplace = inplace_;
  grid1 = new grid(gr1);
  grid2 = new grid(gr2);

  istride = 1;ostride = 1; 
  //  idist,odist;
  isign = trans_type->isign;

  for(int i=0;i<3;i++) {
    dims1[i] = gr1.ldims[i];
    dims2[i] = gr2.ldims[i];
    mo1[i] = gr1.mem_order[i];
    mo2[i] = gr2.mem_order[i];
  }

  if(trans_type->dt1 < trans_type->dt2) { //Real to complex
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = N;
    odist = N/2+1;
  }
  else if(trans_type->dt1 > trans_type->dt2) { //Complex to real
    N=gr2.gdims[d];
    m=dims2[0]*dims2[1]*dims2[2]/N;
    odist = N;
    idist = N/2+1;
  }
  else { // No change in datatype
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = odist = N;
  }

  m = find_m(mo1,mo2,dims1,dims2,trans_dim);

  if(idist <= gr1.gdims[d]) 
    inembed = (int*) &(grid1->gdims[d]);
  else  {
    printf("Error in transplan: dimension too small %d, N=%d\n",gr1.gdims[d],N);
    return;
  }
  if(odist <= gr2.gdims[d]) 
    onembed = (int*) &(grid2->gdims[d]);
  else {
    printf("Error in transplan: dimension too small %d, N=%d\n",dims2[d],N);
    return;
  }

  lib_plan = find_plan(trans_type); 

}

#define TRANS_IN 0
#define TRANS_OUT 1


template <class Type1,class Type2> int transplan<Type1,Type2>::find_m(int *mo1,int *mo2,int *dims1,int *dims2, int trans_dim) {

  int i,m,mc[3],imo1[3],imo2[3],d1[3],d2[3], scheme;    //,rel_change(int [3],int [3],int [3]);
  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  for(i=0;i<3;i++) {
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }

  rel_change(imo1,imo2,mc);

  if(mo1[trans_dim] == 0) 
    scheme = TRANS_IN;
  else if(mo2[trans_dim] == 0) 
    scheme = TRANS_OUT;
  else {
    printf("Error in reorder_trans: expected dimension %d to be the leading dimension for input or output\n",trans_dim);
    return(-1);
  }

  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2                                                          
      if(scheme == TRANS_IN) 
	m = d2[0]; // Need finer grain transform, to reuse cache     
      else
	m = d1[0];
      break;
    case 2: // 1,2,0
      if(scheme == TRANS_IN) {
	int nb31 = CACHE_BL / (sizeof(Type1)*d1[0]*d1[1]);
	if(nb31 < 1) nb31 = 1;
	m = nb31 * d1[1];
      }
      else
	m = d2[1];
      break;
    }
    break;
    
  case 0:
    if(mc[1] == 1) // 012
      m = d1[1]*d1[2];
    else
      m = 1;
    break;

  case 2:
    if(scheme == TRANS_IN)
      m =d1[1]*d1[2];
    else
      m = d2[1]*d2[2];
    break;
  }
	
  return(m);

}

template <class Type1,class Type2> trans_MPIplan<Type1,Type2>::trans_MPIplan(const grid &gr1,const grid &intergrid, const grid &gr2,MPI_Comm mpicomm,int d1,int d2,const gen_trans_type *type,int d,bool inplace_)   
{
  double *C = new double[1024];
  kind = TRANSMPI;
  trplan = new transplan<Type1,Type2>(gr1,intergrid,type,d,inplace_);
  stage_prec = trplan->prec;
  mpiplan = new MPIplan<Type2>(intergrid,gr2,mpicomm,d1,d2,stage_prec);
  is_set = true;
  inplace = inplace_;
  dt1 = trplan->dt1;
  dt2 = trplan->dt2;
  memcpy(dims1,trplan->dims1,3*sizeof(int));
  memcpy(dims2,mpiplan->dims2,3*sizeof(int));

  //  is_trans = is_mpi = true;
}

/*
template <class Type1,class Type2> trans_MPIplan<Type1,Type2>::~trans_MPIplan()
{
  delete [] SndCnts,SndStrt,RcvCnts,RcvStrt;
  delete grid1,grid2;
}
*/

template <class Type> MPIplan<Type>::MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm_,int d1_,int d2_, int prec_)  
{
  int i,d3,l;


  if(!gr1.is_set || !gr2.is_set) {
    printf("Error in MPIplan constr.: grids are not set up\n");
    return;
  }
  for(i=0; i < 3; i++){
    mo1[i] = gr1.mem_order[i];
    mo2[i] = gr2.mem_order[i];
    if(gr1.gdims[i] != gr2.gdims[i]) {
      printf("Error in MPIplan constr.: global grid dimensions dont match: (%d %d %d) vs. (%d %d %d)\n", gr1.gdims[0],gr1.gdims[1],gr1.gdims[2],gr2.gdims[0],gr2.gdims[1],gr2.gdims[2]);
      return;
    }
  }

  int p=gr1.pgrid[d1_];
  if(p != gr2.pgrid[d2_]) {
    cout << "Error in MPIplan constr.: proc. grid dimensions dont match" << gr1.pgrid << gr2.pgrid << endl;
    return;
  }

  stage_prec = prec = prec_;

  if(p >1) {

  MPI_Comm_size(mpicomm_,&numtasks);
  if(p != numtasks) {
    cout << "Error in MPIplan constr.: proc. grid dimension doesnt match communicator size" << p << numtasks << endl;
    return;
  }

  MPI_Comm_rank(gr1.mpi_comm_glob,&taskid);
  //  MPI_Comm_rank(mpicomm_,&commid);



  SndCnts = new int[p];
  SndStrt = new int[p];
  RcvCnts = new int[p];
  RcvStrt = new int[p];
  
  d1 = d1_;
  d2 = d2_;

  for(i=0; i < 3; i++)
    if(i != d1 && i != d2)
      d3 = i;

  for(l=0;l<gr1.nd;l++) 
    if(gr1.D[gr1.proc_order[l]] == d1)
      break;

  comm_id = l;

  memcpy(dims1,gr1.ldims,3*sizeof(int));
  memcpy(dims2,gr2.ldims,3*sizeof(int));

  // int comm_coords[3];
  //memcpy(comm_coords,gr1.grid_id_cart,3*sizeof(int));

  SndStrt[0] = 0;
  RcvStrt[0] = 0;

  int sz=sizeof(Type)/4;
  //int rank;
  for(int j=0; j< p-1;j++) {
    //    comm_coords[l] = j;
    //MPI_Cart_rank(mpicomm_,comm_coords,&rank); 
    SndCnts[j] = gr2.sz[l][j][d2] * dims1[d3] *dims1[d1] * sz;
    RcvCnts[j] = dims2[d2] * gr1.sz[l][j][d1] *dims2[d3] * sz;
    SndStrt[j+1] = SndStrt[j] + SndCnts[j];
    RcvStrt[j+1] = RcvStrt[j] + RcvCnts[j];
  }
  //  comm_coords[l] = p-1;
  //MPI_Cart_rank(mpicomm_,comm_coords,&rank); 
  SndCnts[p-1] = gr2.sz[l][p-1][d2] * dims1[d3]*dims1[d1] * sz;
  RcvCnts[p-1] = dims2[d2] * gr1.sz[l][p-1][d1]*dims2[d3] * sz;

  grid1 = new grid(gr1);
  grid2 = new grid(gr2);
  du = d3;
  //  is_mpi = true;
  }

  is_set = true;
  //  is_trans = false;
  kind = MPI_ONLY;
  mpicomm = mpicomm_;

}

template <class Type> MPIplan<Type>::~MPIplan()
{
  delete [] SndCnts,SndStrt,RcvCnts,RcvStrt;
  delete grid1,grid2;
}


template <class Type1,class Type2> inline long transplan<Type1,Type2>::find_plan(trans_type1D<Type1,Type2> *type)
{
  int i;
  int planID;
  Plantype<Type1,Type2> *pl;
  //  Plan *p;

#ifdef DEBUG
  printf("find_plan: N,m=%d,%d\n",N,m);
#endif

  if(inplace) 
    if(sizeof(Type1) != sizeof(Type2) || istride != ostride || (m > 1 && idist != odist)) {
      cout << "Error in find_plan: inplace transforms should have identical dimensions and types" << endl;
      return(-1);
    }
  planID = 0;
  for(vector<Plan*>::iterator it=Plans.begin(); it < Plans.end();it++,planID++) {
    //    pl = dynamic_cast<Plantype<Type1,Type2> *> (p);
    //if(pl)
    pl = (Plantype<Type1,Type2> *) *it;
    if(type->dt1 == pl->dt1 && type->dt2 == pl->dt2 && prec == pl-> prec &&
       pl->N == N && pl->m == m && pl->inplace == inplace &&	\
       pl->istride == istride && pl->ostride == ostride && pl->isign == isign &&
       pl->fft_flag == fft_flag) {
      if(m > 1) {
	if(pl->idist = idist && pl->odist == odist &&
	   *(pl->inembed) == *inembed && *(pl->onembed) == *onembed) {
	  plan = pl;
	  return((*it)->libplan);
	}
	//	  return(plan->libplan);
      }
      else {
	plan = pl;
	return((*it)->libplan);
      }
    }
    //	return(plan->libplan);
  }

    // If haven't found existing plan with suitable params, define a new one
    //    plan = &DefPlans[i];


#ifdef DEBUG
  cout << "new plantype" << endl;
#endif
 plan = new Plantype<Type1,Type2>(type->doplan,type->exec,N,m,inplace,istride,idist,ostride,odist,inembed,onembed,isign,fft_flag);
  Plans.push_back(plan);
  
  if(inplace) {
    Type1 *A;
    int size=max(sizeof(Type1)*(istride*N+idist*m),sizeof(Type2)*(ostride*N+odist*m));
#ifdef DEBUG
    printf("in-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",istride,ostride,idist,odist);
#endif
#ifdef FFTW
    A = (Type1 *) fftw_malloc(size *sizeof(Type1));
    Type2 *B = (Type2 *) fftw_malloc(size *sizeof(Type2));
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,inembed,istride,idist,B,onembed,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      //      if(isign == 0) 
      //	cout << "Error in find_plan: isign is not set" << endl;
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign);
    }
    else //R2C or C2R
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,inembed,istride,idist,B,onembed,ostride,odist,fft_flag);

    fftw_free(A);
    fftw_free(B);
#endif
  }    
    else { //not inplace

    Type1 *A;
    Type2 *B;
    int size1=(istride*N+idist*m);
    int size2=(ostride*N+odist*m);
    //    printf("size1=%d,size2=%d\n",size1,size2);    
#ifdef DEBUG
    printf("%d: out-of-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",grid1->taskid,istride,ostride,idist,odist);
#endif

#ifdef FFTW
    A = (Type1 *) fftw_malloc(size1*sizeof(Type1));
    //    A = (Type1 *) fftw_malloc(size1);
    B = (Type2 *) fftw_malloc(size2*sizeof(Type2));
    //    A = new Type1[size1];
    //B = new Type2[size2];

#ifdef DEBUG
    cout << "Allocated A and B. Types:" << type-> dt1 << " and " << type->dt2 << endl;
#endif
    //    A = (Type1 *) fftw_malloc(sizeof(Type1)*size1);
    //B = (Type2 *) fftw_malloc(sizeof(Type2)*size2);
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      for(int i=0;i < size1;i++) {
	A[i] = 0.0;
	B[i] = 1.0;
      }       
      //      if(isign == 0) 
      //	cout << "Error in find_plan: isign is not set" << endl;
#ifdef DEBUG
      printf("Calling doplan\n");
#endif
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);

#ifdef DEBUG
      printf("%d: Plan created %ld\n",grid1->taskid,plan->libplan);
#endif
    }
    else { //R2C or C2R
#ifdef DEBUG
      cout << "Calling doplan" << endl;
#endif
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,(double *) A,NULL,istride,idist,(fftw_complex *) B,NULL,ostride,odist,fft_flag);
#ifdef DEBUG
      printf("%d: Plan created %ld\n",grid1->taskid,plan->libplan);
#endif
    }
    //    delete [] A; // fftw_free
    //delete [] B;
    fftw_free(A);
    fftw_free(B);
#else
    //    A = (Type1 *) malloc(sizeof(Type1)*size1);
    //B = (Type2 *) malloc(sizeof(Type2)*size2);
    A = new Type1[size1];
    B = new Type2[size2];
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      //      if(isign == 0) 
      //	cout << "Error in find_plan: isign is not set" << endl;
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);
    }
    else //R2C or C2R
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    free(A);
    free(B);
#endif
  }
  
    return(plan->libplan);
}





int print_type3D(const trans_type3D *type)
{
  printf("trans_type3D values:\n");
  cout << "dt1,dt2=" << type->dt1 << type->dt2 << "prec=" << type->prec << "types=" << type->types[0] << type->types[1] << type->types[2] << "is_set=" << type->is_set << endl;
 
}

}
