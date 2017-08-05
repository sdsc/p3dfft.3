#include "p3dfft.h"

template<class Type1,class Type2> transform3D<Type1,Type2>::transform3D(const grid& grid1_, const grid& grid2_,const trans_type3D *type,bool inplace_)
{

  //#ifdef DEBUG
  cout << "In transform3D" << endl;
  print_type3D(type);
  //#endif

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

  int pgrid1[3],pgrid2[3],nstages,pgrid[3],proc_order[3],gdims[3];
  int ns,d1,d2,nd,dt1,dt2,i,df,L1,L2,L,L3,excl(int,int),dist(int);
  void swap0(int new_mo[3],int mo[3],int L);
  int monext[3];
  int dims1[3],dims2[3];
  gen_trans_type *tmptype;
  MPI_Comm mpicomm;
  bool reverse_steps;

  /*
  if(!grid1.is_set || !grid2.is_set) {
    printf("Error in tran3D_plan: grid is not set up\n");
    return;
  }
  */

  if((nd = grid1_.nd) != grid2_.nd) {
    printf("ERror in tran3D_plan: dimensions of grids don't match %d %d\n",nd,grid2_.nd);
    return;
  }
  for(i=0; i < nd; i++)
    if(grid1_.P[i] != grid2_.P[i] || grid1_.proc_order[i] != grid2_.proc_order[i]) {
      cout << "Error in transform3D: processor grids dont match" <<endl;
      return;
    }


    int tmp;
    MPI_Comm_size(grid1_.mpicomm[0],&tmp);
    printf("size of grid1_ mpicomm=%d\n",tmp);
  grid1 = new grid(grid1_);
    MPI_Comm_size(grid1->mpicomm[0],&tmp);
    printf("size of grid1 mpicomm=%d\n",tmp);
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

  if(nd == 1) {
    d1 = grid1_.D[0];
    if(d1 < 2) d2 = d1+1;
    else d2=d1 - 1;
    
    //  ! First plan local transforms in the 1st local dimensions
    L1 = grid1_.L[0];
    swap0(monext,mocurr,L1);
    tmpgrid0->set_mo(monext);

    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      gdims[i] = grid1_.gdims[i];
      proc_order[i] = grid1_.proc_order[i];
    }
    L2 = grid1_.L[1];
    swap0(monext,mocurr,L2);
    
    tmptype = types1D[type->types[L1]];
    //    grid tmpgrid1 = newgrid(tmpgrid0,tmptype,d);
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[d1] = gdims[d1]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[d1] = (gdims[d1]-1)*2;
    tmpgrid1 = new grid(gdims,pgrid1,proc_order,monext,mpicomm);
      
  //    grid tmpgrid1 = grid(tmpgrid0);
      printf("Calling init_transplan, trans_dim=%d\n",L1);
    Stages = curr_stage = init_transplan(*tmpgrid0,*tmpgrid1,tmptype,L1,inpl,prec);
    prev_stage = curr_stage;
    dt_prev = tmptype->dt2;
    
    //Plan transform of second dimension followed by an exchange of local dimensions with distributed dim. 
    delete tmpgrid0;
    tmptype = types1D[type->types[L2]];
    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      pgrid[i] = tmpgrid1->pgrid[i];
      gdims[i] = tmpgrid1->gdims[i];
      proc_order[i] = tmpgrid1->proc_order[i];
    }


    L1 = grid1_.L[1];
    L2 = d1;
    swap0(monext,mocurr,L2);
    pgrid[d1] = 1;
    pgrid[d2] = pgrid1[d1];
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[d1] = gdims[d1]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[d1] = (gdims[d1]-1)*2;

    tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

      printf("Calling init_trans_MPIplan, trans_dim=%d\n",L1);
    curr_stage = init_trans_MPIplan(*tmpgrid1,*tmpgrid0,grid1_.mpicomm[0],d1,d2,tmptype,L1,inpl,prec);
    curr_stage->kind = TRANSMPI;
    //,dt_prev,dt2
    prev_stage->next = curr_stage;
    //    dt_prev = type.types[0]->dt2;

    //  curr_stage->MPI_plan = stage::MPIplan(tmpgrid0,tmpgrid1,grid1.mpicomm[0],d1,d2,dt2);

    // Plan transform of final dimension
    delete tmpgrid1;
    tmptype = types1D[type->types[d1]];
    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      pgrid[i] = tmpgrid0->pgrid[i];
      gdims[i] = tmpgrid0->gdims[i];
      proc_order[i] = tmpgrid0->proc_order[i];
    }

    L1 = d1;
    swap0(monext,mo2,L1);
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[d2] = gdims[d2]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[d2] = (gdims[d2]-1)*2;

    tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

    prev_stage = curr_stage;
      printf("Calling init_transplan, trans_dim=%d\n",d1);
    curr_stage = init_transplan(*tmpgrid0,*tmpgrid1,tmptype,d1,inpl,prec);
    prev_stage->next = curr_stage;
    dt_prev = tmptype->dt2;

    //If needed, transpose back to the desired layout, specified by grid2
    
    d1 = d2;
    d2 = grid2_.D[0];
    if(d1 != d2) {
      prev_stage = curr_stage;
      curr_stage = init_MPIplan(*tmpgrid1,grid2_,grid1_.mpicomm[0],d1,d2,dt_prev,prec);
      prev_stage->next = curr_stage;
    }
  }
  else if(nd == 2) { //2D decomposition

    //Plan transform of first dimension followed by an exchange of local dimensions with first distributed dim. (row)

    // Determine the first, final and intermediate local dimensions
    d2 = L1 = tmpgrid0->L[0];
    L3 = grid2_.L[0];
    if(L1 == L3) {
      reverse_steps=true;
      L3 = dist(L1);
    }
    L2 = excl(L1,L3);

    swap0(monext,mocurr,L1);
    tmpgrid0->set_mo(monext);
    d1 = L2;
    /*
    switch(L1) {
    case 0: 
    case 2:
      d1 = 1;
      break;
    case 1: d1 = tmpgrid0->D[0];
      break;
    default:
    }
    */

    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      gdims[i] = grid1_.gdims[i];
      proc_order[i] = grid1_.proc_order[i];
      pgrid[i] = grid1_.pgrid[i];
    }
    //Exchange rows with local dimension
    pgrid[d1] = 1;
    pgrid[d2] = grid1_.pgrid[d1];  ////P[0];  // proc_order[0] ???
    //    tmpgrid1->L[0] = d1;
    //tmpgrid1->D[0] = d2;
    swap0(monext,mocurr,d1);
    
    tmptype = types1D[type->types[L1]];
    //    grid tmpgrid1 = newgrid(tmpgrid0,tmptype,d);
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[L1] = gdims[L1]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[L1] = (gdims[L1]-1)*2;
    tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

    //    printf("New grid, D=%d %d\n",tmpgrid1->D[0],tmpgrid1->D[1]);

    //    for(i=0; i < 3; i++) {
      //      tmpgrid1->ldims[i] = tmpgrid1->gdims[i]/tmpgrid1->pgrid[i];
    // mocurr[i] = monext[i];
    // }
    // tmpgrid1->set_mo(monext);
    printf("Calling init_tran_MPIsplan, trans_dim=%d, gdims2=(%d %d %d), ldims2=(%d %d %d), mem_order=(%d %d %d)\n",L1,gdims[0],gdims[1],gdims[2],tmpgrid1->ldims[0],tmpgrid1->ldims[1],tmpgrid1->ldims[2],monext[0],monext[1],monext[2]);

    Stages = curr_stage = init_trans_MPIplan(*tmpgrid0,*tmpgrid1,grid1_.mpicomm[0],d1,d2,tmptype,L1,inpl,prec);
    curr_stage->kind = TRANSMPI;
    //    curr_stage->MPI_plan = stage::MPIplan(tmpgrid0,tmpgrid1,grid1.mpicomm[0],d1,d2,dt_prev);

    //Plan transform of second dimension followed by an exchange of local dimensions with second distributed dim. (column)
    L = d2 = d1;
    d1 = L3;
    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
    }
    swap0(monext,mocurr,d1);

    delete tmpgrid0;
    //Exchange columns with local dimension
    pgrid[d1] = 1;
    pgrid[d2] = grid1_.P[1];
    //tmpgrid2.L[0] = d1;
    //tmpgrid2.D[1] = d2;
    tmptype = types1D[type->types[L]];
    //    grid tmpgrid1 = newgrid(tmpgrid0,tmptype,d);
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[L] = gdims[L]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[L] = (gdims[L]-1)*2;
    tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

    prev_stage = curr_stage;
      printf("Calling init_trans_MPIplan, trans_dim=%d, gdims2=(%d %d %d),ldims2=(%d %d %d), mem_order=(%d %d %d)\n",L,gdims[0],gdims[1],gdims[2],tmpgrid0->ldims[0],tmpgrid0->ldims[1],tmpgrid0->ldims[2],monext[0],monext[1],monext[2]);
    curr_stage = init_trans_MPIplan(*tmpgrid1,*tmpgrid0,grid1_.mpicomm[1],d1,d2,tmptype,L,inpl,prec);
    curr_stage->kind = TRANSMPI;

    prev_stage->next = curr_stage;

    //    curr_stage->MPI_plan = stage::MPIplan(tmpgrid1,tmpgrid2,grid1.mpicomm[1],d1,d2,dt_prev);

    // 1D transform the third dimension (which is now local)
    //    tmpgrid1 = grid(tmpgrid2);
    delete tmpgrid1;
    for(i=0; i < 3; i++) 
      mocurr[i] = monext[i];
    swap0(monext,mocurr,d1);
    //    tmpgrid1->set_mo(monext);
    tmptype = types1D[type->types[d1]];
    //    grid tmpgrid1 = newgrid(tmpgrid0,tmptype,d);
    if(tmptype->dt1 < tmptype->dt2) // Real-to-complex
      gdims[d1] = gdims[d1]/2+1;
    else if(tmptype->dt2 < tmptype->dt1) // Complex-to-real
      gdims[d1] = (gdims[d1]-1)*2;
    tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

    prev_stage = curr_stage;
      printf("Calling init_transplan, trans_dim=%d, gdims2=(%d %d %d),ldims2=(%d %d %d), mem_order=(%d %d %d)\n",d1,gdims[0],gdims[1],gdims[2],tmpgrid1->ldims[0],tmpgrid1->ldims[1],tmpgrid1->ldims[2],monext[0],monext[1],monext[2]);
    curr_stage = init_transplan(*tmpgrid0,*tmpgrid1,tmptype,d1,inpl,prec);
    
    prev_stage->next = curr_stage;
    dt_prev = tmptype->dt2;
    delete tmpgrid0;

    for(i=0;i<3;i++)
      if(tmpgrid1->gdims[i] != grid2_.gdims[i]) {
	printf("Error in transplan3D: global dimensions of final grid don't match. (%d %d %d) vs. (%d %d %d)\n",tmpgrid1->gdims[0],tmpgrid1->gdims[1],tmpgrid1->gdims[2],grid2_.gdims[0],grid2_.gdims[1],grid2_.gdims[2]);
	return;
      }

    if(reverse_steps) {

    d1 = tmpgrid1->D[1];
    df = grid2_.D[1];

    cout << "Return steps" << endl;
    if(d1 != df) {
      // Exchange D1 with L0
      d2 = L1 = tmpgrid1->L[0];
      pgrid[d1] = 1;
      pgrid[d2] = tmpgrid1->P[1];

      tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

      prev_stage = curr_stage;
      curr_stage = init_MPIplan(*tmpgrid1,*tmpgrid0,tmpgrid1->mpicomm[1],d1,d2,dt_prev,prec);
      prev_stage->next = curr_stage;
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
	pgrid[d1] = 1;
	pgrid[d2] = tmpgrid0->P[0];
	
	tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);
	
	prev_stage = curr_stage;
	curr_stage = init_MPIplan(*tmpgrid0,*tmpgrid1,tmpgrid0->mpicomm[0],d1,d2,dt_prev,prec);
	prev_stage->next = curr_stage;
	delete tmpgrid0;	
    }
      
    d1 = tmpgrid1->D[1];
    if(d1 != df) {
      // Again exchange D1 with L0
      d2 = L1 = tmpgrid1->L[0];
      pgrid[d1] = 1;
      pgrid[d2] = tmpgrid1->P[1];

      tmpgrid0 = new grid(gdims,pgrid,proc_order,monext,mpicomm);

      prev_stage = curr_stage;
      curr_stage = init_MPIplan(*tmpgrid1,*tmpgrid0,tmpgrid1->mpicomm[1],d1,d2,dt_prev,prec);
      prev_stage->next = curr_stage;
      delete tmpgrid1;
      L1 = d1;
      if(L1 != L2) {
	// Exchange L0 with D0
	d2 = tmpgrid0->L[0];
	d1 = tmpgrid0->D[0];
	pgrid[d1] = 1;
	pgrid[d2] = tmpgrid0->P[0];
	
	tmpgrid1 = new grid(gdims,pgrid,proc_order,monext,mpicomm);
	
	prev_stage = curr_stage;
	curr_stage = init_MPIplan(*tmpgrid0,*tmpgrid1,tmpgrid0->mpicomm[0],d1,d2,dt_prev,prec);
	prev_stage->next = curr_stage;
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
      printf("Calling init_transplan, trans_dim=%d\n",L2);
      curr_stage = init_transplan(*tmpgrid1,*tmpgrid0,types1D[EMPTY_TYPE],L2,inpl,prec);
      prev_stage->next = curr_stage;
      delete tmpgrid0;
    }
    delete tmpgrid1;

  }
  else if(nd == 3) {
    cout << "Three-dimensional decomposition is presently not supported" << endl;
    return;
  }

  cout << "Done transform3D planning" << endl;
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
    cout << "Error in trans_plan: 1D transform type no set" << endl;
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
  inembed = onembed = (int*) &gr1.gdims[d];
  if(type->dt1 < type->dt2) { //Real to complex
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = N;
    odist = N/2+1;
    if(odist <= gr2.gdims[d]) 
      onembed = (int*) &gr2.gdims[d];
    else {
      printf("Error in transplan: dimension too small %d, N=%d\n",gr2.gdims[d],N);
      return;
    }
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
    if(idist <= gr1.gdims[d]) 
      inembed = (int*) &gr1.gdims[d];
    else  {
      printf("Error in transplan: dimension too small %d, N=%d\n",gr2.gdims[d],N);
	return;
    }
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
    if(odist <= gr2.gdims[d]) 
      onembed = (int*) &gr2.gdims[d];
    else {
      printf("Error in transplan: dimension too small %d, N=%d\n",dims2[d],N);
      return;
    }
  }
  
  int mc[3];    //,rel_change(int [3],int [3],int [3]);
  rel_change(mo1,mo2,mc);
  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2                                                          
      m /= dims1[2]; // Need finer grain transform, to reuse cache     
      break;
    case 1:
      m = 1;
      break;
    }
    break;
    
  case 0:
    switch(mc[1]) {                                                             
    case 2: //2,1,0     
      m = 1;
      break;
    }
  }

	
  //  lib_plan = find_plan(trans_type);

}

template <class Type1,class Type2> transplan<Type1,Type2>::transplan(const grid &gr1,const grid &gr2,int type_ID,int d, bool inplace_) 
{

  is_empty = false;
  lib_plan = 0;plan = NULL;fft_flag = DEF_FFT_FLAGS;
  if(gr1.ldims[d] != gr1.gdims[d] || gr2.ldims[d] != gr2.gdims[d] ) {
    printf("Error in transplan: dimensions dont match %d %d %d\n",gr1.ldims[d],gr2.ldims[d],d);
    return;
  }

  trans_type = types1D[type_ID];

  if(!trans_type || !trans_type->is_set) {
    cout << "Error in trans_plan: 1D transform type no set" << endl;
    return;
  }

  double *A = new double[1024];

  dt1 = trans_type->dt1;
  dt2 = trans_type->dt2;
  stage_prec = prec = trans_type->prec;
  kind = TRANS_ONLY;
  trans_dim = d;
  inplace = inplace_;
  grid1 = new grid(gr1);
  grid2 = new grid(gr2);

  double *B = new double[1024];

  istride = 1;ostride = 1; 
  //  idist,odist;
  isign = trans_type->isign;

  for(int i=0;i<3;i++) {
    dims1[i] = gr1.ldims[i];
    dims2[i] = gr2.ldims[i];
    mo1[i] = gr1.mem_order[i];
    mo2[i] = gr2.mem_order[i];
  }
  inembed = onembed = (int*) &dims1[d];
  if(trans_type->dt1 < trans_type->dt2) { //Real to complex
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = N;
    odist = N/2+1;
    if(odist <= gr2.gdims[d]) 
      onembed = (int*) &(grid2.gdims[d]);
    else {
      cout << "Error in transplan: dimension too small" << dims2[d] << endl;
      return;
    }
  }
  else if(trans_type->dt1 > trans_type->dt2) { //Complex to real
    N=gr2.gdims[d];
    m=dims2[0]*dims2[1]*dims2[2]/N;
    odist = N;
    idist = N/2+1;
    if(idist <= gr1.gdims[d]) 
      inembed = (int*) &grid1.gdims[d];
    else  {
      cout << "Error in transplan: dimension too small" << gr1.gdims[d] << endl;
	return;
    }
  }
  else { // No change in datatype
    N=gr1.gdims[d];
    m=dims1[0]*dims1[1]*dims1[2]/N;
    idist = odist = N;
    if(odist <= gr2.gdims[d]) 
      onembed = (int*) &grid2.gdims[d];
    else {
      cout << "Error in transplan: dimension too small" << gr2.gdims[d] << endl;
      return;
    }
  }

  int mc[3];
  rel_change(mo1,mo2,mc);
  switch(mc1[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2                                                          
      m /= dims1[2]; // Need finer grain transform, to reuse cache     
      break;
    case 1:
      m = 1;
      break;
    }
    break;
    
  case 0:
    switch(mc[1]) {                                                             
    case 2: //2,1,0     
      m = 1;
      break;
    }
  }
  double *C = new double[1024];

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

  printf("find_plan: N,m=%d,%d\n",N,m);

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


  cout << "new plantype" << endl;
 plan = new Plantype<Type1,Type2>(type->doplan,type->exec,N,m,inplace,istride,idist,ostride,odist,inembed,onembed,isign,fft_flag);
  Plans.push_back(plan);
  
  if(inplace) {
    Type1 *A;
    int size=max(sizeof(Type1)*(istride*N+idist*m),sizeof(Type2)*(ostride*N+odist*m));
    printf("in-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",istride,ostride,idist,odist);
#ifdef FFTW
    A = (Type1 *) fftw_malloc(size *sizeof(Type1));
    Type2 *B = (Type2 *) fftw_malloc(size *sizeof(Type2));
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,inembed,istride,idist,B,onembed,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      if(isign == 0) 
	cout << "Error in find_plan: isign is not set" << endl;
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
    printf("%d: out-of-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",grid1->taskid,istride,ostride,idist,odist);
#ifdef FFTW
    A = (Type1 *) fftw_malloc(size1*sizeof(Type1));
    //    A = (Type1 *) fftw_malloc(size1);
    B = (Type2 *) fftw_malloc(size2*sizeof(Type2));
    //    A = new Type1[size1];
    //B = new Type2[size2];

    cout << "Allocated A and B. Types:" << type-> dt1 << " and " << type->dt2 << endl;
    //    A = (Type1 *) fftw_malloc(sizeof(Type1)*size1);
    //B = (Type2 *) fftw_malloc(sizeof(Type2)*size2);
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      for(int i=0;i < size1;i++) {
	A[i] = 0.0;
	B[i] = 1.0;
      }       
      if(isign == 0) 
	cout << "Error in find_plan: isign is not set" << endl;
      printf("Calling doplan\n");
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);

      printf("%d: Plan created %ld\n",grid1->taskid,plan->libplan);
    }
    else { //R2C or C2R
      cout << "Calling doplan" << endl;
      plan->libplan = (long) (*(plan->doplan))(1,&N,m,(double *) A,NULL,istride,idist,(fftw_complex *) B,NULL,ostride,odist,fft_flag);
      printf("%d: Plan created %ld\n",grid1->taskid,plan->libplan);
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
      if(isign == 0) 
	cout << "Error in find_plan: isign is not set" << endl;
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


template class transform3D<float,float>;
template class transform3D<double,double>;
template class transform3D<mycomplex,float>;
template class transform3D<complex_double,double>;
template class transform3D<float,mycomplex>;
template class transform3D<double,complex_double>;
template class transform3D<mycomplex,mycomplex>;
template class transform3D<complex_double,complex_double>;

int print_type3D(const trans_type3D *type)
{
  printf("trans_type3D %s values:\n",type->name);
  cout << "dt1,dt2=" << type->dt1 << type->dt2 << "prec=" << type->prec << "types=" << type->types[0] << type->types[1] << type->types[2] << "is_set=" << type->is_set << endl;
 
}
