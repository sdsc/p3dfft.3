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

void inv_mo(int mo[3],int imo[3]);


  // Set up 3D transform, defined by grid1 and grid2 (before and after grid configurations).
  // 
  template<class Type1,class Type2> transform3D<Type1,Type2>::transform3D(const DataGrid& grid1_, const DataGrid& grid2_,const trans_type3D *type)
{


#ifdef DEBUG
  cout << "In transform3D" << endl;
  print_type3D(type);
#endif

  int prec2;
  Stages = NULL;is_set = false;
  if(typeid(Type1) == type_float) {
    prec = 4;
  }
  else if(typeid(Type1) == type_double) {
    prec = 8;
  }
  else if(typeid(Type1) == type_complex) {
    prec = 4;
  }
  else if(typeid(Type1) == type_complex_double) {
    prec = 8;
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

  int dt = dt1 = sizeof(Type1)/prec;
  dt2 = sizeof(Type2)/prec;

  int dmap1[3],dmap2[3],gdims[3],L[3],st;
  int d1,d2,nd,i;
  void swap0(int new_mo[3],int mo[3],int L);
  int monext[3];
  int dims1[3],dims2[3];
  gen_trans_type *tmptype;
  bool reverse_steps;
  bool orig_input=true;

  if(!(*grid1_.Pgrid == *grid2_.Pgrid)) {
    printf("Error in transform3D: processor grids dont match\n"); 
    MPI_Abort(grid1_.Pgrid->mpi_comm_glob,0);
  }

  nd = grid1_.nd;
  grid1 = new DataGrid(grid1_);
  grid2 = new DataGrid(grid2_);
  Pgrid = new ProcGrid(*grid1_.Pgrid);

  memcpy(dmap1,grid1_.Dmap,sizeof(int)*3);
  memcpy(dmap2,grid2_.Dmap,sizeof(int)*3);

  stage *prev_stage,*curr_stage;
  int dt_prev = dt;
  int mocurr[3],mo1[3],mo2[3],tmp;
  for(i=0; i < 3; i++) {
    mocurr[i] = mo1[i] = grid1_.MemOrder[i];
    mo2[i] = grid2_.MemOrder[i];
  }


  /*
  if(nd == 3) { // Special treament for 3D decomposition: first transpose to get one local dimension
    
    // Find dimension corresponding to rows (adjacent tasks)
    for(i=0;i<3;i++)
      if(proc_order[i] == 0)
	break;
    grid1->nd = 2;
    grid1->L[0] = i;
    grid1->P[0] = grid1->pdims[i] = 1;
    grid1->D[0] = grid1->D[1];
    grid1->D[1] = grid1->D[2];
    grid1->D[2] = -1;
    //    splitcomm = 0;

    curr_stage = init_MPIplan(grid1_,*grid1,d1,d2,dt_prev,prec);
    Stages = prev_stage = curr_stage;
    curr_stage->kind = MPI_ONLY;

  } // nd==3
  */

  DataGrid *tmpgrid0 = new DataGrid(*grid1);
  DataGrid *tmpgrid1,*intgrid;
  int use_type2,dmap[3],dt_1,dt_2;

  /* Find the order of the three transforms, attempting to minimize reordering and transposes */ 

  bool init_steps = find_order(L,type, &grid1_, &grid2_, &reverse_steps);

  if(dt1 != types1D[type->types[L[0]]]->dt1)
    printf("Error in transform3D: input datatypes don't match\n");
  if(dt2 != types1D[type->types[L[2]]]->dt2)
    printf("Error in transform3D: output datatypes don't match\n");

  for(i=0;i<3;i++) 
    monext[i] = mocurr[i];

#ifdef DEBUG
    printf("%d: Planning stages: %d %d %d\n",Pgrid->taskid,L[0],L[1],L[2]);
#endif

    if(init_steps) {
      if(grid1_.Pdims[L[0]] != 1) { // Plan a transpose
	for(i=0;i<3;i++) {
	  gdims[i] = tmpgrid0->Gdims[i];
	  dmap[i] = tmpgrid0->Dmap[i];
	  //	  pdims[i] = tmpgrid0->Pdims[i];
	}
	d2 = L[2];
	d1 = L[0];
	// Swap d1 with d2
	tmp = dmap[d1];
	dmap[d1] = dmap[d2];
	dmap[d2] = tmp;

	tmpgrid1 = new DataGrid(gdims,grid1_.dim_conj_sym,Pgrid,dmap,monext);

	curr_stage = dynamic_cast<stage*> (new MPIplan<Type1>(*tmpgrid0,*tmpgrid1,d1,d2,prec));
	curr_stage->kind = MPI_ONLY;
	Stages = prev_stage = curr_stage;
	delete tmpgrid0;
	tmpgrid0 = tmpgrid1;
      }

    }

  // Plan the stages of the algorithm
  for(int st=0;st < 3;st++) {
    
    for(i=0;i<3;i++) {
      mocurr[i] = monext[i];
      gdims[i] = tmpgrid0->Gdims[i];
      dmap[i] = tmpgrid0->Dmap[i];
      //      pdims[i] = tmpgrid0->Pdims[i];
    }

    // Determine if an MPI transpose is involved in this stage
    d1 = -1;

    int dim_conj_sym = tmpgrid0->dim_conj_sym;
    tmptype = types1D[type->types[L[st]]];
    if(tmptype->dt1 < tmptype->dt2) { // Real-to-complex
      dim_conj_sym = L[st];
      gdims[dim_conj_sym] = gdims[dim_conj_sym]/2+1;
    }
    else if(tmptype->dt2 < tmptype->dt1) { // Complex-to-real
      gdims[dim_conj_sym] = (gdims[dim_conj_sym]-1)*2;
      dim_conj_sym = -1;
    }

    if(st < 2) { 
      if(!tmptype->is_empty)
	swap0(monext,mocurr,L[st]);
      if(tmpgrid0->Pdims[L[st+1]] > 1) {
	d1 = L[st+1];
	d2 = L[st];
      }
   }
    else   {
      for(i=0;i<3;i++) 
	monext[i] = grid2_.MemOrder[i];
      if(monext[L[st]] != 0)
	swap0(monext,mocurr,L[st]);
    }

    dt_1 = tmptype->dt1;
    dt_2 = tmptype->dt2;
    if(dt_1 != dt_prev)
      printf("Error in transform3D: datatypes don't match %d %d\n",dt_1,dt_prev);

    if(dt_1 != dt_2) 
      use_type2 = 12;
    else if(dt_1 == dt1)
      use_type2 = 1;
    else
      use_type2 = 2;

    if(d1 >= 0) { // If a transpose is involved

      // Swap d1 and d2 dimensions
      tmp = dmap[d1];
      dmap[d1] = dmap[d2];
      dmap[d2] = tmp;

      // Set up the new DataGrid (ending grid for this stage)
      //      pdims[d1] = 1;
      //pdims[d2] = tmpgrid0->pdims[d1];  ////P[0];  // proc_order[0] ???

      tmpgrid1 = new DataGrid(gdims,grid1_.dim_conj_sym,Pgrid,dmap,monext);

      if(!tmptype->is_empty) { // Combine transpose with transform

	if(tmpgrid0->MemOrder[L[st]] == 0)
	  intgrid = new DataGrid(*tmpgrid0);
	else if(tmpgrid1->MemOrder[L[st]] == 0) {
	  intgrid = new DataGrid(*tmpgrid0);
	  intgrid->set_mo(tmpgrid1->MemOrder);
	}
	else 
	  printf("Error initializing trans_MPIplan: neither input nor output is lead-dimension local\n");
	
	if(dt_1 != dt_2) {
	  intgrid->get_gdims(gdims);
	  if(dt_1 == REAL) // R2C
	    gdims[L[st]] = gdims[L[st]]/2+1;
	  else
	    gdims[L[st]] = (gdims[L[st]]-1)*2;
	  intgrid->set_gdims(gdims);
	}
	
#ifdef DEBUG
	printf("Calling init_tran_MPIsplan, trans_dim=%d, d1=%d, d2=%d, gdims2=(%d %d %d), Ldims2=(%d %d %d), MemOrder=(%d %d %d)\n",L[st],d1,d2,gdims[0],gdims[1],gdims[2],tmpgrid1->Ldims[0],tmpgrid1->Ldims[1],tmpgrid1->Ldims[2],monext[0],monext[1],monext[2]);
#endif
      // Plan/set up 1D transform combined with MPI exchange for this stage
	switch(use_type2) {
	case 1:
	  curr_stage = dynamic_cast<stage*> (new trans_MPIplan<Type1,Type1>(*tmpgrid0,*intgrid,*tmpgrid1,d1,d2,tmptype,L[st]));
	  break;
	case 2:
	  curr_stage = dynamic_cast<stage*> (new trans_MPIplan<Type2,Type2>(*tmpgrid0,*intgrid,*tmpgrid1,d1,d2,tmptype,L[st]));
	  break;
	case 12:
	  curr_stage = dynamic_cast<stage*> (new trans_MPIplan<Type1,Type2>(*tmpgrid0,*intgrid,*tmpgrid1,d1,d2,tmptype,L[st]));
	  break;
	}

      //      curr_stage = init_trans_MPIplan(*tmpgrid0,*intgrid,*tmpgrid1,d1,d2,tmptype,L[st],prec);
	curr_stage->kind = TRANSMPI;
	delete intgrid;
      //      curr_stage->inplace = inpl1D;
      }
      else { // MPI only
#ifdef DEBUG
	printf("Calling init_MPIsplan, d1=%d, d2=%d, gdims2=(%d %d %d), Ldims2=(%d %d %d), MemOrder=(%d %d %d)\n",d1,d2,gdims[0],gdims[1],gdims[2],tmpgrid1->Ldims[0],tmpgrid1->Ldims[1],tmpgrid1->Ldims[2],monext[0],monext[1],monext[2]);
#endif
	if(dt_1 == dt1)
	  curr_stage = dynamic_cast<stage*> (new MPIplan<Type1>(*tmpgrid0,*tmpgrid1,d1,d2,prec));
	else
	  curr_stage = dynamic_cast<stage*> (new MPIplan<Type2>(*tmpgrid0,*tmpgrid1,d1,d2,prec));
	curr_stage->kind = MPI_ONLY;
      }
    }
    else { // Only transform

      // Set up the ending grid for this stage
      tmpgrid1 = new DataGrid(gdims,dim_conj_sym,Pgrid,dmap,monext);
#ifdef DEBUG
      printf("Calling init_transplan, trans_dim=%d, gdims2=(%d %d %d), Ldims2=(%d %d %d), dmap=(%d %d %d), dim-conj_sym=%d, MemOrder=(%d %d %d)\n",L[st],gdims[0],gdims[1],gdims[2],tmpgrid1->Ldims[0],tmpgrid1->Ldims[1],tmpgrid1->Ldims[2],dmap[0],dmap[1],dmap[2],dim_conj_sym,monext[0],monext[1],monext[2]);
#endif

      // Plan/set up 1D transform, possibly combined with local transpose as needed
      switch(use_type2) {
      case 1:
	curr_stage = dynamic_cast<stage*> (new transplan<Type1,Type1>(*tmpgrid0,*tmpgrid1,tmptype,L[st]));
	break;
      case 2:
	curr_stage = dynamic_cast<stage*> (new transplan<Type2,Type2>(*tmpgrid0,*tmpgrid1,tmptype,L[st]));
	break;
      case 12:
	curr_stage = dynamic_cast<stage*> (new transplan<Type1,Type2>(*tmpgrid0,*tmpgrid1,tmptype,L[st]));
	break;
      }

      //      curr_stage = init_transplan(*tmpgrid0,*tmpgrid1,tmptype,L[st],prec);

      curr_stage->kind = TRANS_ONLY;
            
      }
    
    dt_prev = tmptype->dt2;
    delete tmpgrid0;
    tmpgrid0 = tmpgrid1;
    if(Stages == NULL)
      Stages = prev_stage = curr_stage;
    else {
      prev_stage->next = curr_stage;
      prev_stage = curr_stage;
    }

  }

    //If needed, transpose back to the desired layout, specified by grid2

  if(!reverse_steps) 
    for(i=0;i<3;i++) 
      if(dmap[i] != grid2_.Dmap[i] && grid2_.Pdims[i] > 1)
	reverse_steps = true;

  if(reverse_steps) 
    curr_stage = final_trans<Type2>(&tmpgrid0,grid2_,curr_stage,prec);

  bool iseq = true;
  for(i=0; i < 3; i++) 
    if(monext[i] != grid2_.MemOrder[i]) {
      iseq = false;
      break;
    }
  
  if(!iseq) { //If not in the final memory ordering
    prev_stage = curr_stage;
    curr_stage = final_seq<Type2>(*tmpgrid0,grid2_);
    curr_stage->kind = TRANS_ONLY;
    prev_stage->next = curr_stage;
  }

  delete tmpgrid0;

  if(nd == 3) {
    cout << "Three-dimensional decomposition is presently not supported" << endl;
    return;
  }
  
  //  cout << "Done transform3D planning" << endl;
  is_set = true;
  
}

template <class Type> stage *final_seq(const DataGrid &grid1, const DataGrid &grid2)
{
  gen_trans_type *t=empty_type<Type>();
  int trans_dim = grid2.L[0];
  if(grid2.nd == 1 && grid1.MemOrder[trans_dim] != 0 && grid2.MemOrder[trans_dim] != 0)
    trans_dim = grid2.L[1];
#ifdef DEBUG
  printf("Calling final init_transplan, trans_dim=%d\n",grid2.L[0]);
#endif
  return(dynamic_cast<stage*> (new transplan<Type,Type>(grid1,grid2,t,trans_dim)));
      //    transplan<Type2,Type2> *tr  = new transplan<Type2,Type2>(*tmpgrid1,*tmpgrid0,t,trans_dim);
  //curr_stage = (stage *) tr;
  // = init_transplan(*tmpgrid1,*tmpgrid0,types1D[EMPTY_TYPE],L2,inpl,prec);
}

template <class Type> stage *final_trans(DataGrid **grid1, const DataGrid &grid2,stage *curr,int prec)
  {
#ifdef DEBUG
    cout << "Return steps" << endl;
#endif

    int L1,d1,d2,nd,tmp,dmap[3],i,*mo,*gdims;
    ProcGrid *Pgrid;
    DataGrid *tmpgrid0,*tmpgrid1;
    stage *prev;

    tmpgrid1 = *grid1;
    L1 = tmpgrid1->L[0];
    d1 = tmpgrid1->D[0];
    nd = tmpgrid1->nd;
    for(i=0;i<3;i++)
      dmap[i] = tmpgrid1->Dmap[i];
    mo = tmpgrid1->MemOrder;
    Pgrid = tmpgrid1->Pgrid;
    gdims = tmpgrid1->Gdims;


    if(nd == 2) {
      d2 = tmpgrid1->D[1];
      if(dmap[d1] == grid2.Dmap[L1]) {
	//Exchange d1 and L
	  tmp = dmap[d1];
	  dmap[d1] = dmap[L1];
	  dmap[L1] = tmp;
	  
#ifdef DEBUG
	  printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d1,L1);
#endif
	  tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	  
	  prev = curr;
	  curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d1,L1,prec));
	  prev->next = curr;
	  curr->kind = MPI_ONLY;
	  
	  tmp = L1;
	  L1 = d1;
	  d1 = tmp;

	  if(dmap[d2] != grid2.Dmap[d2]) {
	    //Exchange d2 and L
	    tmp = dmap[d2];
	    dmap[d2] = dmap[L1];
	    dmap[L1] = tmp;
	    
#ifdef DEBUG
	    printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d2,L1);
#endif
	    delete tmpgrid1;
	    tmpgrid1 = tmpgrid0;
	    tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	    
	    prev = curr;
	    curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d2,L1,prec));
	    prev->next = curr;
	    curr->kind = MPI_ONLY;
	  }
      }
    
    else

      if(dmap[d2] == grid2.Dmap[L1]) {
	//Exchange d2 and L
	tmp = dmap[d2];
	dmap[d2] = dmap[L1];
	dmap[L1] = tmp;
       
#ifdef DEBUG
	printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d2,L1);
#endif
	tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);

	prev = curr;
	curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d2,L1,prec));
	prev->next = curr;
	curr->kind = MPI_ONLY;

	tmp = L1;
	L1 = d2;
	d2 = tmp;
	
	if(dmap[d1] != grid2.Dmap[d1]) {
	  //Exchange d1 and L
	  tmp = dmap[d1];
	  dmap[d1] = dmap[L1];
	  dmap[L1] = tmp;
	  
#ifdef DEBUG
	    printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d1,L1);
#endif
	    delete tmpgrid1;
	    tmpgrid1 = tmpgrid0;
	    tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	    
	    prev = curr;
	    curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d1,L1,prec));
	    prev->next = curr;
	    curr->kind = MPI_ONLY;
	  }
      }
	
	else if(dmap[d1] == grid2.Dmap[d2]) {
	  //Exchange d1 and L
	  tmp = dmap[d1];
	  dmap[d1] = dmap[L1];
	  dmap[L1] = tmp;
	  
#ifdef DEBUG
	  printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d1,L1);
#endif
	  tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	  
	  prev = curr;
	  curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d1,L1,prec));
	  prev->next = curr;
	  curr->kind = MPI_ONLY;
	  delete tmpgrid1;
	  tmpgrid1 = tmpgrid0;
	  tmp = L1;
	  L1 = d1;
	  d1 = tmp;


	  //Exchange d2 and L
	  tmp = dmap[d2];
	  dmap[d2] = dmap[L1];
	  dmap[L1] = tmp;
	  
#ifdef DEBUG
	  printf("MPI plan 2 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d2,L1);
#endif
	  tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	  
	  prev = curr;
	  curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d2,L1,prec));
	  prev->next = curr;
	  curr->kind = MPI_ONLY;
	  delete tmpgrid1;
	  tmpgrid1 = tmpgrid0;
	  tmp = L1;
	  L1 = d2;
	  d2 = tmp;

	  if(dmap[d1] == grid2.Dmap[L1]) { // difficult case: 3 steps
	  //Exchange d1 and L
	    tmp = dmap[d1];
	    dmap[d1] = dmap[L1];
	    dmap[L1] = tmp;
	  
#ifdef DEBUG
	    printf("MPI plan 3 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d1,L1);
#endif
	    tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	  
	    prev = curr;
	    curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d1,L1,prec));
	    prev->next = curr;
	    curr->kind = MPI_ONLY;
	  }
	}
	}
	else if(nd == 1) {
	  d2 = grid2.D[0];

	  //Exchange d2 and d1
	  tmp = dmap[d2];
	  dmap[d2] = dmap[d1];
	  dmap[d1] = tmp;
	  
#ifdef DEBUG
	  printf("MPI plan 1 from dmap %d %d %d to %d %d %d, d1=%d,d2=%d\n",tmpgrid1->Dmap[0],tmpgrid1->Dmap[1],tmpgrid1->Dmap[2],dmap[0],dmap[1],dmap[2],d1,d2);
#endif
	  tmpgrid0 = new DataGrid(gdims,grid2.dim_conj_sym,Pgrid,dmap,mo);
	  
	  prev = curr;
	  curr = dynamic_cast<stage*> (new MPIplan<Type>(*tmpgrid1,*tmpgrid0,d1,d2,prec));
	  prev->next = curr;
	  curr->kind = MPI_ONLY;
	}
	
    delete tmpgrid1;
    *grid1 = tmpgrid0;
    return curr;
}

bool find_order(int L[3],const trans_type3D *tp,const DataGrid *gr1,const DataGrid *gr2,bool *reverse_steps)
{
  int i,d1,d2;
  bool init_steps=false;
  gen_trans_type *tmptype; 
  int excl(int,int),dist(int);
  *reverse_steps = false;

  for(i=0;i<3;i++) 
    L[i] = -1;

  for(i=0;i<3;i++) {
    tmptype = types1D[tp->types[i]];
  // First, see if we have a R2C 1D transform, and start there
    if(tmptype->dt1 < tmptype->dt2) { // Real-to-complex
      if(L[0] >= 0)
	printf("ERror in transform3D: more than one real-to-complex 1D transform\n");
      else
	L[0] = i;
    }
    // End on C2R
    else if(tmptype->dt1 > tmptype->dt2) { // Complex-to-real
      if(L[2] >= 0)
	printf("ERror in transform3D: more than one complex-to-real 1D transforms\n");
      else
	L[2] = i;
    }
  }

  if(L[0] >= 0 && L[2] >= 0)
    printf("Error in transform3D: can't have both R2C and C2R transforms\n");

// This is the local dimension as we start
  if(L[0] <0) { // If initial local dimension is not assigned yet
    if(L[2] < 0 || gr1->Pdims[L[2]] > 1) {
      //    if(L[2] != gr1.L[0]) {

      L[0] = gr1->L[0]; 
      if(gr1->nd == 1 && gr1->MemOrder[L[0]] > gr1->MemOrder[gr1->L[1]]) // Make sure the first local dimension is also leading dimension in storage order (if we have a choice, as in 1D case)
	//	if(gr1.MemOrder[gr1.L[1]] == 0)
	L[0] = gr1->L[1];
    }
    else if(gr1->nd==1) {
      L[0] = gr1->L[0];
      if(L[0] == L[2])
	L[0] = gr1->L[1];
    }
    else {  // 2D grid and the final dimension is local - have to start with distributed dimension, and transpose data an extra time
      //      L[0] = gr1.D[0];
      //L[0] = (L[2]+1)%3;
      int m=2;
      for(i=0;i<3;i++)
	if(i != L[2])
	  if(gr1->MemOrder[i] < m) {
	    m = gr1->MemOrder[i];
	    L[0] = i;
	  }
      if(m != 0)
	init_steps = true;
    }
  }


      // Next choose intermediate and final local dimensions, L[1] and L[2]
  if(L[2] < 0) {// Final dimension not assigned yet 
    if(gr1->nd == 1) { // Special case of 1D proc grid
      int d1=gr1->D[0];
      int d2=gr2->D[0];
      if(d1 == d2) *reverse_steps = true;
      if(d2 == L[0])
	L[1] = d1;
      else
	L[1] = excl(d1,L[0]);
    }
    else {  // nd=2
      if(gr1->D[0] == gr2->D[1])
	L[1] = gr1->D[0];
      else if(gr1->D[1] == gr2->D[0])
	L[1] = gr1->D[1];
      else { *reverse_steps = true;
	if(gr1->D[0] == gr2->D[0])
	  L[1] = gr1->D[1];
	else
	  L[1] = gr1->D[0];

      }
      if(gr1->L[0] == gr2->L[0])
	*reverse_steps = true;
      /*	

      if(gr1->MemOrder[L[0]] != 0) { //      lead=false;
	for(i=0;i<3;i++)
	  if(gr1->MemOrder[i] == 0) { // Choose next dimension as current stride-1
	    L[1] = i;
	    L[2] = excl(L[0],L[1]);
	    break;
	  }
      */
    }
    L[2] = excl(L[0],L[1]);
  }
  /*
  
  else {

       if(gr1->nd == 1) {
	 
	 if(L[0] == gr2->L[0])
	   L[2] = gr2->L[1];
	 else if(L[0] == gr2->L[1])
	   L[2] = gr2->L[0];
	 else {
	   L[1] = gr1->D[0];
	   L[2] = excl(L[0],L[1]);
	 }
  */
	   /*	   

	   int l1=gr2.MemOrder[gr2.L[0]];
	   int l2=gr2.MemOrder[gr2.L[1]];
	   if(l1 < l2)
	     L[2] = gr2.L[0];
	   else
	     L[2] = gr2.L[1];
	 }
	d1 = gr1.D[0];
	d2 = gr2.D[0];
	if(d1 == d2 && d1 > 1) 
	  *reverse_steps=true;

	if(L[0] == d2) {
	  L[1] = (L[0]+1)%3;
	  if(gr1.MemOrder[L[1]] == 2)
	    L[1] = (L[0]+2)%3;
	  Lfin = L[2] = excl(L[0],L[1]);
	}
	else 
	  if(gr2.MemOrder[d2] <= 1){ // In case of 1D we have a choice, so choose the one with most favorable transmutation  of indices in local reordering
	    L[1] = d2;
	    Lfin = L[2] = excl(L[0],L[1]);
	  }
	  else {
	    Lfin = L[2] = d2;
	    L[1] = excl(L[0],L[2]);
	    }*/
  /*
       }
      else { // nd=2
      
	L[2] = gr2->L[0]; // This dimension will be local upon finish
	if(L[2] == L[0]) {
	  *reverse_steps=true; // If start and finish local dimensions are the same, there is no way to complete execution of 3D transform in 3 steps with max. 2 exchanges. Additional MPI exchanges will be necessary to make the original local dimension local again in the end.
	  L[2] = dist(L[0]); // Choose target local dimension for the third step.
	}
	L[1] = excl(L[0],L[2]); // When first and third local dimensions are known, the second local dimension is found by exclusion.
      }
  

  */
  L[1] = excl(L[0],L[2]);

  if(gr1->Pdims[L[0]] != 1)
    init_steps = true;

  return init_steps;
}

void swap0(int newmo[3],int mo[3],int L)
{
  int i,j;
  for(i=0; i< 3; i++)
    newmo[i] = mo[i];

  if(newmo[L] == 0)
    return;
  if(newmo[L] == 1) {
    for(j=0;j<3;j++)
      if(newmo[j] == 0)
	break;
    newmo[j] = 1;
    newmo[L] = 0;
    return;
  }
  // [L]=2
  for(i=0;i<3;i++)
    newmo[i] = (newmo[i] +1)%3;

  /*
  if(newmo[L] !=0) {
    for(i=0; i < 3; i++)
      if(newmo[i] == 0)
	break;
    switch(abs(L-i)) {
    case 1:
      newmo[i] = L;
      break;
    case 2:
      if(L > i) {
	newmo[i] = mo[i+1];
	newmo[i+1] = mo[L];
      }
      else {
	newmo[i] = mo[i-1];
	newmo[i-1] = mo[L];
      }
    }
    newmo[L] = 0;
  }
  */
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


  template <class Type1,class Type2> transplan<Type1,Type2>::transplan(const DataGrid &gr1,const DataGrid &gr2, const gen_trans_type *type,int d) // bool inplace_) 
  {   init_tr(gr1,gr2, type,d);
  }

template <class Type1,class Type2> void transplan<Type1,Type2>::init_tr(const DataGrid &gr1,const DataGrid &gr2, const gen_trans_type *type,int d) // bool inplace_) 
{
  plan = NULL;fft_flag = DEF_FFT_FLAGS;
  if(!type->is_set) {
    cout << "Error in trans_plan: 1D transform type not set" << endl;
    return;
  }
  if(gr1.Ldims[d] != gr1.Gdims[d] || gr2.Ldims[d] != gr2.Gdims[d] ) {
    printf("Error in transplan: dimensions dont match %d, %d, %d\n",gr1.Ldims[d],gr2.Ldims[d],d);
    return;
  }
  stage_prec = prec = type->prec;
  trans_dim = d; 
  trans_type = (trans_type1D<Type1,Type2> *) type;


  dt1 = type->dt1;
  dt2 = type->dt2;
  kind = TRANS_ONLY;
  //  inplace = inplace_;
  grid1 = new DataGrid(gr1);
  grid2 = new DataGrid(gr2);

  istride = 1;ostride = 1; 
  //  idist,odist;
  isign = type->isign;

  for(int i=0;i<3;i++) {
    dims1[i] = grid1->Ldims[i];
    dims2[i] = grid2->Ldims[i];
    mo1[i] = grid1->MemOrder[i];
    mo2[i] = grid2->MemOrder[i];
  }

  if(trans_type->is_empty)
    is_empty = true;
  else {
    //  inembed = onembed = (int*) &grid1.Gdims[d];
    if(type->dt1 < type->dt2) { //Real to complex
      N=gr1.Gdims[d];
    //   m=dims1[0]*dims1[1]*dims1[2]/N;
      idist = N;
      odist = N/2+1;
      /*
    dims2[0] = dims2[0]/2+1;
    mygrid->Ldims[d] = (mygrid->Ldims[d]+2)/2;
    mygrid->Gdims[d] = (mygrid->Gdims[d]+2)/2;
    mygrid->sz[0][d] = (mygrid->sz[0][d]+2)/2;  
    mygrid->en[0][d] = (mygrid->en[0][d]+2)/2;  
    mygrid->st[0][d] = (mygrid->st[0][d]+2)/2; 
    */
  }
    else if(type->dt1 > type->dt2) { //Complex to real
      N=gr2.Gdims[d];
      //   m=dims2[0]*dims2[1]*dims2[2]/N;
      odist = N;
      idist = N/2+1;
    /*
    dims2[0] = dims2[0]*2-2;
    mygrid->Ldims[d] = mygrid->Ldims[d]*2-2;
    mygrid->Gdims[d] = mygrid->Gdims[d]*2-2;
    mygrid->sz[0][d] = mygrid->sz[0][d]*2-2;
    mygrid->st[0][d] = mygrid->st[0][d]*2-2;
    mygrid->en[0][d] = mygrid->en[0][d]*2-2;
    */
    }
    else { // No change in datatype
      N=gr1.Gdims[d];
      //    m=dims1[0]*dims1[1]*dims1[2]/N;
      idist = odist = N;
    }
    
    m = find_m(mo1,mo2,dims1,dims2,trans_dim);
    
    if(idist <= gr1.Gdims[d]) 
      inembed = (int*) &(grid1->Gdims[d]);
    else  {
      printf("Error in transplan: dimension too small %d, N=%d\n",gr1.Gdims[d],N);
      return;
    }
    if(odist <= gr2.Gdims[d]) 
      onembed = (int*) &(grid2->Gdims[d]);
    else {
      printf("Error in transplan: dimension too small %d, N=%d\n",dims2[d],N);
      return;
    }
    
    Pgrid = gr1.Pgrid;
    if(!(*Pgrid == *gr2.Pgrid))
      printf("Error in transplan: processor grids don't match\n");
    find_plan(trans_type); 
  }
}  

template <class Type1,class Type2> transplan<Type1,Type2>::transplan(const DataGrid &gr1,const DataGrid &gr2,int type_ID,int d) // bool inplace_) 
{   
  trans_type = (trans_type1D<Type1,Type2> *) types1D[type_ID];

  if(!trans_type || !trans_type->is_set) {
    cout << "Error in trans_plan: 1D transform type no set" << endl;
    return;
  }
  init_tr(gr1,gr2, trans_type,d);
}


#define TRANS_IN 0
#define TRANS_OUT 1


template <class Type1,class Type2> int transplan<Type1,Type2>::find_m(int *mo1,int *mo2,int *dims1,int *dims2, int trans_dim) {

  int i,m,mc[3],imo1[3],imo2[3],d1[3],d2[3], scheme;    //,rel_change(int [3],int [3],int [3]);
  void swap0(int new_mo[3],int mo[3],int L);

  if(mo1[trans_dim] == 0) 
    scheme = TRANS_IN;
  else if(mo2[trans_dim] == 0) 
    scheme = TRANS_OUT;
  else {  // Have to use a two-step scheme
    swap0(mo1,mo1,trans_dim); // Find an intermediate ordering with leading dimension we need 
    scheme = TRANS_IN;
  }

  inv_mo(mo1,imo1);
  inv_mo(mo2,imo2);
  for(i=0;i<3;i++) {
    d1[i] = dims1[imo1[i]];
    d2[i] = dims2[imo2[i]];
  }

  rel_change(imo1,imo2,mc);

  switch(mc[0]) {
  case 1:
    switch(mc[1]) {
    case 0: //1,0,2                                                          
      if(scheme == TRANS_IN) 
	m = d1[1]; // Need finer grain transform, to reuse cache     
      else
	m = d2[1];
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
      m = d2[1]*d2[2];
    else // 021
      m = 1;
    break;
    
  case 2:
    //    if(mc[1] == 0) { //2 0 1
      if(scheme == TRANS_IN) 
	m = d1[1]*d1[2];
      else
	m = d2[1]*d2[2];
      //    }
      //else //210
      //m = d2[1]*d2[2];
    break;
  }
  
  return(m);

}

  template <class Type1,class Type2> trans_MPIplan<Type1,Type2>::trans_MPIplan(const DataGrid &gr1,const DataGrid &intergrid, const DataGrid &gr2,int d1,int d2,const gen_trans_type *type,int d) //bool inplace_)   
{
  kind = TRANSMPI;
  trplan = new transplan<Type1,Type2>(gr1,intergrid,type,d); //,inplace_);
  if(trplan->plan == NULL && !trplan->is_empty)
    cout << "Error in trans_MPIplan: null plan" << endl;

  stage_prec = trplan->prec;
  mpiplan = new MPIplan<Type2>(intergrid,gr2,d1,d2,stage_prec);
  is_set = true;
  //  inplace = inplace_;
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

template <class Type> MPIplan<Type>::MPIplan(const DataGrid &gr1,const DataGrid &gr2, int d1_,int d2_, int prec_)  
{
  int i,d3,l;


  if(!gr1.is_set || !gr2.is_set) {
    printf("Error in MPIplan constr.: grids are not set up\n");
    return;
  }

  comm_id = gr1.Dmap[d1_]; // Dimension of processor grid where transpose takes place
  if(comm_id != gr2.Dmap[d2_])
    printf("Error in MPIplan: swap dimensions don't match\n");

  for(i=0; i < 3; i++){
    mo1[i] = gr1.MemOrder[i];
    mo2[i] = gr2.MemOrder[i];
    if(gr1.Gdims[i] != gr2.Gdims[i]) {
      printf("Error in MPIplan constr.: global grid dimensions dont match: (%d %d %d) vs. (%d %d %d)\n", gr1.Gdims[0],gr1.Gdims[1],gr1.Gdims[2],gr2.Gdims[0],gr2.Gdims[1],gr2.Gdims[2]);
      return;
    }
  }

  int p=gr1.Pdims[d1_];
  if(p != gr2.Pdims[d2_]) {
    cout << "Error in MPIplan constr.: proc. grid dimensions dont match" << gr1.Pdims[d1_] << ", " << gr2.Pdims[d2_] << d1_ << ", " << d2_ << endl;
    return;
  }

  stage_prec = prec = prec_;

  Pgrid = gr1.Pgrid;
  if(!(*Pgrid == *gr2.Pgrid)) 
    printf("Error in MPIplan: ProcGrids don't match\n");
  MPI_Comm mpicomm = Pgrid->mpicomm[comm_id];

  if(p >1) {

 
  MPI_Comm_size(mpicomm,&numtasks);
  if(p != numtasks) {
    cout << "Error in MPIplan constr.: proc. grid dimension doesnt match communicator size" << p << numtasks << endl;
    return;
  }

  MPI_Comm_rank(mpicomm,&taskid);
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

  /*
  for(l=0;l<gr1.nd;l++) 
    if(gr1.D[l] == gr1.proc_order[d1])
      break;
  */


  memcpy(dims1,gr1.Ldims,3*sizeof(int));
  memcpy(dims2,gr2.Ldims,3*sizeof(int));

  // int comm_coords[3];
  //memcpy(comm_coords,gr1.grid_id_cart,3*sizeof(int));

  SndStrt[0] = 0;
  RcvStrt[0] = 0;

  int sz=sizeof(Type)/4;
  
  //int rank;
  for(int j=0; j< p-1;j++) {
    //    comm_coords[l] = j;
    //MPI_Cart_rank(mpicomm_,comm_coords,&rank); 
    SndCnts[j] = gr2.sz[d2][j] * dims1[d3] *dims1[d1] * sz;
    RcvCnts[j] = dims2[d2] * gr1.sz[d1][j] *dims2[d3] * sz;
    SndStrt[j+1] = SndStrt[j] + SndCnts[j];
    RcvStrt[j+1] = RcvStrt[j] + RcvCnts[j];
  }
  //  comm_coords[l] = p-1;
  //MPI_Cart_rank(mpicomm_,comm_coords,&rank); 
  SndCnts[p-1] = gr2.sz[d2][p-1] * dims1[d3]*dims1[d1] * sz;
  RcvCnts[p-1] = dims2[d2] * gr1.sz[d1][p-1]*dims2[d3] * sz;

  grid1 = new DataGrid(gr1);
  grid2 = new DataGrid(gr2);
  du = d3;
  //  is_mpi = true;
  }

  is_set = true;
  //  is_trans = false;
  kind = MPI_ONLY;
  dt1 = dt2 = sizeof(Type)/prec;

}



  template <class Type1,class Type2> inline void transplan<Type1,Type2>::find_plan(trans_type1D<Type1,Type2> *type)
{
  int i;
  int planID;
  Plantype<Type1,Type2> *pl;
  //  Plan *p;

#ifdef DEBUG
  printf("find_plan: N,m=%d,%d\n",N,m);
#endif

  /*
  if(inplace) 
    if(sizeof(Type1) != sizeof(Type2) || istride != ostride || (m > 1 && idist != odist)) {
      cout << "Error in find_plan: inplace transforms should have identical dimensions and types" << endl;
      return(-1);
    }
  
  if(inplace)
    if(type->dt1 >= type->dt2 && !arcmp(mo1,mo2,3))
      inplace = true;
    else
      inplace=false;

  if(inplace)
    for(i=0;i<3;i++)
      if(mo1[i] != mo2[i]) {
	inplace = false;
	break;
      }
  */

  planID = 0;
  for(vector<Plan*>::iterator it=Plans.begin(); it < Plans.end();it++,planID++) {
    //    pl = dynamic_cast<Plantype<Type1,Type2> *> (p);
    //if(pl)
    pl = (Plantype<Type1,Type2> *) *it;
    if(type->dt1 == pl->dt1 && type->dt2 == pl->dt2 && prec == pl-> prec &&
       pl->N == N && pl->m == m && // pl->inplace == inplace &&		
       pl->istride == istride && pl->ostride == ostride && pl->isign == isign &&
       pl->fft_flag == fft_flag) {
      if(m > 1) {
	if(pl->idist == idist && pl->odist == odist &&
	   *(pl->inembed) == *inembed && *(pl->onembed) == *onembed) {
	  plan = pl;
#ifdef DEBUG
	  printf("find_plan: plan found, isign=%d %d\n",isign,plan->isign);
#endif 
	  return;
	  //plan = pl;
	  // return((Plantype<Type1,Type2> *) it);
	}
	//	  return(plan->libplan);
      }
      else {
	plan = pl;
#ifdef DEBUG
	printf("find_plan: plan found, isign=%d %d\n",isign,plan->isign);
#endif 
	return;
      }
    }
    //	return(plan->libplan);
  }

    // If haven't found existing plan with suitable params, define a new one
    //    plan = &DefPlans[i];


#ifdef DEBUG
  cout << "new plantype" << endl;
#endif
 plan = new Plantype<Type1,Type2>(type->doplan,type->exec,N,m,istride,idist,ostride,odist,inembed,onembed,isign,fft_flag);
  Plans.push_back(plan);
  
  // Inplace
    Type1 *A;
    int size=max(sizeof(Type1)*(istride*N+idist*m),sizeof(Type2)*(ostride*N+odist*m));
#ifdef DEBUG
    printf("in-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",istride,ostride,idist,odist);
#endif
#ifdef FFTW
    A = (Type1 *) fftw_malloc(size);
    Type2 *B = (Type2 *) A; //fftw_malloc(size *sizeof(Type2));
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan_in = (long) (*(plan->doplan))(1,&N,m,A,inembed,istride,idist,B,onembed,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      //      if(isign == 0) 
      //	cout << "Error in find_plan: isign is not set" << endl;
      plan->libplan_in = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);
    }
    else //R2C or C2R
      plan->libplan_in = (long) (*(plan->doplan))(1,&N,m,A,inembed,istride,idist,B,onembed,ostride,odist,fft_flag);

    if(plan->libplan_in == NULL) 
      printf("ERror: NULL plan in find_plan: N=%d,m=%d,istride=%d,idist=%d,ostride=%d,odist=%d,fft_flag=%d\n",N,m,istride,idist,ostride,odist,fft_flag);

    fftw_free(A);
    //    fftw_free(B);
#endif

    // Out of place
    int size1=(istride*N+idist*m);
    int size2=(ostride*N+odist*m);
    //    printf("size1=%d,size2=%d\n",size1,size2);    
#ifdef DEBUG
    printf("%d: out-of-place: istride=%d,ostride=%d,idist=%d,odist=%d\n",grid1->Pgrid->taskid,istride,ostride,idist,odist);
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
      plan->libplan_out = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
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
      plan->libplan_out = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);

#ifdef DEBUG
      printf("%d: Plan created %ld\n",grid1->Pgrid->taskid,plan->libplan_out);
#endif
    }
    else { //R2C or C2R
#ifdef DEBUG
      cout << "Calling doplan" << endl;
#endif
      plan->libplan_out = (long) (*(plan->doplan))(1,&N,m,(double *) A,NULL,istride,idist,(fftw_complex *) B,NULL,ostride,odist,fft_flag);
#ifdef DEBUG
      printf("%d: Plan created %ld\n",grid1->Pgrid->taskid,plan->libplan_out);
#endif
    }
    //    delete [] A; // fftw_free
    //delete [] B;

    if(plan->libplan_out == NULL) 
      printf("ERror: NULL plan in find_plan: N=%d,m=%d,istride=%d,idist=%d,ostride=%d,odist=%d,fft_flag=%d\n",N,m,istride,idist,ostride,odist,fft_flag);

    fftw_free(A);
    fftw_free(B);
#else // non-FFTW
    //    A = (Type1 *) malloc(sizeof(Type1)*size1);
    //B = (Type2 *) malloc(sizeof(Type2)*size2);
    A = new Type1[size1];
    B = new Type2[size2];
    if(type->dt1 == 1 && type->dt2 == 1) //Real-to-real
      plan->libplan_inout = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    else if(type->dt1 == 2 && type->dt2 == 2) { //Complex-to-complex
      //      if(isign == 0) 
      //	cout << "Error in find_plan: isign is not set" << endl;
      plan->libplan_inout = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,isign,fft_flag);
    }
    else //R2C or C2R
      plan->libplan_inout = (long) (*(plan->doplan))(1,&N,m,A,NULL,istride,idist,B,NULL,ostride,odist,fft_flag);
    free(A);
    free(B);
#endif
  
  return;
}





int print_type3D(const trans_type3D *type)
{
  printf("trans_type3D values:\n");
  cout  << "prec=" << type->prec << "types=" << type->types[0] << type->types[1] << type->types[2] << "is_set=" << type->is_set << endl;
 
}


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

