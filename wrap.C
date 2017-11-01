#include "p3dfft.h"

using namespace std;
using namespace p3dfft;

extern "C" {

void p3dfft_setup() {
    p3dfft::setup();

P3DFFT_EMPTY_TYPE=p3dfft::EMPTY_TYPE;
P3DFFT_R2CFFT_S=p3dfft::R2CFFT_S;
P3DFFT_R2CFFT_D=p3dfft::R2CFFT_D;
P3DFFT_C2RFFT_S=p3dfft::C2RFFT_S;
P3DFFT_C2RFFT_D=p3dfft::C2RFFT_D;
P3DFFT_CFFT_FORWARD_S=p3dfft::CFFT_FORWARD_S;
P3DFFT_CFFT_FORWARD_D=p3dfft::CFFT_FORWARD_D;
P3DFFT_CFFT_BACKWARD_S=p3dfft::CFFT_BACKWARD_S;
P3DFFT_CFFT_BACKWARD_D=p3dfft::CFFT_BACKWARD_D;
P3DFFT_COSTRAN_REAL_S=p3dfft::COSTRAN_REAL_S;
P3DFFT_COSTRAN_REAL_D=p3dfft::COSTRAN_REAL_D; 
P3DFFT_SINTRAN_REAL_S=p3dfft::SINTRAN_REAL_S;
P3DFFT_SINTRAN_REAL_D=p3dfft::SINTRAN_REAL_D;
}

void p3dfft_cleanup() {
    p3dfft::cleanup();
}

  Type3D p3dfft_init_3Dtype(int types[3]) //,char *name)
{
  trans_type3D tp = trans_type3D(types);
  int count = types3D.size();
  types3D.push_back(tp);
  return(count);
}

Plan3D p3dfft_plan_3Dtrans_f(int *Fgr1,int *Fgr2,Type3D *tp,int *inplace){

  /* 
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: type3D=%d.  initiating gr1\n",*tp);
#endif
 grid *gr1 = new grid(Fgr1->gdims,Fgr1->pgrid,Fgr1->proc_order,Fgr1->mem_order,MPI_Comm_f2c(Fgr1->mpi_comm_glob));
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: initiating gr2\n");
#endif
  */

  grid *gr1 = &stored_grids[*Fgr1];
  grid *gr2 = &stored_grids[*Fgr2];

  //  grid *gr2 = new grid(Fgr2->gdims,Fgr2->pgrid,Fgr2->proc_order,Fgr2->mem_order,MPI_Comm_f2c(Fgr2->mpi_comm_glob));
  trans_type3D *type3D = &types3D[*tp];
  gen_transform3D *tr3D;
    
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: new transform3D\n");
#endif

  if(type3D->prec == 4) 
    if(type3D->dt1 == 1)
      if(type3D->dt2 == 1)
	tr3D = new transform3D<float,float>(*gr1,*gr2,type3D,*inplace);
      else
	tr3D = new transform3D<float,mycomplex>(*gr1,*gr2,type3D,*inplace);
    else
      if(type3D->dt2 == 1)
	tr3D = new transform3D<mycomplex,float>(*gr1,*gr2,type3D,*inplace);
      else
	tr3D = new transform3D<mycomplex,mycomplex>(*gr1,*gr2,type3D,*inplace);
    else
    if(type3D->dt1 == 1)
      if(type3D->dt2 == 1)
	tr3D = new transform3D<double,double>(*gr1,*gr2,type3D,*inplace);
      else
	tr3D = new transform3D<double,complex_double>(*gr1,*gr2,type3D,*inplace);
    else
      if(type3D->dt2 == 1)
	tr3D = new transform3D<complex_double,double>(*gr1,*gr2,type3D,*inplace);
      else
	tr3D = new transform3D<complex_double,complex_double>(*gr1,*gr2,type3D,*inplace);

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: count\n");
#endif

  int count = stored_trans3D.size();

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: push back\n");
#endif

  stored_trans3D.push_back(tr3D);
  //  delete gr1,gr2;
  return count;
    
  /*
  gen_trans_type *tp1D[3];
  int i;

  for(i=0;i<3;i++) 
    tp1D[i] = types1D[type3D->types[i]];

  int dt1 = tp1D[0]->dt1;
  int dt2 = tp1D[2]->dt2;
  */

}
Plan3D p3dfft_plan_3Dtrans(Grid *Cgr1,Grid *Cgr2,Type3D tp,int inplace){

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: type3D=%d.  initiating gr1\n",tp);
#endif

  grid *gr1 = new grid(Cgr1->gdims,Cgr1->pgrid,Cgr1->proc_order,Cgr1->mem_order,Cgr1->mpi_comm_glob);
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: initiating gr2\n");
#endif

  grid *gr2 = new grid(Cgr2->gdims,Cgr2->pgrid,Cgr2->proc_order,Cgr2->mem_order,Cgr2->mpi_comm_glob);
  trans_type3D *type3D = &types3D[tp];
  gen_transform3D *tr3D;
    
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: new transform3D\n");
#endif

  if(type3D->prec == 4) 
    if(type3D->dt1 == 1)
      if(type3D->dt2 == 1)
	tr3D = new transform3D<float,float>(*gr1,*gr2,type3D,inplace);
      else
	tr3D = new transform3D<float,mycomplex>(*gr1,*gr2,type3D,inplace);
    else
      if(type3D->dt2 == 1)
	tr3D = new transform3D<mycomplex,float>(*gr1,*gr2,type3D,inplace);
      else
	tr3D = new transform3D<mycomplex,mycomplex>(*gr1,*gr2,type3D,inplace);
    else
    if(type3D->dt1 == 1)
      if(type3D->dt2 == 1)
	tr3D = new transform3D<double,double>(*gr1,*gr2,type3D,inplace);
      else
	tr3D = new transform3D<double,complex_double>(*gr1,*gr2,type3D,inplace);
    else
      if(type3D->dt2 == 1)
	tr3D = new transform3D<complex_double,double>(*gr1,*gr2,type3D,inplace);
      else
	tr3D = new transform3D<complex_double,complex_double>(*gr1,*gr2,type3D,inplace);

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: count\n");
#endif

  int count = stored_trans3D.size();

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: push back\n");
#endif

  stored_trans3D.push_back(tr3D);
  delete gr1,gr2;
  return count;
    
  /*
  gen_trans_type *tp1D[3];
  int i;

  for(i=0;i<3;i++) 
    tp1D[i] = types1D[type3D->types[i]];

  int dt1 = tp1D[0]->dt1;
  int dt2 = tp1D[2]->dt2;
  */

}

  int p3dfft_init_grid_f(int *ldims,int *glob_start,int *gdims,int *pgrid,int *proc_order,int *mem_order,int *mpicomm) {
    
    int num=find_grid(gdims,pgrid,proc_order,mem_order,MPI_Comm_f2c(*mpicomm));
    if(num >= 0) 
      return(num);
    else {

    grid *gr1;
  gr1 = new grid(gdims,pgrid,proc_order,mem_order,MPI_Comm_f2c(*mpicomm));
  memcpy(ldims,gr1->ldims,3*sizeof(int));
  memcpy(glob_start,gr1->glob_start,3*sizeof(int));
  num = stored_grids.size();
  stored_grids.push_back(*gr1);
  return(num);
    }
  /*
  memcpy(&gr->mem_order,mem_order,3*sizeof(int));
  memcpy(&gr->pgrid,pgrid,3*sizeof(int));
  memcpy(&gr->proc_order,proc_order,3*sizeof(int));
  memcpy(&gr->P,gr1->P,3*sizeof(int));
  memcpy(&gr->D,gr1->D,3*sizeof(int));
  memcpy(&gr->L,gr1->L,3*sizeof(int));
  memcpy(&gr->grid_id,gr1->grid_id,3*sizeof(int));
  gr->mpi_comm_cart = MPI_Comm_f2c(gr1->mpi_comm_cart);  
  for(int i=0;i<3;i++)
    gr->mpicomm[i] = MPI_Comm_c2f(gr1->mpicomm[i]);  
  */
  //  MPI_Comm_dup(MPI_Comm_f2c(gr1->mpi_comm_glob),&gr->mpi_comm_glob);  
  //  delete gr1;
}

  int find_grid(int gdims[3],int pgrid[3],int proc_order[3],int mem_order[3],MPI_Comm mpicomm) {
    
    int cnt=stored_grids.size();
    int i;
    int res;

    for(i=0;i<cnt;i++) {
      grid *gr=&stored_grids[i];
      if((arcmp(gr->gdims,gdims,3) == 0) && (arcmp(gr->pgrid,pgrid,3) == 0))
	if((arcmp(gr->proc_order,proc_order,3) == 0) && (arcmp(gr->mem_order,mem_order,3)== 0)) {
	MPI_Comm_compare(gr->mpi_comm_glob,mpicomm,&res);
	if(res == MPI_IDENT)
	  return(i);
      }
    }
    return(-1);
      
}

Grid *p3dfft_init_grid(int gdims[3],int pgrid[3],int proc_order[3],int mem_order[3],MPI_Comm mpicomm) {

  grid gr = grid(gdims,pgrid,proc_order,mem_order,mpicomm);
  Grid *Cgr = new Grid;
  Cgr->numtasks = gr.numtasks;
  Cgr->taskid = gr.taskid;
  Cgr->nd = gr.nd;
  memcpy(&Cgr->gdims,gdims,3*sizeof(int));
  memcpy(&Cgr->ldims,gr.ldims,3*sizeof(int));
  memcpy(&Cgr->mem_order,mem_order,3*sizeof(int));
  memcpy(&Cgr->pgrid,pgrid,3*sizeof(int));
  memcpy(&Cgr->proc_order,proc_order,3*sizeof(int));
  memcpy(&Cgr->P,gr.P,3*sizeof(int));
  memcpy(&Cgr->D,gr.D,3*sizeof(int));
  memcpy(&Cgr->L,gr.L,3*sizeof(int));
  memcpy(&Cgr->grid_id,gr.grid_id,3*sizeof(int));
  memcpy(&Cgr->glob_start,gr.glob_start,3*sizeof(int));
  memcpy(&Cgr->mpicomm,gr.mpicomm,3*sizeof(int));
  Cgr->mpi_comm_glob = mpicomm;

  return Cgr;
}

  void p3dfft_free_grid_f(Grid_fort *gr)
  {
    delete gr;
  }

  void p3dfft_free_grid(Grid *gr)
  {
    delete gr;
  }

  void p3dfft_inv_mo(int mo[3],int imo[3]) {
    p3dfft::inv_mo(mo,imo);
  }

  void p3dfft_write_buf(double *ar,char *label,int dims[3],int imo[3]) {
    p3dfft::write_buf<double>(ar,label,dims,imo);
  }

  void p3dfft_exec_3Dtrans_double_f(Plan3D *plan,double *in,double *out,int *OW) {
    return(p3dfft_exec_3Dtrans_double(*plan,in,out,*OW));
  }

  void p3dfft_exec_3Dtrans_double(Plan3D plan,double *in,double *out,int OW) {

    gen_transform3D *trans3D = stored_trans3D[plan];
    
    
    if(trans3D->dt1 == 1)
      if(trans3D->dt2 == 1) {
	transform3D<double,double> *tr3D = (transform3D<double,double> *) trans3D;
	tr3D->exec(in,out,OW);
      }
      else {
	transform3D<double,complex_double> *tr3D = (transform3D<double,complex_double> *) trans3D;
	tr3D->exec(in,(complex_double *) out,OW);
      }
    else
      if(trans3D->dt2 == 1) {
	transform3D<complex_double,double> *tr3D = (transform3D<complex_double,double> *) trans3D;
	tr3D->exec((complex_double *)in,out,OW);
      }
      else {
	transform3D<complex_double,complex_double> *tr3D = (transform3D<complex_double,complex_double> *) trans3D;
	tr3D->exec((complex_double *)in,(complex_double *)out,OW);
      }
    
    if(trans3D->prec != 8) {
      printf("ERror in p3dfft_exec_3Dtrans_double: expecting double precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }

    
}

void p3dfft_exec_3Dtrans_single_f(Plan3D *plan,float *in,float *out,int *OW) 
 {
    return(p3dfft_exec_3Dtrans_single(*plan,in,out,*OW));
  }
void p3dfft_exec_3Dtrans_single(Plan3D plan,float *in,float *out,int OW) {

    gen_transform3D *trans3D = stored_trans3D[plan];
    

    if(trans3D->dt1 == 1)
      if(trans3D->dt2 == 1) {
	transform3D<float,float> *tr3D = (transform3D<float,float> *) trans3D;
	tr3D->exec(in,out,OW);
      }
      else {
	transform3D<float,mycomplex> *tr3D = (transform3D<float,mycomplex> *) trans3D;
	tr3D->exec(in,(mycomplex *) out,OW);
      }
    else
      if(trans3D->dt2 == 1) {
	transform3D<mycomplex,float> *tr3D = (transform3D<mycomplex,float> *) trans3D;
	tr3D->exec((mycomplex *) in,out,OW);
      }
      else {
	transform3D<mycomplex,mycomplex> *tr3D = (transform3D<mycomplex,mycomplex> *) trans3D;
	tr3D->exec((mycomplex *) in,(mycomplex *) out,OW);
      }
    
    if(trans3D->prec != 4) {
      printf("ERror in p3dfft_exec_3Dtrans_single: expecting single precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
    
    
}
  
}
