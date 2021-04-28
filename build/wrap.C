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

using namespace std;
using namespace p3dfft;



extern "C" {

void p3dfft_setup() {
    p3dfft::setup();

P3DFFT_EMPTY_TYPE_SINGLE=p3dfft::EMPTY_TYPE_SINGLE;
P3DFFT_EMPTY_TYPE_DOUBLE=p3dfft::EMPTY_TYPE_DOUBLE;
P3DFFT_EMPTY_TYPE_SINGLE_COMPLEX=p3dfft::EMPTY_TYPE_SINGLE_COMPLEX;
P3DFFT_EMPTY_TYPE_DOUBLE_COMPLEX=p3dfft::EMPTY_TYPE_DOUBLE_COMPLEX;
P3DFFT_R2CFFT_S=p3dfft::R2CFFT_S;
P3DFFT_R2CFFT_D=p3dfft::R2CFFT_D;
P3DFFT_C2RFFT_S=p3dfft::C2RFFT_S;
P3DFFT_C2RFFT_D=p3dfft::C2RFFT_D;
P3DFFT_CFFT_FORWARD_S=p3dfft::CFFT_FORWARD_S;
P3DFFT_CFFT_FORWARD_D=p3dfft::CFFT_FORWARD_D;
P3DFFT_CFFT_BACKWARD_S=p3dfft::CFFT_BACKWARD_S;
P3DFFT_CFFT_BACKWARD_D=p3dfft::CFFT_BACKWARD_D;
P3DFFT_DCT1_REAL_S=p3dfft::DCT1_REAL_S;
P3DFFT_DCT1_REAL_D=p3dfft::DCT1_REAL_D; 
P3DFFT_DST1_REAL_S=p3dfft::DST1_REAL_S;
P3DFFT_DST1_REAL_D=p3dfft::DST1_REAL_D;
P3DFFT_DCT2_REAL_S=p3dfft::DCT2_REAL_S;
P3DFFT_DCT2_REAL_D=p3dfft::DCT2_REAL_D; 
P3DFFT_DST2_REAL_S=p3dfft::DST2_REAL_S;
P3DFFT_DST2_REAL_D=p3dfft::DST2_REAL_D;
P3DFFT_DCT3_REAL_S=p3dfft::DCT3_REAL_S;
P3DFFT_DCT3_REAL_D=p3dfft::DCT3_REAL_D; 
P3DFFT_DST3_REAL_S=p3dfft::DST3_REAL_S;
P3DFFT_DST3_REAL_D=p3dfft::DST3_REAL_D;
P3DFFT_DCT4_REAL_S=p3dfft::DCT4_REAL_S;
P3DFFT_DCT4_REAL_D=p3dfft::DCT4_REAL_D; 
P3DFFT_DST4_REAL_S=p3dfft::DST4_REAL_S;
P3DFFT_DST4_REAL_D=p3dfft::DST4_REAL_D;

P3DFFT_DCT1_COMPLEX_S=p3dfft::DCT1_COMPLEX_S;
P3DFFT_DCT1_COMPLEX_D=p3dfft::DCT1_COMPLEX_D; 
P3DFFT_DST1_COMPLEX_S=p3dfft::DST1_COMPLEX_S;
P3DFFT_DST1_COMPLEX_D=p3dfft::DST1_COMPLEX_D;
P3DFFT_DCT2_COMPLEX_S=p3dfft::DCT2_COMPLEX_S;
P3DFFT_DCT2_COMPLEX_D=p3dfft::DCT2_COMPLEX_D; 
P3DFFT_DST2_COMPLEX_S=p3dfft::DST2_COMPLEX_S;
P3DFFT_DST2_COMPLEX_D=p3dfft::DST2_COMPLEX_D;
P3DFFT_DCT3_COMPLEX_S=p3dfft::DCT3_COMPLEX_S;
P3DFFT_DCT3_COMPLEX_D=p3dfft::DCT3_COMPLEX_D; 
P3DFFT_DST3_COMPLEX_S=p3dfft::DST3_COMPLEX_S;
P3DFFT_DST3_COMPLEX_D=p3dfft::DST3_COMPLEX_D;
P3DFFT_DCT4_COMPLEX_S=p3dfft::DCT4_COMPLEX_S;
P3DFFT_DCT4_COMPLEX_D=p3dfft::DCT4_COMPLEX_D; 
P3DFFT_DST4_COMPLEX_S=p3dfft::DST4_COMPLEX_S;
P3DFFT_DST4_COMPLEX_D=p3dfft::DST4_COMPLEX_D;
}

void p3dfft_cleanup() {
    p3dfft::cleanup();
}

  /////////////// C wrap functions /////////////////////////////


  Type3D p3dfft_init_3Dtype(int types[3]) //,char *name)
{
  trans_type3D tp = trans_type3D(types);
  int count = types3D.size();
  types3D.push_back(tp);
  return(count);
}



  Plan3D p3dfft_plan_3Dtrans(Grid *Cgr1,Grid *Cgr2,Type3D tp){


#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: type3D=%d.  initiating gr1\n",tp);
#endif
  
  ProcGrid *pgrid = stored_proc_grids[Cgr1->pgrid];
  DataGrid *gr1 = new DataGrid(Cgr1->Gdims,Cgr1->dim_conj_sym,pgrid,Cgr1->Dmap,Cgr1->MemOrder);
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: initiating gr2\n");
#endif

  pgrid = stored_proc_grids[Cgr2->pgrid];
  DataGrid *gr2 = new DataGrid(Cgr2->Gdims,Cgr2->dim_conj_sym,pgrid,Cgr2->Dmap,Cgr2->MemOrder);
  trans_type3D *type3D = &types3D[tp];
  gen_transform3D *tr3D;
    
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: new transform3D\n");
#endif

  int L[3];
  bool reverse_steps;
  bool init_steps = find_order(L,type3D, gr1, gr2, &reverse_steps);
  int dt1 = types1D[type3D->types[L[0]]]->dt1;
  int dt2 = types1D[type3D->types[L[2]]]->dt2;

  if(type3D->prec == 4) 
    if(dt1 == 1)
      if(dt2 == 1)
	tr3D = new transform3D<float,float>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<float,mycomplex>(*gr1,*gr2,type3D);
    else
      if(dt2 == 1)
	tr3D = new transform3D<mycomplex,float>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<mycomplex,mycomplex>(*gr1,*gr2,type3D);
    else
    if(dt1 == 1)
      if(dt2 == 1)
	tr3D = new transform3D<double,double>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<double,complex_double>(*gr1,*gr2,type3D);
    else
      if(dt2 == 1)
	tr3D = new transform3D<complex_double,double>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<complex_double,complex_double>(*gr1,*gr2,type3D);

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

int p3dfft_plan_1Dtrans(Grid *Cgr1,Grid *Cgr2,int type_ID,int d)
{
    stage *tr;

    ProcGrid *pgrid = stored_proc_grids[Cgr1->pgrid];
  DataGrid *gr1 = new DataGrid(Cgr1->Gdims,Cgr1->dim_conj_sym,pgrid,Cgr1->Dmap,Cgr1->MemOrder);
pgrid = stored_proc_grids[Cgr2->pgrid];
  DataGrid *gr2 = new DataGrid(Cgr2->Gdims,Cgr2->dim_conj_sym,pgrid,Cgr2->Dmap,Cgr2->MemOrder);

  //  printf("TypeID=%d, type is %s\n",type_ID,types1D[type_ID]->name);

  gen_trans_type *tp = types1D[type_ID];

  if(tp->prec == 4) 
    if(tp->dt1 == 1)
      if(tp->dt2 == 1)
	tr = (stage *) new transplan<float,float>(*gr1,*gr2,tp,d);
      else
	tr = (stage *) new transplan<float,mycomplex>(*gr1,*gr2,tp,d);
    else
      if(tp->dt2 == 1)
	tr = (stage *) new transplan<mycomplex,float>(*gr1,*gr2,tp,d);
      else
	tr = (stage *) new transplan<mycomplex,mycomplex>(*gr1,*gr2,tp,d);
    else
    if(tp->dt1 == 1)
      if(tp->dt2 == 1)
	tr = (stage *) new transplan<double,double>(*gr1,*gr2,tp,d);
      else
	tr = (stage *) new transplan<double,complex_double>(*gr1,*gr2,tp,d);
    else
      if(tp->dt2 == 1)
	tr = (stage *) new transplan<complex_double,double>(*gr1,*gr2,tp,d);
      else
	tr = (stage *) new transplan<complex_double,complex_double>(*gr1,*gr2,tp,d);

  int count = stored_trans1D.size();

  stored_trans1D.push_back(tr);
  delete gr1,gr2;
  return count;

  }

  int find_grid(int gdims[3],int pgrid,int dmap[3],int mem_order[3]) {
    
    int cnt=stored_data_grids.size();
    int i;
    int res;

    for(i=0;i<cnt;i++) {
      DataGrid *gr=stored_data_grids[i];
      if((arcmp(gr->Gdims,gdims,3) == 0) && *(gr->Pgrid) == *stored_proc_grids[pgrid])
	if((arcmp(gr->Dmap,dmap,3) == 0) && (arcmp(gr->MemOrder,mem_order,3)== 0)) {
	  return(i);
      }
    }
    return(-1);
      
}

  int p3dfft_init_proc_grid(int pdims[3],MPI_Comm mpicomm)
  {
    ProcGrid *pgrid = new ProcGrid(pdims,mpicomm);
    stored_proc_grids.push_back(pgrid);
    return(stored_proc_grids.size()-1);
  }


  Grid *p3dfft_init_data_grid(int gdims[3],int dim_conj_sym,int pgrid_id,int dmap[3],int mem_order[3]) {

    ProcGrid *pgrid = stored_proc_grids[pgrid_id];
    DataGrid gr = DataGrid(gdims,dim_conj_sym,pgrid,dmap,mem_order);
    Grid *Cgr = new Grid;
    Cgr->numtasks = gr.Pgrid->numtasks;
    Cgr->taskid = gr.Pgrid->taskid;
    Cgr->nd = gr.nd;
    Cgr->pgrid = pgrid_id;
    Cgr->dim_conj_sym = dim_conj_sym;
    memcpy(&Cgr->Gdims,gdims,3*sizeof(int));
    memcpy(&Cgr->Ldims,gr.Ldims,3*sizeof(int));
    memcpy(&Cgr->MemOrder,mem_order,3*sizeof(int));
    memcpy(&Cgr->Dmap,dmap,3*sizeof(int));
    memcpy(&Cgr->ProcDims,pgrid->ProcDims,3*sizeof(int));
    memcpy(&Cgr->grid_id,gr.grid_id,3*sizeof(int));
    memcpy(&Cgr->GlobStart,gr.GlobStart,3*sizeof(int));
    Cgr->mpi_comm_glob = gr.Pgrid->mpi_comm_glob;

    return Cgr;
}

/*
  void p3dfft_free_grid_f(int *gr)
  {
    stored_grids[*gr].~grid();
//    stored_grids.erase(stored_grids.begin()+ *gr);
  }
*/

  void p3dfft_free_data_grid(Grid *gr)
  {
    delete gr;
  }

  void p3dfft_free_proc_grid(int pgrid_id)
  {
    ProcGrid *pgrid = stored_proc_grids[pgrid_id];
    delete pgrid;
    stored_proc_grids.erase(stored_proc_grids.begin() + pgrid_id);
  }

  void p3dfft_inv_mo(int mo[3],int imo[3]) {
    p3dfft::inv_mo(mo,imo);
  }

  //  void p3dfft_write_buf(double *ar,char *label,int dims[3],int imo[3]) {
  //  p3dfft::write_buf<double>(ar,label,dims,imo);
  // }

  void p3dfft_exec_3Dtrans_single(Plan3D plan,float *in,float *out, int OW) {

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
		//trans3D->grid1->mpi_comm_glob,0);
    }
    
    
}

  void p3dfft_exec_3Dtrans_double(Plan3D plan,double *in,double *out, int OW) {

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

  void p3dfft_exec_3Dderiv_single(Plan3D plan,float *in,float *out,int idir, int OW) {

    gen_transform3D *trans3D = stored_trans3D[plan];
    

    if(trans3D->dt1 == 1)
      if(trans3D->dt2 == 1) {
	transform3D<float,float> *tr3D = (transform3D<float,float> *) trans3D;
	tr3D->exec_deriv(in,out,idir,OW);
      }
      else {
	transform3D<float,mycomplex> *tr3D = (transform3D<float,mycomplex> *) trans3D;
	tr3D->exec_deriv(in,(mycomplex *) out,idir,OW);
      }
    else
      if(trans3D->dt2 == 1) {
	transform3D<mycomplex,float> *tr3D = (transform3D<mycomplex,float> *) trans3D;
	tr3D->exec_deriv((mycomplex *) in,out,idir,OW);
      }
      else {
	transform3D<mycomplex,mycomplex> *tr3D = (transform3D<mycomplex,mycomplex> *) trans3D;
	tr3D->exec_deriv((mycomplex *) in,(mycomplex *) out,idir,OW);
      }
    
    if(trans3D->prec != 4) {
      printf("ERror in p3dfft_exec_3Dderiv_single: expecting single precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
    
    
}

  void p3dfft_exec_3Dderiv_double(Plan3D plan,double *in,double *out,int idir, int OW) {

    gen_transform3D *trans3D = stored_trans3D[plan];
    
    
    if(trans3D->dt1 == 1)
      if(trans3D->dt2 == 1) {
	transform3D<double,double> *tr3D = (transform3D<double,double> *) trans3D;
	tr3D->exec_deriv(in,out,idir,OW);
      }
      else {
	transform3D<double,complex_double> *tr3D = (transform3D<double,complex_double> *) trans3D;
	tr3D->exec_deriv(in,(complex_double *) out,idir,OW);
      }
    else
      if(trans3D->dt2 == 1) {
	transform3D<complex_double,double> *tr3D = (transform3D<complex_double,double> *) trans3D;
	tr3D->exec_deriv((complex_double *)in,out,idir,OW);
      }
      else {
	transform3D<complex_double,complex_double> *tr3D = (transform3D<complex_double,complex_double> *) trans3D;
	tr3D->exec_deriv((complex_double *)in,(complex_double *)out,idir,OW);
      }
    
    if(trans3D->prec != 8) {
      printf("ERror in p3dfft_exec_3Dderiv_double: expecting double precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
		//trans3D->grid1->mpi_comm_glob,0);
    }

    
}

  void p3dfft_exec_1Dtrans_double(int plan,double *in,double *out,int OW) {

    stage *trans = stored_trans1D[plan];
    
    
    if(trans->dt1 == 1)
      if(trans->dt2 == 1) {
	transplan<double,double> *tr = (transplan<double,double> *) trans;
	tr->exec((char *) in,(char *) out,OW);
      }
      else {
	transplan<double,complex_double> *tr = (transplan<double,complex_double> *) trans;
	tr->exec((char *) in,(char *) out,OW);
      }
    else
      if(trans->dt2 == 1) {
	transplan<complex_double,double> *tr = (transplan<complex_double,double> *) trans;
	tr->exec((char *) in,(char *) out,OW);
      }
      else {
	transplan<complex_double,complex_double> *tr = (transplan<complex_double,complex_double> *) trans;
	tr->exec((char *) in,(char *) out,OW);
      }

    /*    
    if(tr->prec != 8) {
      printf("ERror in p3dfft_exec_1Dtrans_double: expecting double precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
    */
    
  }

  void p3dfft_exec_1Dtrans_single(int plan,float *in,float *out, int OW) {

    stage *trans = stored_trans1D[plan];
    int prec;
    
    if(trans->dt1 == 1)
      if(trans->dt2 == 1) {
	transplan<float,float> *tr = (transplan<float,float> *) trans;
	tr->exec((char *) in,(char *) out, OW);
      }
      else {
	transplan<float,mycomplex> *tr = (transplan<float,mycomplex> *) trans;
	tr->exec((char *) in,(char *)  out, OW);
      }
    else
      if(trans->dt2 == 1) {
	transplan<mycomplex,float> *tr = (transplan<mycomplex,float> *) trans;
	tr->exec((char *) in,(char *) out, OW);
      }
      else {
	transplan<mycomplex,mycomplex> *tr = (transplan<mycomplex,mycomplex> *) trans;
	tr->exec((char *) in,(char *) out, OW);
      }
    
    /*
    if(tr->prec != 4) {
      printf("ERror in p3dfft_exec_1Dtrans_single: expecting single precision data\n");
      MPI_Abort(MPI_COMM_WORLD,0);
    }
    */
    
  }

  
  void p3dfft_compute_deriv_single(float *in,float *out,Grid *Cgrid,int idir) {

    ProcGrid *pgrid = stored_proc_grids[Cgrid->pgrid];
    DataGrid *grid1 = new DataGrid(Cgrid->Gdims,Cgrid->dim_conj_sym,pgrid,Cgrid->Dmap,Cgrid->MemOrder);

     compute_deriv<mycomplex>((mycomplex *) in,(mycomplex *) out,grid1,idir);
    delete grid1;

  }

  void p3dfft_compute_deriv_double(double *in,double *out,Grid *Cgrid,int idir) {

    ProcGrid *pgrid = stored_proc_grids[Cgrid->pgrid];
    DataGrid *grid1 = new DataGrid(Cgrid->Gdims,Cgrid->dim_conj_sym,pgrid,Cgrid->Dmap,Cgrid->MemOrder);

    compute_deriv<complex_double>((complex_double *) in,(complex_double *) out,grid1,idir);
    delete grid1;

  }
 


  ///////////////// Fortran wrap functions ////////////////////////

  void p3dfft_init_3Dtype_f(int *type,int types[3]) //,char *name)
{
  trans_type3D tp = trans_type3D(types);
  int count = types3D.size();
  types3D.push_back(tp);
  *type = count;
  //  return(count);
}

  void p3dfft_plan_1Dtrans_f(int *plan,int *Fgr1,int *Fgr2,int *type_ID,int *d)
{
  DataGrid *gr1 = stored_data_grids[*Fgr1];
  DataGrid *gr2 = stored_data_grids[*Fgr2];
  stage *tr;
  gen_trans_type *tp = types1D[*type_ID];

  if(tp->prec == 4) 
    if(tp->dt1 == 1)
      if(tp->dt2 == 1)
	tr = new transplan<float,float>(*gr1,*gr2,tp,*d);
      else
	tr = new transplan<float,mycomplex>(*gr1,*gr2,tp,*d);
    else
      if(tp->dt2 == 1)
	tr = new transplan<mycomplex,float>(*gr1,*gr2,tp,*d);
      else
	tr = new transplan<mycomplex,mycomplex>(*gr1,*gr2,tp,*d);
    else
    if(tp->dt1 == 1)
      if(tp->dt2 == 1)
	tr = new transplan<double,double>(*gr1,*gr2,tp,*d);
      else
	tr = new transplan<double,complex_double>(*gr1,*gr2,tp,*d);
    else
      if(tp->dt2 == 1)
	tr = new transplan<complex_double,double>(*gr1,*gr2,tp,*d);
      else
	tr = new transplan<complex_double,complex_double>(*gr1,*gr2,tp,*d);


  int count = stored_trans1D.size();

  stored_trans1D.push_back(tr);
  //  delete gr1,gr2;
  *plan = count;
  //  return count;

  }

  void p3dfft_compute_deriv_single_f(float *in,float *out,int *igrid,int *idir) {
    DataGrid *grid1 = stored_data_grids[*igrid];

    compute_deriv<mycomplex>((mycomplex *) in,(mycomplex *) out,grid1,*idir-1);
  }

  void p3dfft_compute_deriv_double_f(double *in,double *out,int *igrid,int *idir) {
    DataGrid *grid1 = stored_data_grids[*igrid];

    compute_deriv<complex_double>((complex_double *) in,(complex_double *) out,grid1,*idir-1);
  }

  void p3dfft_plan_3Dtrans_f(Plan3D *plan,int *Fgr1,int *Fgr2,Type3D *tp){

  /* 
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: type3D=%d.  initiating gr1\n",*tp);
#endif
 grid *gr1 = new grid(Fgr1->gdims,Fgr1->pgrid,Fgr1->proc_order,Fgr1->mem_order,MPI_Comm_f2c(Fgr1->mpi_comm_glob));
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: initiating gr2\n");
#endif
  */


    DataGrid *gr1 = stored_data_grids[*Fgr1];
    DataGrid *gr2 = stored_data_grids[*Fgr2];

  //  grid *gr2 = new grid(Fgr2->gdims,Fgr2->pgrid,Fgr2->proc_order,Fgr2->mem_order,MPI_Comm_f2c(Fgr2->mpi_comm_glob));
  trans_type3D *type3D = &types3D[*tp];
  gen_transform3D *tr3D;
    
#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: new transform3D\n");
#endif

  int L[3];
  bool reverse_steps;
  bool init_steps = find_order(L,type3D, gr1, gr2, &reverse_steps);
  
  int dt1 = types1D[type3D->types[L[0]]]->dt1;
  int dt2 = types1D[type3D->types[L[2]]]->dt2;

  if(type3D->prec == 4) 
    if(dt1 == 1)
      if(dt2 == 1)
	tr3D = new transform3D<float,float>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<float,mycomplex>(*gr1,*gr2,type3D);
    else
      if(dt2 == 1)
	tr3D = new transform3D<mycomplex,float>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<mycomplex,mycomplex>(*gr1,*gr2,type3D);
    else
    if(dt1 == 1)
      if(dt2 == 1)
	tr3D = new transform3D<double,double>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<double,complex_double>(*gr1,*gr2,type3D);
    else
      if(dt2 == 1)
	tr3D = new transform3D<complex_double,double>(*gr1,*gr2,type3D);
      else
	tr3D = new transform3D<complex_double,complex_double>(*gr1,*gr2,type3D);

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: count\n");
#endif

  int count = stored_trans3D.size();

#ifdef DEBUG
  printf("p3dfft_plan_3Dtrans: push back\n");
#endif

  stored_trans3D.push_back(tr3D);
  //  delete gr1,gr2;
  *plan = count;
  //return count;
    
  /*
  gen_trans_type *tp1D[3];
  int i;

  for(i=0;i<3;i++) 
    tp1D[i] = types1D[type3D->types[i]];

  int dt1 = tp1D[0]->dt1;
  int dt2 = tp1D[2]->dt2;
  */

}

  int p3dfft_init_proc_grid_f(int *pdims,int *mpicomm)
  {
    ProcGrid *pgrid = new ProcGrid(pdims,MPI_Comm_f2c(*mpicomm));
    stored_proc_grids.push_back(pgrid);
    return(stored_proc_grids.size()-1);
  }

  void p3dfft_init_data_grid_f(int *mygrid,int *ldims,int *glob_start,int *gdims,int *dim_conj_sym,int *pgrid_id,int *dmap,int *mem_order) {
    
    int num=find_grid(gdims,*pgrid_id,dmap,mem_order); //,MPI_Comm_f2c(*mpicomm));
    if(num >= 0) 
      *mygrid = num;
    else {

    DataGrid *gr1;
    ProcGrid *pgrid=stored_proc_grids[*pgrid_id];
    gr1 = new DataGrid(gdims,*dim_conj_sym,pgrid,dmap,mem_order);
    memcpy(ldims,gr1->Ldims,3*sizeof(int));
    memcpy(glob_start,gr1->GlobStart,3*sizeof(int));
    num = stored_data_grids.size();
    stored_data_grids.push_back(gr1);
    *mygrid = num;
//  return(num);
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


  void p3dfft_exec_3Dtrans_double_f(Plan3D *plan,double *in,double *out, int *OW) {
    return(p3dfft_exec_3Dtrans_double(*plan,in,out,*OW));
  }


  void p3dfft_exec_3Dtrans_single_f(Plan3D *plan,float *in,float *out, int *OW) 
 {
   return(p3dfft_exec_3Dtrans_single(*plan,in,out, *OW));
  }
  void p3dfft_exec_3Dderiv_double_f(Plan3D *plan,double *in,double *out,int *idir, int *OW) {
    return(p3dfft_exec_3Dderiv_double(*plan,in,out,*idir-1, *OW));
  }


  void p3dfft_exec_3Dderiv_single_f(Plan3D *plan,float *in,float *out,int *idir, int *OW) 
 {
   return(p3dfft_exec_3Dderiv_single(*plan,in,out,*idir-1,*OW));
  }
 
  void p3dfft_exec_1Dtrans_double_f(int *plan,double *in,double *out, int *OW) {
    return(p3dfft_exec_1Dtrans_double(*plan,in,out,*OW));
  }


  void p3dfft_exec_1Dtrans_single_f(int *plan,float *in,float *out,int *OW) {
    return(p3dfft_exec_1Dtrans_single(*plan,in,out, *OW));
  }


 
}


