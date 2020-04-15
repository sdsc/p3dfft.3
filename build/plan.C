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
					    
#include "templ.C"

namespace p3dfft {

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

  




stage *init_MPIplan(const grid &gr1,const grid &gr2,int mpicomm,int d1,int d2,int dt, int prec)
{
  if(dt == REAL) {
    if(prec == 4) 
      return(new MPIplan<float>(gr1,gr2,mpicomm,d1,d2,prec));
    else if(prec == 8) 
      return(new MPIplan<double>(gr1,gr2,mpicomm,d1,d2,prec));
  }
  else if(dt == COMPLEX) {
    if(prec == 4) 
      return(new MPIplan<mycomplex>(gr1,gr2,mpicomm,d1,d2,prec));
    else if(prec == 8) 
      return(new MPIplan<complex_double>(gr1,gr2,mpicomm,d1,d2,prec));
  }
  else {
    cout << "Unknown datatype" << dt << endl;
    return(NULL);
    }
}

stage *init_trans_MPIplan(const grid &gr1,const grid &gr2,int mpicomm,int d1,int d2,const gen_trans_type *type,int trans_dim_, int prec)
{


  if(type->dt1 == REAL)
    if(type->dt2 == REAL) {
      if(prec == 4) {
	trans_MPIplan<float,float> *p;
	if(gr1.mem_order[trans_dim_] == 0)
	  p = new trans_MPIplan<float,float>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_);
	else if(gr2.mem_order[trans_dim_] == 0) {
	  grid tmpgrid(gr1);
	  for(int i=0;i<3;i++)
	    tmpgrid.mem_order[i] = gr2.mem_order[i];
	  p = new trans_MPIplan<float,float>(gr1,tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	}
  else 
    printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	//	(p->trplan)->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
    //    p->is_trans = true;
    //p->is_mpi = false;
	}
	return(p);
      }
      else if(prec == 8) {
	trans_MPIplan<double,double> *p;
  if(gr1.mem_order[trans_dim_] == 0)
	p=new trans_MPIplan<double,double>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_);
  else if(gr2.mem_order[trans_dim_] == 0) {
	  grid tmpgrid(gr1);
	  for(int i=0;i<3;i++)
	    tmpgrid.mem_order[i] = gr2.mem_order[i];
	  p=new trans_MPIplan<double,double>(gr1,tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
  }
  else 
    printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
  if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  //p->is_trans = true;
	  //p->is_mpi = false;
	}
	return(p);
      }
    }
    else if(type->dt2 == COMPLEX) {
	int gdims[3],pgrid[3],mo[3],proc_order[3],d,i;
	MPI_Comm mpicomm_glob = gr1.mpi_comm_glob;
	
	for(i=0;i<3;i++) {
	  gdims[i] = gr1.gdims[i];
	  pgrid[i] = gr1.pgrid[i];
	  proc_order[i] = gr1.proc_order[i];
	  if(gr1.mem_order[trans_dim_] == 0)
	    mo[i] = gr1.mem_order[i];
	  else if(gr2.mem_order[trans_dim_] == 0)
	    mo[i] = gr2.mem_order[i];
	  else 
	    printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	}
	gdims[trans_dim_] = gdims[trans_dim_]/2+1;
	grid *tmpgrid = new grid(gdims,gr1.dim_conj_sym,pgrid,proc_order,mo,mpicomm_glob);

      if(prec == 4) {
	trans_MPIplan<float,mycomplex> *p= new trans_MPIplan<float,mycomplex>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);

	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  //	  p->is_trans = true;
	  //p->is_mpi = false;
	}
	return(p);
      }
      else if(prec == 8) {
	trans_MPIplan<double,complex_double> *p= new trans_MPIplan<double,complex_double>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  // p->is_trans = true;
	  // p->is_mpi = false;
	}
	return(p);
      }
    }
    else {
      cout << "Unknown datatype" << type->dt2 << endl;
      return NULL;
    }
  else if(type->dt1 == COMPLEX)
    if(type->dt2 == REAL) {
	int gdims[3],pgrid[3],mo[3],proc_order[3],d,i;
	MPI_Comm mpicomm_glob = gr1.mpi_comm_glob;
	
	for(i=0;i<3;i++) {
	  gdims[i] = gr1.gdims[i];
	  pgrid[i] = gr1.pgrid[i];
	  proc_order[i] = gr1.proc_order[i];
	  if(gr1.mem_order[trans_dim_] == 0)
	    mo[i] = gr1.mem_order[i];
	  else if(gr2.mem_order[trans_dim_] == 0)
	    mo[i] = gr2.mem_order[i];
	  else 
	    printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	}
	gdims[trans_dim_] = (gdims[trans_dim_]-1)*2;
	grid *tmpgrid = new grid(gdims,gr1.dim_conj_sym,pgrid,proc_order,mo,mpicomm_glob);

      if(prec == 4) {
	trans_MPIplan<mycomplex,float> *p = new trans_MPIplan<mycomplex,float>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  // p->is_trans = true;
	  // p->is_mpi = false;
	}
	return(p);
      }
      else if(prec == 8) {
	trans_MPIplan<complex_double,double> *p= new trans_MPIplan<complex_double,double>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  //  p->is_trans = true;
	  // p->is_mpi = false;
	}
	return(p);
      }
    }
    else if(type->dt2 == COMPLEX) {
      if(prec == 4) {
	trans_MPIplan<mycomplex,mycomplex> *p;
	  if(gr1.mem_order[trans_dim_] == 0)
	p= new trans_MPIplan<mycomplex,mycomplex>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_);
	  else if(gr2.mem_order[trans_dim_] == 0) {
	  grid tmpgrid(gr1);
	  for(int i=0;i<3;i++)
	    tmpgrid.mem_order[i] = gr2.mem_order[i];
	p= new trans_MPIplan<mycomplex,mycomplex>(gr1,tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	  }
	  else 
	    printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	  if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
	  //  p->is_trans = true;
	  //p->is_mpi = false;
	}
	return(p);
      }
      else if(prec == 8) {
	trans_MPIplan<complex_double,complex_double> *p;
	if(gr1.mem_order[trans_dim_] == 0)
	  p= new trans_MPIplan<complex_double,complex_double>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_);
	else if(gr2.mem_order[trans_dim_] == 0) {
	  grid tmpgrid(gr1);
	  for(int i=0;i<3;i++)
	    tmpgrid.mem_order[i] = gr2.mem_order[i];
	  p= new trans_MPIplan<complex_double,complex_double>(gr1,tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_);
	}
	else 
	  printf("Error in init_trans_MPIplan: neither input %d nor output %d has current dimension as leading\n",gr1.mem_order[trans_dim_],gr2.mem_order[trans_dim_]);
	//	p->trplan->libplan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->plan == NULL) {
	  cout << "Error in trans_plan: null plan" << endl;
	  return(NULL);
	}
	else {
	  p->is_set = true;
    // p->is_trans = true;
    //p->is_mpi = false;
	}
	return(p);
      }
    }
    else {
      cout << "Unknown datatype" << type->dt2 << endl;
      return(NULL);
      
    }
  else {
    cout << "Unknown datatype" << type->dt1 << endl;
    return(NULL);
    
  }
}

stage *init_transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d, int prec)
{
  if(type->dt1 == REAL)
    if(type->dt2 == REAL) {
      if(prec == 4) {
	transplan<float,float> *p = new transplan<float,float>(gr1,gr2,type,d);
	
	if(!p->trans_type->is_empty) {
	  //	  p->libplan = p->find_plan(p->trans_type); 
	  if(p->plan == NULL) {
	    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
	  else {
	    p->is_set = true;
	    //	    p->is_trans = true;
	    // p->is_mpi = false;
	  }
	}
	return(p);
      }
      else if(prec == 8) {
	transplan<double,double> *p =new transplan<double,double>(gr1,gr2,type,d);
	if(!p->trans_type->is_empty) {
	  //	  p->libplan = p->find_plan(p->trans_type); 
	  if(p->plan == NULL) {
	    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
	  else {
	    p->is_set = true;
	    //   p->is_trans = true;
	    // p->is_mpi = false;
	  }
	}
	return(p);
      }
    }
    else if(type->dt2 == COMPLEX) {
      if(prec == 4) {
	transplan<float,mycomplex> *p=new transplan<float,mycomplex>(gr1,gr2,type,d);
	//  p->libplan = p->find_plan(p->trans_type); 
  if(p->plan == NULL) {
    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
  else {
    p->is_set = true;
    //  p->is_trans = true;
    // p->is_mpi = false;
  }
  return(p);
      }
      else if(prec == 8) {
	transplan<double,complex_double> *p=new transplan<double,complex_double>(gr1,gr2,type,d);
	//  p->libplan = p->find_plan(p->trans_type); 
  if(p->plan == NULL) {
    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
  else {
    p->is_set = true;
    // p->is_trans = true;
    //p->is_mpi = false;
  }
  return(p);
      }
    }
    else {
      cout << "Unknown datatype" << type->dt2 << endl;
	    return(NULL);
	 
    }
  else if(type->dt1 == COMPLEX)
    if(type->dt2 == REAL) {
      if(prec == 4) {
	transplan<mycomplex,float> *p=new transplan<mycomplex,float>(gr1,gr2,type,d);
	//  p->libplan = p->find_plan(p->trans_type); 
  if(p->plan == NULL) {
    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
  else {
    p->is_set = true;
    // p->is_trans = true;
    //p->is_mpi = false;
  }
  return(p);
      }
      else if(prec == 8) {
	transplan<complex_double,double> *p=new transplan<complex_double,double>(gr1,gr2,type,d);
	//  p->libplan = p->find_plan(p->trans_type); 
  if(p->plan == NULL) {
    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
  else {
    p->is_set = true;
    // p->is_trans = true;
    //p->is_mpi = false;
  }
  return(p);
      }
    }
    else if(type->dt2 == COMPLEX) {
      if(prec == 4) {
	transplan<mycomplex,mycomplex> *p=new transplan<mycomplex,mycomplex>(gr1,gr2,type,d);
	if(!p->trans_type->is_empty) {
	  //	  p->libplan = p->find_plan(p->trans_type); 
	  if(p->plan == NULL) {
	    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
	  else {
	    p->is_set = true;
	    //	    p->is_trans = true;
	    //p->is_mpi = false;
	  }
	}
	return(p);
      }
      else if(prec == 8) {
	transplan<complex_double,complex_double> *p=new transplan<complex_double,complex_double>(gr1,gr2,type,d);
	if(!p->trans_type->is_empty) {
	  //	  p->libplan = p->find_plan(p->trans_type); 
	  if(p->plan == NULL) {
	    cout << "Error in trans_plan: null plan" << endl;
	    return(NULL);
	  }
	  else {
	    p->is_set = true;
	    //  p->is_trans = true;
	    // p->is_mpi = false;
	  }
	}
	return(p);
      }
    }
    else {
      cout << "Unknown datatype" << type->dt2 << endl;
	    return(NULL);
    }
  else {
    cout << "Unknown datatype" << type->dt1 << endl;
	    return(NULL);
    }

}


}

