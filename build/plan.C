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
					    
#if defined(HAVE_CONFIG_H)
#include "config.h"
#endif

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

  




stage *init_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2,int dt, int prec)
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

stage *init_trans_MPIplan(const grid &gr1,const grid &gr2,MPI_Comm mpicomm,int d1,int d2,const gen_trans_type *type,int trans_dim_, bool inplace,int prec)
{
  if(type->dt1 == REAL)
    if(type->dt2 == REAL) {
      if(prec == 4) {
	trans_MPIplan<float,float> *p = new trans_MPIplan<float,float>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	(p->trplan)->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->lib_plan == NULL) {
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
	trans_MPIplan<double,double> *p=new trans_MPIplan<double,double>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->lib_plan == NULL) {
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
	  mo[i] = gr1.mem_order[i];
	}
	gdims[trans_dim_] = gdims[trans_dim_]/2+1;
	double *A=new double[1024];
	grid *tmpgrid = new grid(gdims,pgrid,proc_order,mo,mpicomm_glob);

      if(prec == 4) {
	trans_MPIplan<float,mycomplex> *p= new trans_MPIplan<float,mycomplex>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->lib_plan == NULL) {
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
	trans_MPIplan<double,complex_double> *p= new trans_MPIplan<double,complex_double>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->lib_plan == NULL) {
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
	  mo[i] = gr1.mem_order[i];
	}
	gdims[trans_dim_] = (gdims[trans_dim_]-1)*2;
	grid *tmpgrid = new grid(gdims,pgrid,proc_order,mo,mpicomm_glob);

      if(prec == 4) {
	trans_MPIplan<mycomplex,float> *p = new trans_MPIplan<mycomplex,float>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->lib_plan == NULL) {
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
	trans_MPIplan<complex_double,double> *p= new trans_MPIplan<complex_double,double>(gr1,*tmpgrid,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	delete tmpgrid;
	if(p->trplan->lib_plan == NULL) {
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
	trans_MPIplan<mycomplex,mycomplex> *p= new trans_MPIplan<mycomplex,mycomplex>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->lib_plan == NULL) {
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
	trans_MPIplan<complex_double,complex_double> *p= new trans_MPIplan<complex_double,complex_double>(gr1,gr1,gr2,mpicomm,d1,d2,type,trans_dim_,inplace);
	//	p->trplan->lib_plan = p->trplan->find_plan(p->trplan->trans_type); 
	if(p->trplan->lib_plan == NULL) {
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

stage *init_transplan(const grid &gr1,const grid &gr2,const gen_trans_type *type,int d, bool inplace,int prec)
{
  if(type->dt1 == REAL)
    if(type->dt2 == REAL) {
      if(prec == 4) {
	transplan<float,float> *p = new transplan<float,float>(gr1,gr2,type,d,inplace);
	
	if(!p->trans_type->is_empty) {
	  //	  p->lib_plan = p->find_plan(p->trans_type); 
	  if(p->lib_plan == NULL) {
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
	transplan<double,double> *p =new transplan<double,double>(gr1,gr2,type,d,inplace);
	if(!p->trans_type->is_empty) {
	  //	  p->lib_plan = p->find_plan(p->trans_type); 
	  if(p->lib_plan == NULL) {
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
	transplan<float,mycomplex> *p=new transplan<float,mycomplex>(gr1,gr2,type,d,inplace);
	//  p->lib_plan = p->find_plan(p->trans_type); 
  if(p->lib_plan == NULL) {
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
	transplan<double,complex_double> *p=new transplan<double,complex_double>(gr1,gr2,type,d,inplace);
	//  p->lib_plan = p->find_plan(p->trans_type); 
  if(p->lib_plan == NULL) {
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
	transplan<mycomplex,float> *p=new transplan<mycomplex,float>(gr1,gr2,type,d,inplace);
	//  p->lib_plan = p->find_plan(p->trans_type); 
  if(p->lib_plan == NULL) {
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
	transplan<complex_double,double> *p=new transplan<complex_double,double>(gr1,gr2,type,d,inplace);
	//  p->lib_plan = p->find_plan(p->trans_type); 
  if(p->lib_plan == NULL) {
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
	transplan<mycomplex,mycomplex> *p=new transplan<mycomplex,mycomplex>(gr1,gr2,type,d,inplace);
	if(!p->trans_type->is_empty) {
	  //	  p->lib_plan = p->find_plan(p->trans_type); 
	  if(p->lib_plan == NULL) {
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
	transplan<complex_double,complex_double> *p=new transplan<complex_double,complex_double>(gr1,gr2,type,d,inplace);
	if(!p->trans_type->is_empty) {
	  //	  p->lib_plan = p->find_plan(p->trans_type); 
	  if(p->lib_plan == NULL) {
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
