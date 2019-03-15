
! Title: P3DFFT++ library

! Authors: Dmitry Pekurovsky

! Copyright (c) 2006-2019 

! The Regents of the University of California.

! All Rights Reserved.                        

 

!    Permission to use, copy, modify and  distribute  any part

!    of this software for  educational,  research  and  non-profit

!    purposes, by individuals or non-profit organizations,

!    without fee,  and  without a written  agreement is

!    hereby granted,  provided  that the  above  copyright notice,

!    this paragraph  and the following  three  paragraphs appear in

!    all copies.       

 

!    For-profit organizations desiring to use this software and others

!    wishing to incorporate this  software into commercial

!    products or use it for  commercial  purposes should contact the:    

!          Office of Innovation & Commercialization 

!          University of California San Diego

!          9500 Gilman Drive,  La Jolla,  California, 92093-0910        

!          Phone: (858) 534-5815

!          E-mail: innovation@ucsd.edu

 

!    IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE

!    TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR    

!    CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT

!    OF THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF

!    CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH

!    DAMAGE.

 

!    THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND

!    THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE        

!    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 

!    THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND    

!    EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR

!    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES

!    OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR

!    THAT THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,        

!    TRADEMARK OR OTHER RIGHTS.

module p3dfft_plus_plus

use iso_c_binding
implicit none
!include 'mpif.h'


integer, bind(C,name='P3DFFT_EMPTY_TYPE') :: P3DFFT_EMPTY_TYPE
integer, bind(C,name='P3DFFT_R2CFFT_S') :: P3DFFT_R2CFFT_S
integer, bind(C,name='P3DFFT_R2CFFT_D') :: P3DFFT_R2CFFT_D
integer, bind(C,name='P3DFFT_C2RFFT_S') :: P3DFFT_C2RFFT_S
integer, bind(C,name='P3DFFT_C2RFFT_D') :: P3DFFT_C2RFFT_D
integer, bind(C,name='P3DFFT_CFFT_FORWARD_S') :: P3DFFT_CFFT_FORWARD_S
integer, bind(C,name='P3DFFT_CFFT_FORWARD_D') :: P3DFFT_CFFT_FORWARD_D
integer, bind(C,name='P3DFFT_CFFT_BACKWARD_S') :: P3DFFT_CFFT_BACKWARD_S
integer, bind(C,name='P3DFFT_CFFT_BACKWARD_D') :: P3DFFT_CFFT_BACKWARD_D
integer, bind(C,name='P3DFFT_DCT1_REAL_S') :: P3DFFT_DCT1_REAL_S
integer, bind(C,name='P3DFFT_DCT1_REAL_D') :: P3DFFT_DCT1_REAL_D
integer, bind(C,name='P3DFFT_DST1_REAL_S') :: P3DFFT_DST1_REAL_S
integer, bind(C,name='P3DFFT_DST1_REAL_D') :: P3DFFT_DST1_REAL_D
integer, bind(C,name='P3DFFT_DCT1_COMPLEX_S') :: P3DFFT_DCT1_COMPLEX_S
integer, bind(C,name='P3DFFT_DCT1_COMPLEX_D') :: P3DFFT_DCT1_COMPLEX_D
integer, bind(C,name='P3DFFT_DST1_COMPLEX_S') :: P3DFFT_DST1_COMPLEX_S
integer, bind(C,name='P3DFFT_DST1_COMPLEX_D') :: P3DFFT_DST1_COMPLEX_D


type, public, bind(C) :: grid
      integer(C_INT) :: taskid
      integer(C_INT) :: numtasks
      integer(C_INT) :: nd
      integer(C_INT) :: gdims(3)
      integer(C_INT) :: mem_order(3)
      integer(C_INT) :: ldims(3)
      integer(C_INT) :: pgrid(3)
      integer(C_INT) :: proc_order(3)
!      integer(C_INT) :: P(3)
!      integer(C_INT) :: D(3)
!      integer(C_INT) :: L(3)
!      integer(C_INT) :: grid_id(3)
!      integer(C_INT) :: grid_id_cart(3)
     integer(C_INT) :: glob_start(3)
      integer(C_INT) :: mpi_comm_glob
!
!      integer(C_INT) :: mpi_comm_start
!      integer(C_INT) :: mpicomm(3)
end type grid 

interface

      ! integer(C_INT) 
      subroutine p3dfft_init_3Dtype(mytype,types) bind(C, name='p3dfft_init_3Dtype_f')
      import 
      integer(C_INT) :: mytype
      integer(C_INT) :: types(3)

      end subroutine

      subroutine p3dfft_setup() bind(C,name="p3dfft_setup")
      end subroutine

      subroutine p3dfft_cleanup() bind(C,name="p3dfft_cleanup")
      end subroutine

!      integer(C_INT) 
      subroutine p3dfft_init_grid(mygrid,ldims,glob_start,gdims,dim_conj_sym,pgrid,proc_order,mem_order,mpicomm) bind(C,name='p3dfft_init_grid_f')
!        use iso_c_binding
      import 
!      type(grid) :: gr
      integer(C_INT), dimension(3) :: gdims,pgrid,proc_order,mem_order,ldims,glob_start
      integer(C_INT) :: mygrid,mpicomm,dim_conj_sym
    end subroutine p3dfft_init_grid

!      subroutine p3dfft_free_grid(gr) bind(C,name='p3dfft_free_grid_f')
!      import
!      integer(C_INT) :: gr
!      end subroutine

      subroutine p3dffft_inv_mo(in,out) bind(C,name='p3dfft_inv_mo')
      import
      integer(C_INT) :: in(3),out(3)
      end subroutine

      subroutine p3dfft_exec_1Dtrans_double(plan,in,out) bind(C,name='p3dfft_exec_1Dtrans_double_f')
      import
      integer(C_INT) plan
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_1Dtrans_single(plan,in,out) bind(C,name='p3dfft_exec_1Dtrans_single_f')
      import
      integer(C_INT) plan
      real(C_FLOAT), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dtrans_double(plan,in,out,OW) bind(C,name='p3dfft_exec_3Dtrans_double_f')
      import
      integer(C_INT) plan,OW;
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dderiv_single(plan,in,out,idir,OW) bind(C,name='p3dfft_exec_3Dderiv_single_f')
      import
      integer(C_INT) plan,OW,idir;
      real(C_FLOAT), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dderiv_double(plan,in,out,idir,OW) bind(C,name='p3dfft_exec_3Dderiv_double_f')
      import
      integer(C_INT) plan,OW,idir;
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dtrans_single(plan,in,out,OW) bind(C,name='p3dfft_exec_3Dtrans_single_f')
      import
      integer(C_INT) plan,OW;
      real(C_FLOAT), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_plan_1Dtrans(myplan,grid1,grid2,trans_ID,dim,inplace) bind(C, name='p3dfft_plan_1Dtrans_f')
      import
      integer(C_INT) :: myplan,grid1,grid2,dim
      integer(C_INT) :: trans_ID,inplace

      end subroutine

      subroutine p3dfft_plan_3Dtrans(myplan,grid1,grid2,trans_ID,inplace) bind(C, name='p3dfft_plan_3Dtrans_f')
      import
      integer(C_INT) :: myplan,grid1,grid2
      integer(C_INT) :: trans_ID,inplace

      end subroutine

      subroutine p3dfft_compute_deriv_single(in,out,grid,idir) bind(C,name='p3dfft_compute_deriv_single_f')
      import
      integer(C_INT) grid,idir
      real(C_FLOAT ), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_compute_deriv_double(in,out,grid,idir) bind(C,name='p3dfft_compute_deriv_double_f')
      import
      integer(C_INT) grid,idir
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine


end interface

end module
      

