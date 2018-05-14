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
integer, bind(C,name='P3DFFT_COSTRAN_REAL_S') :: P3DFFT_COSTRAN_REAL_S
integer, bind(C,name='P3DFFT_COSTRAN_REAL_D') :: P3DFFT_COSTRAN_REAL_D
integer, bind(C,name='P3DFFT_SINTRAN_REAL_S') :: P3DFFT_SINTRAN_REAL_S
integer, bind(C,name='P3DFFT_SINTRAN_REAL_D') :: P3DFFT_SINTRAN_REAL_D


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
      integer(C_INT) function p3dfft_init_3Dtype(types) bind(C, name='p3dfft_init_3Dtype')
      import
      integer(C_INT) :: types(3)

      end function

      subroutine p3dfft_setup() bind(C,name="p3dfft_setup")
      end subroutine

      subroutine p3dfft_cleanup() bind(C,name="p3dfft_cleanup")
      end subroutine

      integer(C_INT) function p3dfft_init_grid(ldims,glob_start,gdims,pgrid,proc_order,mem_order,mpicomm) bind(C,name='p3dfft_init_grid_f')
        use iso_c_binding
      import 
!      type(grid) :: gr
      integer(C_INT), dimension(3) :: gdims,pgrid,proc_order,mem_order,ldims,glob_start
      integer(C_INT) :: mpicomm
    end function p3dfft_init_grid

!      subroutine p3dfft_free_grid(gr) bind(C,name='p3dfft_free_grid_f')
!      import
!      integer(C_INT) :: gr
!      end subroutine

      subroutine p3dffft_inv_mo(in,out) bind(C,name='p3dfft_inv_mo')
      import
      integer(C_INT) :: in(3),out(3)
      end subroutine

      subroutine p3dfft_exec_1Dtrans_double_f(plan,in,out) bind(C,name='p3dfft_exec_1Dtrans_double_f')
      import
      integer(C_INT) plan
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_1Dtrans_single_f(plan,in,out) bind(C,name='p3dfft_exec_1Dtrans_single_f')
      import
      integer(C_INT) plan
      real(C_FLOAT), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dtrans_double_f(plan,in,out,OW) bind(C,name='p3dfft_exec_3Dtrans_double_f')
      import
      integer(C_INT) plan,OW;
      real(C_DOUBLE), dimension(*) :: in,out
      end subroutine

      subroutine p3dfft_exec_3Dtrans_single_f(plan,in,out,OW) bind(C,name='p3dfft_exec_3Dtrans_single_f')
      import
      integer(C_INT) plan,OW;
      real(C_FLOAT), dimension(*) :: in,out
      end subroutine

      integer(C_INT) function p3dfft_plan_1Dtrans_f(grid1,grid2,trans_ID,dim,inplace) bind(C, name='p3dfft_plan_1Dtrans_f')
      import
      integer(C_INT) :: grid1,grid2,dim
      integer(C_INT) :: trans_ID,inplace

      end function

      integer(C_INT) function p3dfft_plan_3Dtrans_f(grid1,grid2,trans_ID,inplace) bind(C, name='p3dfft_plan_3Dtrans_f')
      import
      integer(C_INT) :: grid1,grid2
      integer(C_INT) :: trans_ID,inplace

      end function
end interface

end module
      

