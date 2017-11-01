      subroutine p3dfft_3Dtrans_double_f(plan,IN,OUT,OW)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_3Dtrans_double_f(plan,IN,OUT,OW)

      return
      end subroutine

      subroutine p3dfft_3Dtrans_single_f(plan,IN,OUT,OW)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_3Dtrans_single_f(plan,IN,OUT,OW)

      return
      end subroutine
