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

      subroutine p3dfft_1Dtrans_double(plan,IN,OUT)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_1Dtrans_double(plan,IN,OUT)

      return
      end subroutine

      subroutine p3dfft_1Dtrans_single(plan,IN,OUT)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_1Dtrans_single(plan,IN,OUT)

      return
      end subroutine

      subroutine p3dfft_3Dtrans_double(plan,IN,OUT,OW)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_3Dtrans_double(plan,IN,OUT,OW)

      return
      end subroutine

      subroutine p3dfft_3Dtrans_single(plan,IN,OUT,OW)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer OW,plan

      call p3dfft_exec_3Dtrans_single(plan,IN,OUT,OW)

      return
      end subroutine
