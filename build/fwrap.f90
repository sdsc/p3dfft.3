
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


      subroutine p3dfft_1Dtrans_double(plan,IN,OUT,dim_deriv,OW)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer OW,plan,dim_deriv

      call p3dfft_exec_1Dtrans_double(plan,IN,OUT,dim_deriv,OW)

      return
      end subroutine

      subroutine p3dfft_1Dtrans_single(plan,IN,OUT,dim_deriv,OW)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer OW,plan,dim_deriv

      call p3dfft_exec_1Dtrans_single(plan,IN,OUT,dim_deriv,OW)

      return
      end subroutine

      subroutine p3dfft_3Dtrans_double(plan,IN,OUT,OW)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer plan,OW

      call p3dfft_exec_3Dtrans_double(plan,IN,OUT,OW)

      return
      end subroutine

      subroutine p3dfft_3Dtrans_single(plan,IN,OUT,OW)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer plan,OW

      call p3dfft_exec_3Dtrans_single(plan,IN,OUT,OW)

      return
      end subroutine

      subroutine p3dfft_3Dderiv_single(plan,IN,OUT,idir,OW)
      use p3dfft_plus_plus


      real, TARGET :: IN(1,1,*)
      real, TARGET :: OUT(1,1,*)
      integer plan,idir,OW

      call p3dfft_exec_3Dderiv_single(plan,IN,OUT,idir,OW)

      return
      end subroutine

      subroutine p3dfft_3Dderiv_double(plan,IN,OUT,idir,OW)
      use p3dfft_plus_plus

      double precision, TARGET :: IN(1,1,*)
      double precision, TARGET :: OUT(1,1,*)
      integer plan,idir,OW

      call p3dfft_exec_3Dderiv_double(plan,IN,OUT,idir,OW)

      return
      end subroutine

