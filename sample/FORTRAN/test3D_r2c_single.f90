! This sample program illustrates the
! use of P3DFFT++ library for highly scalable parallel 3D FFT,
! for a 3D real-to-complex FFT.
!
! This program initializes a 3D array with a 3D sine wave, then
! performs forward transform, backward transform, and checks that
! the results are correct, namely the same as in the start except
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions.
!
! The program expects 'stdin' file in the working directory, with
! a single line of numbers : Nx,Ny,Nz,Ndim,Nrep. Here Nx,Ny,Nz
! are box dimensions, Ndim is the dimentionality of processor grid
! (1 or 2), and Nrep is the number of repititions. Optionally
! a file named 'dims' can also be provided to guide in the choice
! of processor geometry in case of 2D decomposition. It should contain
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying
! iproc and jproc. In many cases, minimizing iproc gives best results.
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program fft3d

      use p3dfft_plus_plus
      implicit none
      include 'mpif.h'

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus
      logical flg_inplace

      real, dimension(:,:,:),  allocatable :: BEG,C
      complex, dimension(:,:,:),  allocatable :: AEND

      integer(8) Ntot
      real factor
      real(8) rtime1,rtime2
      real Nglob,prec
      real(8) gt(12,3),gtcomm(3),tc
      integer ierr,nu,ndim,dims(2),nproc,proc_id
      integer istart(3),iend(3),isize(3)
      integer fstart(3),fend(3),fsize(3)
      integer iproc,jproc,nxc,nyc,nzc
      logical iex
      integer type_ids1(3),type_ids2(3),trans_f,trans_b,pdims(2)
      integer type_rcc,type_ccr,glob_start1(3),glob_start2(3)
      integer gstart1(3),gstart2(3)
      integer(8) size1,size2
      integer(C_INT) ldims1(3),ldims2(3),mem_order1(3),mem_order2(3),proc_order(3),pgrid1(3),pgrid2(3),gdims1(3),gdims2(3)
      integer(C_INT) grid1,grid2
      integer mpicomm,myid
      integer a(3)
      integer mydims1(3),mydims2(3)

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)


      gt=0.0
      gtcomm=0.0

! Read input parameters from 'stdin' file

      if (proc_id.eq.0) then
         open (unit=3,file='stdin',status='old', &
               access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file stdin'
         endif
         ndim = 2

        read (3,*) nx, ny, nz, ndim,n
	print *,'P3DFFT test, Single Precision, 3D wave input, 3D Real-to-complex FFT'
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"ndim=",ndim," repeat=", n
       endif

! Broadcast parameters

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ndim,1, MPI_INTEGER,0,mpi_comm_world,ierr)

! Establish 2D processor grid decomposition, either by reading from file 'dims' or by an MPI default

!    nproc is devided into a iproc x jproc stencil
!

      if(ndim .eq. 1) then
         dims(1) = 1
         dims(2) = nproc
      else if(ndim .eq. 2) then
	inquire(file='dims',exist=iex)
	if (iex) then
           if (proc_id.eq.0) print *, 'Reading proc. grid from file dims'
           open (999,file='dims')
           read (999,*) dims(1), dims(2)
           close (999)
           if(dims(1) * dims(2) .ne. nproc) then
              dims(2) = nproc / dims(1)
           endif
	else
           if (proc_id.eq.0) print *, 'Creating proc. grid with mpi_dims_create'
           dims(1) = 0
           dims(2) = 0
           call MPI_Dims_create(nproc,2,dims,ierr)
           if(dims(1) .gt. dims(2)) then
              dims(1) = dims(2)
              dims(2) = nproc / dims(1)
           endif
        endif
      endif

      iproc = dims(1)
      jproc = dims(2)

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',iproc,' x ',jproc
      endif

! Set up work structures for P3DFFT

      call p3dfft_setup

! Set up 2 transform types for 3D transforms

      type_ids1(1) = P3DFFT_R2CFFT_S;
      type_ids1(2) = P3DFFT_CFFT_FORWARD_S;
      type_ids1(3) = P3DFFT_CFFT_FORWARD_S;

      type_ids2(1) = P3DFFT_C2RFFT_S;
      type_ids2(2) = P3DFFT_CFFT_BACKWARD_S;
      type_ids2(3) = P3DFFT_CFFT_BACKWARD_S;

! Now initialize 3D transforms (forward and backward) with these types
      call p3dfft_init_3Dtype(type_rcc,type_ids1)
      call p3dfft_init_3Dtype(type_ccr,type_ids2)

! Set up global dimensions of the grid

      gdims1(1)=nx
      gdims1(2)=ny
      gdims1(3)=nz

! Set up processor order and memory ordering, as well as the final global grid dimensions (these will be different from the original dimensions in one dimension due to conjugate symmetry, since we are doing real-to-complex transform)

      do i=1,3
         proc_order(i) = i-1
         mem_order1(i) = i-1
         gdims2(i) = gdims1(i)
      enddo
      gdims2(1) = gdims2(1)/2+1

! Set up memory order for the final grid layout (for complex array in Fourier space). It is more convenient to have the storage order of the array reversed, this helps save on memory access 
!bandwidth, and shouldn't affect the operations in the Fourier space very much, requiring basically a change in the loop order. However it is possible to define the memory ordering the same 
!as default (0,1,2). Note that the memory ordering is specified in C indeices, i.e. starting from 0

      mem_order2(1) = 2
      mem_order2(2) = 1
      mem_order2(3) = 0

! Define the initial processor grid. In this case, it's a 2D pencil, with 1st dimension local and the 2nd and 3rd split by iproc and jproc tasks respectively

      pgrid1(1) = 1
      pgrid1(2) = iproc
      pgrid1(3) = jproc

! Set up memory order for the final grid layout (for complex array in Fourier space). It is more convenient to have the storage order of the array reversed, this helps save on memory access 
!bandwidth, and shouldn't affect the operations in the Fourier space very much, requiring basically a change in the loop order. However, note that as an alternative, it is possible to 
!define the memory ordering the same as default (0,1,2). Note that the memory ordering is specified in C indices, i.e. starting from 0.

      pgrid2(1) = iproc
      pgrid2(2) = jproc
      pgrid2(3) = 1

! Specify the default communicator for P3DFFT++. This can be different from your program default communicator if you wish to keep P3DFFT++ communications separate from yours

      mpicomm = MPI_COMM_WORLD

! Initialize initial and final grids, based on the above information
! Initial grid has no conjugate symmetry (-1)
      call p3dfft_init_grid(grid1,ldims1, glob_start1,gdims1,-1,pgrid1,proc_order,mem_order1,MPI_COMM_WORLD)
! Final grid has conjugate symmetry in X dimension (0)
      call p3dfft_init_grid(grid2,ldims2,glob_start2,gdims2,0,pgrid2,proc_order,mem_order2,MPI_COMM_WORLD)

! Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization.

      call p3dfft_plan_3Dtrans(trans_f,grid1,grid2,type_rcc)

! Now set up the backward transform
      call p3dfft_plan_3Dtrans(trans_b,grid2,grid1,type_ccr)

! Determine local array dimensions. These are defined taking into account memory ordering. 

      do i=1,3
         mydims1(mem_order1(i)+1) = ldims1(i)
         mydims2(mem_order2(i)+1) = ldims2(i)
         gstart1(mem_order1(i)+1) = glob_start1(i)
         gstart2(mem_order2(i)+1) = glob_start2(i)
      enddo

! Now allocate initial and final arrays in physical space (real-valued)
      allocate(BEG(mydims1(1),mydims1(2),mydims1(3)))
      allocate(C(mydims1(1),mydims1(2),mydims1(3)))

! Initialize the BEG array with a sine wave in 3D

      call init_wave(BEG,gdims1,mydims1,gstart1)

! Now allocate the complex array for holding Fourier space data

      allocate(AEND(mydims2(1),mydims2(2),mydims2(3)))

! Warm-up call to execute forward 3D FFT transform
      call p3dfft_3Dtrans_single(trans_f,BEG,AEND,0)

      Ntot = ldims2(1)*ldims2(2)*ldims2(3)
      Nglob = nx * ny
      Nglob = Nglob * nz
      factor = 1.0d0/Nglob

      rtime1 = 0.0

! Start the timing loop

      do  m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Forward transform
         call p3dfft_3Dtrans_single(trans_f,BEG,AEND,0)

         rtime1 = rtime1 + MPI_wtime()

         if(proc_id .eq. 0) then
            print *,'Result of forward transform:'
         endif
         call print_all(AEND,Ntot,proc_id,Nglob,mydims2,gstart2)

! normalize
         call mult_array(AEND, Ntot,factor)

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Backward transform
         call p3dfft_3Dtrans_single(trans_b,AEND,C,1)
         rtime1 = rtime1 + MPI_wtime()

      end do

! Free work space
!      call p3dfft_cleanup

! Check results
      call check_res(C,gdims1,mydims1,gstart1,Nglob)

! Gather timing statistics
      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop', &
         proc_id,rtime2/dble(n)

      call MPI_FINALIZE (ierr)

    end program fft3d

    subroutine intcpy(in,out,n)

      integer in(3),out(3),n,i

      do i=1,n
         out(i) = in(i)
      enddo
      return
    end subroutine intcpy

!=========================================================
	subroutine check_res(C,gdims,ldims,glob_start,Nglob)
!=========================================================

        implicit none
        include 'mpif.h'

        integer gdims(3),ldims(3),glob_start(3)
	real C(glob_start(1)+1:glob_start(1)+ldims(1), &
     glob_start(2)+1:glob_start(2)+ldims(2),glob_start(3)+1:glob_start(3)+ldims(3))
	real cdiff,ccdiff,sinyz,ans,prec,twopi,Nglob
	integer x,y,z,ierr,myid
        real sinx(gdims(1))
        real siny(gdims(2))
        real sinz(gdims(3))
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)


      twopi=atan(1.0)*8.0

      do z=1,gdims(3)
         sinz(z)=sin((z-1)*twopi/gdims(3))
      enddo
      do y=1,gdims(2)
         siny(y)=sin((y-1)*twopi/gdims(2))
      enddo
      do x=1,gdims(1)
         sinx(x)=sin((x-1)*twopi/gdims(1))
      enddo

      cdiff=0.0d0
      do 20 z=glob_start(3)+1,glob_start(3)+ldims(3)
         do 20 y=glob_start(2)+1,glob_start(2)+ldims(2)
            sinyz=siny(y)*sinz(z)
            do 20 x=glob_start(1)+1,glob_start(1)+ldims(1)
            ans=sinx(x)*sinyz
            if(cdiff .lt. abs(C(x,y,z)-ans)) then
               cdiff = abs(C(x,y,z)-ans)
!               print *,'x,y,z,cdiff=',x,y,z,cdiff
            endif
 20   continue
      call MPI_Reduce(cdiff,ccdiff,1,MPI_REAL,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if(myid .eq. 0) then
         prec = 1e-6
         if(ccdiff .gt. prec * Nglob*0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
         write (6,*) 'max diff =',ccdiff
      endif

      return
      end subroutine

!=========================================================
	subroutine init_wave(A,gdims,ldims,glob_start)
!=========================================================

        implicit none

        integer gdims(3),ldims(3),glob_start(3)
	real A(glob_start(1)+1:glob_start(1)+ldims(1), &
        glob_start(2)+1:glob_start(2)+ldims(2),glob_start(3)+1:glob_start(3)+ldims(3))
	integer x,y,z
	real sinyz,twopi
        real sinx(gdims(1))
        real siny(gdims(2))
        real sinz(gdims(3))

      twopi=atan(1.0)*8.0

      do z=1,gdims(3)
         sinz(z)=sin((z-1)*twopi/gdims(3))
      enddo
      do y=1,gdims(2)
         siny(y)=sin((y-1)*twopi/gdims(2))
      enddo
      do x=1,gdims(1)
         sinx(x)=sin((x-1)*twopi/gdims(1))
      enddo

! Initialize with 3D sine wave

      do z=glob_start(3)+1,glob_start(3)+ldims(3)
         do y=glob_start(2)+1,glob_start(2)+ldims(2)
            sinyz=siny(y)*sinz(z)
            do x=glob_start(1)+1,glob_start(1)+ldims(1)
               A(x,y,z)=sinx(x)*sinyz
            enddo
         enddo
      enddo

      return
      end subroutine


!=========================================================
      subroutine mult_array(X,nar,f)

      implicit none
 
      integer(8) nar,i
      complex X(nar)
      real f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
      end subroutine

!=========================================================
! Translate one-dimensional index into three dimensions,
!    print out significantly non-zero values
!

      subroutine print_all(Ar,Nar,proc_id,Nglob,ldims,gstart)

      implicit none
 
      integer x,y,z,proc_id,ldims(3),gstart(3)
      integer(8) i,Nar
      complex Ar(ldims(1),ldims(2),ldims(3))
      real Nglob

      do z=1,ldims(3)
         do y=1,ldims(2)
            do x=1,ldims(1)
               if(abs(Ar(x,y,z)) .gt. Nglob *1.25e-4) then 
                  print *,'(',x+gstart(1),y+gstart(2),z+gstart(3),') ',Ar(x,y,z)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine
