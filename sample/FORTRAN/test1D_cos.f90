!
!This program exemplifies the use of 1D transforms in P3DFFT++, for a 1D cosine transform, for real-valued arrays. 1D transforms are performed on 3D arrays, in the dimension specified as an 
!argument. This could be an isolated 1D transform or a stage in a multidimensional transform. This function can do local transposition, i.e. arbitrary input and output memory ordering. 
!However it does not do an inter-processor transpose (see test_transMPI for that). 
!
! This program initializes a 3D array with a 3D sine wave, then
! performs cosine transform twice, and checks that
! the results are correct, namely the same as in the start except
! for a normalization factor. It can be used both as a correctness
! test and for timing the library functions.
!
! The program expects 'trans.in' file in the working directory, with
! a single line of numbers : Nx,Ny,Nz,dim,Nrep,MOIN(1)-(3),MOOUT(1)-(3). 
! Here 
!    Nx,Ny,Nz are 3D grid dimensions (note: dimension fo cosine transform must be odd)
!    dim is the dimension of 1D transform (valid values are 0 through 2, and the logical dimension si specified, i.e. actual storage dimension may be different as specified by MOIN mapping)
!    Nrep is the number of repititions. 
!    MOIN are 3 values for the memory order of the input grid, valid values of each is 0 - 2, not repeating. 
!    MOOUT is the memory order of the output grid. 
!
! Optionally a file named 'dims' can also be provided to guide in the choice
! of processor geometry in case of 2D decomposition. It should contain
! two numbers in a line, with their product equal to the total number
! of tasks. Otherwise processor grid geometry is chosen automatically.
! For better performance, experiment with this setting, varying
! iproc and jproc. In many cases, minimizing iproc gives best results.
! Setting it to 1 corresponds to one-dimensional decomposition.
!
! If you have questions please contact Dmitry Pekurovsky, dmitry@sdsc.edu

      program fft1d

      use p3dfft_plus_plus
      implicit none
      include 'mpif.h'

      double precision, dimension(:,:,:),  allocatable :: BEG,C,AEND

      integer i,n,nx,ny,nz
      integer m,x,y,z
      integer fstatus

      integer(8) Ntot
      double precision factor,rtime1,rtime2
      double precision gt(12,3),gtcomm(3),tc
      logical iex
      integer ierr,dim,dims(2),nproc,proc_id,cnt
      integer type_ids1,type_ids2,trans_f,trans_b,pdims(3),glob_start1(3),glob_start2(3)
      integer(8) size1,size2
      integer(C_INT) ldims1(3),ldims2(3),mem_order1(3),mem_order2(3),dmap(3),pgrid,gdims(3)
      integer(C_INT) grid1,grid2,workspace
      integer mpicomm,myid
      integer mydims1(3),mydims2(3),ar_dim,ar_dim2,gstart2(3)

      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD,nproc,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,proc_id,ierr)

! Read input parameters from 'stdin' file

      if (proc_id.eq.0) then
         open (unit=3,file='trans.in',status='old', &
               access='sequential',form='formatted', iostat=fstatus)
         if (fstatus .eq. 0) then
            write(*, *) ' Reading from input file trans.in'
         endif
         dim = 1

        read (3,*) nx, ny, nz, dim,n, mem_order1(1:3), mem_order2(1:3)
        dim = dim +1
	print *,'P3DFFT++ test, Fortran, 1D cosine transform, 3D wave input'
        write (*,*) "procs=",nproc," nx=",nx, &
                " ny=", ny," nz=", nz,"dim=",dim," repeat=", n
        write (*,*) "mem_order on input: ",mem_order1
        write (*,*) "mem_order on output: ",mem_order2
       endif

! Broadcast parameters

      call MPI_Bcast(nx,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(ny,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(nz,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(n,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(dim,1, MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(mem_order1,3,MPI_INTEGER,0,mpi_comm_world,ierr)
      call MPI_Bcast(mem_order2,3,MPI_INTEGER,0,mpi_comm_world,ierr)

! Establish 2D processor grid decomposition, either by readin from file 'dims' or by an MPI default

!    nproc is devided into a iproc x jproc stencil
!

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

      if(proc_id .eq. 0) then
         print *,'Using processor grid ',dims(1),' x ',dims(2)
      endif

! Set up work structures for P3DFFT

      call p3dfft_setup

! Set up 2 transform types for 1D transforms

      type_ids1 = P3DFFT_DCT1_REAL_D;

      type_ids2 = P3DFFT_DCT1_REAL_D;

! Set up global dimensions of the grid

      gdims(1)=nx
      gdims(2)=ny
      gdims(3)=nz

! Set up processor order and memory ordering, as well as the final global grid dimensions (these will be different from the original dimensions in one dimension since we are doing real-to-complex transform)

      factor = 0.5d0 /(gdims(dim)-1.0d0)

      ar_dim = mem_order1(dim)+1
      ar_dim2 = mem_order2(dim)+1

! Define processor grid. Make the direction of transform local

      pdims(1) = 1
      pdims(2) = dims(1)
      pdims(3) = dims(2)

      cnt =1
      do i=1,3
         if(i .eq. dim) then
            dmap(i)=0
         else
            dmap(i) = cnt
            cnt = cnt+1
         endif
      enddo

! Specify the default communicator for P3DFFT++. This can be different from your program default communicator if you wish to keep P3DFFT++ communications separate from yours

      pgrid = p3dfft_init_proc_grid(pdims,MPI_COMM_WORLD)

! Initialize initial and final grids, based on the above information
! No conjugate symmetry (-1)
      call p3dfft_init_data_grid(grid1,ldims1, glob_start1,gdims,-1,pgrid,dmap,mem_order1)

      call p3dfft_init_data_grid(grid2,ldims2,glob_start2,gdims,-1,pgrid,dmap,mem_order2)

! Set up the forward transform, based on the predefined 3D transform type and grid1 and grid2. This is the planning stage, needed once as initialization.

      call p3dfft_plan_1Dtrans(trans_f,grid1,grid2,type_ids1,dim-1,workspace)

! Now set up the backward transform
      call p3dfft_plan_1Dtrans(trans_b,grid2,grid1,type_ids2,dim-1,workspace)

! Determine local array dimensions. These are defined taking into account memory ordering. 

      do i=1,3
         mydims1(mem_order1(i)+1) = ldims1(i)
         mydims2(mem_order2(i)+1) = ldims2(i)
         gstart2(mem_order2(i)+1) = glob_start2(i)
      enddo

! Now allocate initial and final arrays in physical space (real-valued)
      allocate(BEG(mydims1(1),mydims1(2),mydims1(3)))
      allocate(C(mydims1(1),mydims1(2),mydims1(3)))

! Initialize the BEG array with a sine wave in 3D

      call init_wave(BEG,mydims1,ar_dim)

! Now allocate the array for holding Fourier space data

      allocate(AEND(mydims2(1),mydims2(2),mydims2(3)))

! Warm-up call to execute forward 3D FFT transform
      call p3dfft_1Dtrans_double(trans_f,BEG,AEND,0)

      Ntot = ldims2(1)*ldims2(2)
      Ntot = Ntot*ldims2(3)
      rtime1 = 0.0

! Start the timing loop
      do m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Forward transform
         call p3dfft_1Dtrans_double(trans_f,BEG,AEND,0)

         rtime1 = rtime1 + MPI_wtime()

!         if(proc_id .eq. 0) then
!            print *,'Result of forward transform:'
!         endif
!         call print_all(AEND,mydims2,glob_start2,mem_order2,ar_dim2)

! normalize
         call mult_array(AEND, Ntot,factor)

      call check_res_forward(AEND,mydims2,ar_dim2,gstart2)

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Backward transform
         call p3dfft_1Dtrans_double(trans_b,AEND,C,1)
         rtime1 = rtime1 + MPI_wtime()

      enddo

! Free work space
      
      call p3dfft_cleanup


! Check results
      call check_res(BEG,C,mydims1,ar_dim)
! Gather timing statistics
      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop', &
         proc_id,rtime2/dble(n)

      call MPI_FINALIZE (ierr)

    end program fft1d

    subroutine intcpy(in,out,n)

      integer in(3),out(3),n,i

      do i=1,n
         out(i) = in(i)
      enddo
      return
    end subroutine intcpy

!=========================================================
	subroutine check_res_forward(A,mydims,dim,globstart)
!=========================================================

        implicit none
        include 'mpif.h'

        integer dim,mydims(3),globstart(3),myid,ierr
        double precision A(mydims(1),mydims(2),mydims(3)),ans
        integer x,y,z
        double precision d,cdiff,ccdiff,prec

        cdiff = 0.0
        ccdiff = 0.0

        do z=1,mydims(3)
        do y=1,mydims(2)
        do x=1,mydims(1)
           if(dim .eq. 1) then
              if(x+globstart(dim) .eq. 2) then
                 ans = 0.5
              else
                 ans = 0.0
              endif
           else if(dim .eq. 2) then
              if(y+globstart(dim) .eq. 2) then
                 ans = 0.5
              else
                 ans = 0.0
              endif
           else if(dim .eq. 3) then
              if(z+globstart(dim) .eq. 2) then
                 ans = 0.5
              else
                 ans = 0.0
              endif
           endif


           d = abs(A(x,y,z) - ans) 
           if(cdiff .lt. d) then
              cdiff = d
           endif
           enddo
        enddo
     enddo

      call MPI_Reduce(cdiff,ccdiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      if(myid .eq. 0) then
         prec = 1e-14
         if(ccdiff .gt. prec * mydims(dim) *0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
         write (6,*) 'max diff =',ccdiff
      endif

      return
      end subroutine

!=========================================================
	subroutine check_res(A,B,mydims,ar_dim)
!=========================================================

        implicit none
        include 'mpif.h'

        integer gdims(3),mydims(3),glob_start(3)
	double precision A(mydims(1),mydims(2),mydims(3))
	double precision B(mydims(1),mydims(2),mydims(3))
	double precision cdiff,ccdiff,sinyz,ans,prec
	integer x,y,z,ierr,myid,ar_dim
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)

      cdiff=0.0d0
      do 20 z=1,mydims(3)
         do 20 y=1,mydims(2)
            do 20 x=1,mydims(1)
            if(cdiff .lt. abs(A(x,y,z)-B(x,y,z))) then
               cdiff = abs(A(x,y,z)-B(x,y,z))
!               print *,'x,y,z,cdiff=',x,y,z,cdiff
            endif
 20   continue
            
      print *,myid,': My diff=',cdiff
      
      call MPI_Reduce(cdiff,ccdiff,1,MPI_DOUBLE_PRECISION,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if(myid .eq. 0) then
         prec = 1e-14
         if(ccdiff .gt. prec * mydims(ar_dim) *0.25) then
            print *,'Results are incorrect'
         else
            print *,'Results are correct'
         endif
         write (6,*) 'max diff =',ccdiff
      endif

      return
      end subroutine

!=========================================================
	subroutine init_wave(A,mydims,ar_dim)
!=========================================================
   
        implicit none
   
        integer mydims(3)
	double precision A(mydims(1),mydims(2),mydims(3))
	integer x,y,z,ar_dim
	double precision pi,cos_coords(mydims(ar_dim))

        pi=atan(1.0d0)*4.0d0
        
        do z=1,mydims(ar_dim)
           cos_coords(z)=cos((z-1)*pi/(mydims(ar_dim)-1.0d0))
        enddo

! Initialize with 3D sine wave
      
      if(ar_dim .eq. 1) then

         do z=1,mydims(3)
            do y=1,mydims(2)
               do x=1,mydims(1)
                  A(x,y,z)=cos_coords(x)
               enddo
            enddo
         enddo

      else if(ar_dim .eq. 2) then

         do z=1,mydims(3)
            do y=1,mydims(2)
               do x=1,mydims(1)
                  A(x,y,z)=cos_coords(y)
               enddo
            enddo
         enddo

      else

         do z=1,mydims(3)
            do y=1,mydims(2)
               do x=1,mydims(1)
                  A(x,y,z)=cos_coords(z)
               enddo
            enddo
         enddo

      endif

      return
      end subroutine


!=========================================================
      subroutine mult_array(X,nar,f)

      implicit none
 
      integer(8) nar,i
      double precision X(nar)
      double precision f

      do i=1,nar
         X(i) = X(i) * f
      enddo

      return
      end subroutine

!=========================================================
! Translate one-dimensional index into three dimensions,
!    print out significantly non-zero values
!

      subroutine print_all(Ar,mydims,gstart,mo,ar_dim)

      implicit none
 
      integer x,y,z,proc_id,mydims(3),gstart(3),ar_dim,mo(3),imo(3)
      integer(8) i,Nar
      double precision Ar(mydims(1),mydims(2),mydims(3))
      real(8) Nglob

      do i=1,3
         imo(mo(i)+1)=i
      enddo

      do z=1,mydims(3)
         do y=1,mydims(2)
            do x=1,mydims(1)
               if(abs(Ar(x,y,z)) .gt. mydims(ar_dim) *1.25e-4) then 
                  print *,'(',x+gstart(imo(1)),y+gstart(imo(2)),z+gstart(imo(3)),') ',Ar(x,y,z)
               endif
            enddo
         enddo
      enddo

      return
      end subroutine
