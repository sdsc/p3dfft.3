
      Ntot = fsize(1)*fsize(2)*fsize(3)
      Nglob = nx * ny
      Nglob = Nglob * nz
      factor = 1.0d0/Nglob

      rtime1 = 0.0
!
! Repeat n times

      do  m=1,n
         if(proc_id .eq. 0) then
            print *,'Iteration ',m
         endif

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Forward transform
         call p3dfft_ftran_r2c (BEG,AEND,'fft')

         rtime1 = rtime1 + MPI_wtime()

         if(proc_id .eq. 0) then
            print *,'Result of forward transform:'
         endif
         call print_all(AEND,Ntot,proc_id,Nglob)

! normalize
         call mult_array(AEND, Ntot,factor)

! Barrier for correct timing
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         rtime1 = rtime1 - MPI_wtime()
! Backward transform
         call p3dfft_btran_c2r (AEND,C,'tff')
         rtime1 = rtime1 + MPI_wtime()

      end do

! Free work space
      call p3dfft_clean

! Check results
      call check_res(C)

! Gather timing statistics
      call MPI_Reduce(rtime1,rtime2,1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if (proc_id.eq.0) write(6,*)'proc_id, cpu time per loop', &
         proc_id,rtime2/dble(n)

      timers = timers / dble(n)

      call MPI_Reduce(timers,gt(1,1),12,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,2),12,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      call MPI_Reduce(timers,gt(1,3),12,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      tc = (timers(1)+timers(2)+timers(3)+timers(4))
      call MPI_Reduce(tc,gtcomm(1),1,mpi_real8,MPI_SUM,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(2),1,mpi_real8,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)
      call MPI_Reduce(tc,gtcomm(3),1,mpi_real8,MPI_MIN,0, &
        MPI_COMM_WORLD,ierr)

      gt(1:12,1) = gt(1:12,1) / dble(nproc)
      gtcomm(1) = gtcomm(1) / dble(nproc)

      if(proc_id .eq. 0) then
         do i=1,12
            print *,'timer',i,' (avg/max/min): ',gt(i,:)
         enddo
         print *,'Total comm (avg/max/min): ',gtcomm
      endif


      call MPI_FINALIZE (ierr)

      contains

!=========================================================
	subroutine check_res(C)
!=========================================================

	double precision C(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
	double precision cdiff,ccdiff,sinyz,ans,prec
	integer x,y,z


      cdiff=0.0d0
      do 20 z=istart(3),iend(3)
         do 20 y=istart(2),iend(2)
            sinyz=siny(y)*sinz(z)
            do 20 x=istart(1),iend(1)
            ans=sinx(x)*sinyz
            if(cdiff .lt. abs(C(x,y,z)-ans)) then
               cdiff = abs(C(x,y,z)-ans)
!               print *,'x,y,z,cdiff=',x,y,z,cdiff
            endif
 20   continue
      call MPI_Reduce(cdiff,ccdiff,1,p3dfft_mpireal,MPI_MAX,0, &
        MPI_COMM_WORLD,ierr)

      if(proc_id .eq. 0) then
         if(p3dfft_type .eq. 8) then
            prec = 1e-14
         else
            prec = 1e-5
         endif
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
	subroutine init_ar_sine(A)
!=========================================================

	double precision A(istart(1):iend(1),istart(2):iend(2),istart(3):iend(3))
	integer x,y,z
	double precision sinyz

      allocate (sinx(nx))
      allocate (siny(ny))
      allocate (sinz(nz))

      do z=istart(3),iend(3)
         sinz(z)=sin((z-1)*twopi/nz)
      enddo
      do y=istart(2),iend(2)
         siny(y)=sin((y-1)*twopi/ny)
      enddo
      do x=istart(1),iend(1)
         sinx(x)=sin((x-1)*twopi/nx)
      enddo

! Initialize with 3D sine wave

      do z=istart(3),iend(3)
         do y=istart(2),iend(2)
            sinyz=siny(y)*sinz(z)
            do x=istart(1),iend(1)
               A(x,y,z)=sinx(x)*sinyz
            enddo
         enddo
      enddo

      return
      end subroutine


!=========================================================
      subroutine mult_array(X,nar,f)

      use p3dfft

      integer(i8) nar,i
      complex(p3dfft_type) X(nar)
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

      subroutine print_all(Ar,Nar,proc_id,Nglob)

      use p3dfft

      integer x,y,z,proc_id
      integer(i8) i,Nar
      complex(p3dfft_type) Ar(1,1,*)
      integer Fstart(3),Fend(3),Fsize(3)
      real(r8) Nglob

      call p3dfft_get_dims(Fstart,Fend,Fsize,2)

      do i=1,Nar
         if(abs(Ar(1,1,i)) .gt. Nglob *1.25e-4) then
            z = (i-1)/(Fsize(1)*Fsize(2))
            y = (i-1 - z * Fsize(1)*Fsize(2))/Fsize(1)
            x = i-1-z*Fsize(1)*Fsize(2) - y*Fsize(1)
            print *,'(',x+Fstart(1),y+Fstart(2),z+Fstart(3),') ',Ar(1,1,i)
         endif
      enddo

      return
      end subroutine
