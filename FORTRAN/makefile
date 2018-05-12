#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O3 -DFFTW 
CC = mpicxx
CFLAGS = -O3 -DFFTW
FF = mpif90
FFLAGS = -O3
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f
LDFLAGS= $(FFT_LIB) -cxxlib  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(TACC_FFTW3_INC) -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------

.SUFFIXES: .f90 .mod

all: test3D test_trans1D 
test3D_r2c.o: test3D_r2c.f90 p3dfft++.mod
wrap.o:   wrap.f90 p3dfft++.mod
test3D: test3D_r2c.o wrap.o 
	$(FF) test3D_r2c.o wrap.o -o test3D_r2c_f -L.. -lp3dfft.3 $(LDFLAGS) 

test1D: test1D.o wrap.o
	$(FF) test1D.o wrap.o -o test1D_f -L.. -lp3dfft.3 $(LDFLAGS) 
test1D.o: test1D.f90 p3dfft++.mod

test_trans1D: test_trans1D.o wrap.o
	$(FF) test_trans1D.o wrap.o -o test_trans1D_f -L.. -lp3dfft.3 $(LDFLAGS) 
test_trans1D.o: test_trans1D.f90 p3dfft++.mod

.C.o:   
	$(CPP) -c $(CPPFLAGS) $(INCL) $<
.c.o: 
	$(CC) -c $(CFLAGS) $(INCL) *.c
.f90.o:
	$(FF) -c $(FFLAGS) $(INCL) $<
.f90.mod:
	$(FF) -c $(FFLAGS) $(INCL) $<
.F.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o:
	$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm *.o *.mod test1D_f test3D_f
