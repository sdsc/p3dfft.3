#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O3 -DFFTW
CC = mpicxx
CFLAGS = -O3 -DFFTW
FF = mpif90
#FFLAGS = -O3 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(HOME)/fftw-3.3.5/lib -lfftw3 -lfftw3f
LDFLAGS= -g $(FFT_LIB)  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(FFTWHOME)/include -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------


all: test3D test_trans1D
test3D: test3D.o 
	$(CPP) test3D.o -o test3D_c -L.. -lp3dfft.3 $(LDFLAGS) 
test3D.o: test3D.c ../p3dfft.h 

test_trans1D: test_trans1D.o 
	$(CPP) test_trans1D.o -o test_trans1D_c -L.. -lp3dfft.3 $(LDFLAGS) 
test_trans1D.o: test_trans1D.c ../p3dfft.h



.C.o:   
	$(CPP) -c $(CPPFLAGS) $(INCL) $<
.c.o: 
	$(CC) -c $(CFLAGS) $(INCL) $<
.F90.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.F.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o:
	$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm *.o test1D_c test2D_c