#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O3 -DFFTW -DTIMERS -xCORE-AVX2 -axCORE-AVX512,MIC-AVX512
CC = mpicxx
CFLAGS = -O3 -DFFTW
FF = mpif90
#FFLAGS = -O3 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f
LDFLAGS= $(FFT_LIB) -lm
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(TACC_FFTW3_INC) -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------


all: test3D test_transplan test_cos test_sin test_cos_complex

test3D: test3D_r2c.o 
	$(CPP) test3D_r2c.o -o test3D_r2c_cpp.ipo -L.. -lp3dfft.3 $(LDFLAGS) 
test3D_r2c.o: test3D_r2c.C ../p3dfft.h 

test_transplan: test_transplan.o 
	$(CPP) test_transplan.o -o test_transplan_cpp -L.. -lp3dfft.3 $(LDFLAGS)
test_transplan.o: test_transplan.C ../p3dfft.h

test_cos: test_cos.o 
	$(CPP) test_cos.o -o test_cos_cpp -L.. -lp3dfft.3 $(LDFLAGS)
test_cos.o: test_cos.C ../p3dfft.h

test_cos_complex: test_cos_complex.o 
	$(CPP) test_cos_complex.o -o test_cos_complex_cpp -L.. -lp3dfft.3 $(LDFLAGS)
test_cos_complex.o: test_cos_complex.C ../p3dfft.h

test_sin: test_sin.o 
	$(CPP) test_sin.o -o test_sin_cpp -L.. -lp3dfft.3 $(LDFLAGS)
test_sin.o: test_sin.C ../p3dfft.h
 

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
	/bin/rm *.o 
