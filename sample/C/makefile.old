#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O3 -DFFTW
CC = mpicxx
CFLAGS = -O0 -g -DFFTW
FF = mpif90
#FFLAGS = -O3 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f
LDFLAGS= -g $(FFT_LIB) -lm
# -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(TACC_FFTW3_INC) -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------


all: test3D test1D test_deriv
test3D: test3D_r2c.o test3D_c2c.o test3D_c2c_inplace.o
	$(CPP) test3D_r2c.o -o test3D_r2c_c -L.. -lp3dfft.3 $(LDFLAGS) 
	$(CPP) test3D_c2c.o -o test3D_c2c_c -L.. -lp3dfft.3 $(LDFLAGS) 
	$(CPP) test3D_c2c_inplace.o -o test3D_c2c_inplace_c -L.. -lp3dfft.3 $(LDFLAGS) 
test3D_r2c.o: test3D_r2c.c ../p3dfft.h 
test3D_c2c.o: test3D_c2c.c ../p3dfft.h 
test3D_c2c_inplace.o: test3D_c2c_inplace.c ../p3dfft.h 

test_deriv: test_deriv.o
	$(CPP) test_deriv.o -o test_deriv_c -L.. -lp3dfft.3 $(LDFLAGS)

test1D: test1D_r2c.o test1D_cos.o test1D_cos_complex.o 
	$(CPP) test1D_r2c.o -o test1D_r2c_c -L.. -lp3dfft.3 $(LDFLAGS) 
	$(CPP) test1D_cos.o -o test1D_cos_c -L.. -lp3dfft.3 $(LDFLAGS) 
	$(CPP) test1D_cos_complex.o -o test1D_cos_complex_c -L.. -lp3dfft.3 $(LDFLAGS) 
test1D_r2c.o: test1D_r2c.c ../p3dfft.h
test1D_cos.o: test1D_cos.c ../p3dfft.h
test1D_cos_complex.o: test1D_cos_complex.c ../p3dfft.h


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
