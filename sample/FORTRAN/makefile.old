#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O3 -DFFTW 
CC = mpicxx
CFLAGS = -O3 -DFFTW
FF = mpif90
FFLAGS = -O0 -g -fno-stack-protector 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(TACC_FFTW3_LIB) -lfftw3 -lfftw3f
LDFLAGS= -g $(FFT_LIB) -cxxlib  
#-Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
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

all: test3D test1D test_deriv
test3D_r2c.o: test3D_r2c.f90 p3dfft_plus_plus.mod
test3D_c2c.o: test3D_c2c.f90 p3dfft_plus_plus.mod
test_deriv.o: test_deriv.f90 p3dfft_plus_plus.mod
p3dfft_plus_plus.mod: p3dfft++.o 
wrap.o:   wrap.f90 p3dfft_plus_plus.mod
test3D: test3D_r2c.o test3D_c2c.o test3D_c2c_inplace.o wrap.o 
	$(FF) test3D_r2c.o wrap.o -o test3D_r2c_f -L.. -lp3dfft.3 $(LDFLAGS) 
	$(FF) test3D_c2c.o wrap.o -o test3D_c2c_f -L.. -lp3dfft.3 $(LDFLAGS) 
	$(FF) test3D_c2c_inplace.o wrap.o -o test3D_c2c_inplace_f -L.. -lp3dfft.3 $(LDFLAGS) 

test3D_r2c.o: test3D_r2c.f90 wrap.f90 p3dfft++.f90
test3D_c2c.o: test3D_c2c.f90 wrap.f90 p3dfft++.f90
test3D_c2c_inplace.o: test3D_c2c_inplace.f90 wrap.f90 p3dfft++.f90

test_deriv: test_deriv.o wrap.o
	$(FF) test_deriv.o wrap.o -o test_deriv_f -L.. -lp3dfft.3 $(LDFLAGS) 

test1D: test1D_r2c.o test1D_cos.o test1D_cos_complex.o wrap.o
	$(FF) test1D_r2c.o wrap.o -o test1D_r2c_f -L.. -lp3dfft.3 $(LDFLAGS) 
	$(FF) test1D_cos.o wrap.o -o test1D_cos_f -L.. -lp3dfft.3 $(LDFLAGS) 
	$(FF) test1D_cos_complex.o wrap.o -o test1D_cos_complex_f -L.. -lp3dfft.3 $(LDFLAGS) 
test1D_r2c.o: test1D_r2c.f90 p3dfft_plus_plus.mod
test1D_cos.o: test1D_cos.f90 p3dfft_plus_plus.mod
test1D_cos_complex.o: test1D_cos_complex.f90 p3dfft_plus_plus.mod

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
