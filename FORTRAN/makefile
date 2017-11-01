#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O0 -g -DFFTW -DDEBUG 
CC = mpicxx
CFLAGS = -O0 -g -DFFTW
FF = mpif90
FFLAGS = -g -O0
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(HOME)/fftw-3.3.5/lib -lfftw3 -lfftw3f
LDFLAGS= $(FFT_LIB) -cxxlib
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(FFTWHOME)/include -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------

.SUFFIXES: .f90 .mod

all: test2D 
test2D.o: test2D.f90 p3dfft++.mod
wrap.o:   wrap.f90 p3dfft++.mod
test2D: test2D.o wrap.o 
	$(FF) test2D.o wrap.o -o test2D_f -L.. -lp3dfft.3 $(LDFLAGS) 

test1D: test1D.o wrap.o
	$(CPP) test1D.o wrap.o -o test1D_f -L.. -lp3dfft.3 $(LDFLAGS) 
test1D.o: test1D.f90 p3dfft++.mod

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
	/bin/rm *.o *.mod test1D_f test2D_f
