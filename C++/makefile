#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O0 -g -DFFTW -DDEBUG 
CC = mpicxx
CFLAGS = -O0 -g -DFFTW
FF = mpif90
#FFLAGS = -O3 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(HOME)/fftw-3.3.5/lib -lfftw3 -lfftw3f
LDFLAGS= $(FFT_LIB) -lm
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(FFTWHOME)/include -I.. 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------


all: test2D test1D
test2D: test2D.o 
	$(CPP) test2D.o -o test2D_cpp -L.. -lp3dfft.3 $(LDFLAGS) 
test2D.o: test2D.C ../p3dfft.h 

test1D: test1D.o 
	$(CPP) test1D.o -o test1D_cpp -L.. -lp3dfft.3 $(LDFLAGS) 
test1D.o: test1D.C ../p3dfft.h

.C.o:   
	$(CPP) -c $(CPPFLAGS) $(INCL) $<
.c.o: 
	$(CC) -c $(CFLAGS) $(INCL) *.c
.F90.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.F.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o:
	$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm *.o test1D_cpp test2D_cpp