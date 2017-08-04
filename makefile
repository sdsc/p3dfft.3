#SHELL=/bin/csh

CPP = mpicxx
CPPFLAGS = -O0 -g -DFFTW -DDEBUG 
CC = mpicc
FF = mpif90
#FFLAGS = -O3 
AR = ar
ARFLAGS = -v -r -u
FFT_LIB = -L$(FFTWHOME)/lib -lfftw3 -lfftw3f
LDFLAGS= $(FFT_LIB) -lm
# For FFTW use path to the installed FFTW library:
# -L/usr/local/apps/fftw301s/lib -lfftw3f 

INCL = -I$(FFTWHOME)/include 
# For FFTW add include file location, for example: 
# INCL = -I/usr/local/apps/fftw312s/include 

DFLAGS = -DFFTW

P3DFFT_ROOT = .
P3DFFT_LIB = libp3dfft.3.a

# ----------------

FFT3DLIB = init.o plan.o exec.o templ.o

all: lib test
lib: $(FFT3DLIB)
	$(AR) $(ARFLAGS) $(P3DFFT_LIB) $(FFT3DLIB)	
test: $(FFT3DLIB) test1.o 
	$(CPP) test1.o -L. -lp3dfft.3 $(LDFLAGS) 
install: 
	if(!(-e $(P3DFFT_ROOT))) mkdir $(P3DFFT_ROOT)
	if (!(-e $(P3DFFT_ROOT)/lib)) mkdir $(P3DFFT_ROOT)/lib	
	cp $(P3DFFT_LIB) $(P3DFFT_ROOT)/lib
	if (!(-e $(P3DFFT_ROOT)/include)) mkdir $(P3DFFT_ROOT)/include
	cp p3dfft.mod $(P3DFFT_ROOT)/include


init.o: init.C p3dfft.h
plan.o: plan.C templ.C p3dfft.h
exec.o: exec.C p3dfft.h
templ.o: templ.C p3dfft.h
test1.o: test1.C p3dfft.h templ.o

.C.o:   
	$(CPP) -c $(CPPFLAGS) $(INCL) $<
.c.o: 
	$(CC) -c $(CFLAGS) *.c
.F90.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.F.o:
	$(FF) $(DFLAGS) -c $(FFLAGS) $(INCL) $<
.f.o:
	$(FF) -c $(FFLAGS) $(INCL) $<
clean:
	/bin/rm $(FFT3DLIB)