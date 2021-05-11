
GIT_DATE=2021-5-11
GIT_VERSION=0
CC = mpicc
LD = $(CC)
FFTW_LIB = -L$(FFTWHOME)/lib -lfftw3 -lfftw3f
p3dfft_topdir = ../..
CFLAGS = -I$(p3dfft_topdir)/include -DGIT_VERSION=$(GIT_VERSION) -DGIT_DATE=$(GIT_DATE)
LDFLAGS = -lm -L$(p3dfft_topdir)/build -lp3dfft.3 $(FFTW_LIB)

test3D_r2c.o: test3D_r2c.c
	      $(CC) -c test3D_r2c.c $(CFLAGS)

all: test3D_r2c.o
	$(LD) test3D_r2c.o $(LDFLAGS) -o test3D_r2c_c

	      