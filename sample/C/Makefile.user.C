# Sample makefile for user programs in C

# C compiler
CC = mpicc
# Loader/linker
LD = $(CC)
# FFTW library flags
FFTW_LIB = -L$(FFTWHOME)/lib -lfftw3 -lfftw3f
# Top level directory where P3DFFT++ is installed
p3dfft_topdir = ../../../install
# Compiler flags
CFLAGS = -I$(p3dfft_topdir)/include -DGIT_VERSION=$(GIT_VERSION) -DGIT_DATE=$(GIT_DATE)
# Linker flags
LDFLAGS = -lm -L$(p3dfft_topdir)/lib -lp3dfft.3 $(FFTW_LIB)
# Ignore these two variables - they are needed for the package examples
GIT_DATE='"2021-5-11"'
GIT_VERSION='"0"'

test3D_r2c.o: test3D_r2c.c
	      $(CC) -c test3D_r2c.c $(CFLAGS)

all: test3D_r2c.o
	$(LD) test3D_r2c.o $(LDFLAGS) -o test3D_r2c_c

	      
