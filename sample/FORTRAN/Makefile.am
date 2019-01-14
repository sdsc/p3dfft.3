AM_FCFLAGS = -I$(top_builddir)/include -O3

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTWF_LIB) $(FFTW_LIB)

AM_LDFLAGS = -cxxlib

fsampledir = $(datadir)/p3dfft++_samples/

fsample_PROGRAMS = test3D_r2c_f test3D_c2c_f test_deriv_f test_trans1D_f

p3dfft++.o: p3dfft++.f90

test3D_r2c.o: test3D_r2c.f90 p3dfft++.o

test3D_c2c.o: test3D_c2c.f90 p3dfft++.o

test_deriv.o: test_deriv.f90 p3dfft++.o

test_trans1D.o: test_trans1D.f90 p3dfft++.o

test3D_r2c_f_SOURCES = test3D_r2c.f90 wrap.f90

test3D_c2c_f_SOURCES = test3D_c2c.f90 wrap.f90

test_deriv_f_SOURCES = test_deriv.f90 wrap.f90

test_trans1D_f_SOURCES = test_trans1D.f90 wrap.f90
