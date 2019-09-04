AM_FCFLAGS = -I$(top_builddir)/include $(FFTW_INC) $(ESSL_INC)

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTW_LIB) $(FFTWF_LIB) $(ESSL_LIB) -lstdc++

AM_LDFLAGS = #-cxxlib

fsampledir = $(datadir)/p3dfft++_samples/

fsample_PROGRAMS = test1D_cos_f test1D_cos_complex_f test1D_r2c_f test3D_c2c_f test3D_r2c_f  testIDIR_deriv_f test3D_c2c_inplace_f testIDIR_deriv2_f# test2D_f test3D_f

p3dfft_plus_plus.mod: p3dfft++.o

wrap.o: wrap.f90 p3dfft_plus_plus.mod

test1D_cos.o: test1D_cos.f90 p3dfft_plus_plus.mod

test1D_cos_complex.o: test1D_cos_complex.f90 p3dfft_plus_plus.mod

test1D_r2c.o: test1D_r2c.f90 p3dfft_plus_plus.mod

#test2D.o: test2D.f90 p3dfft++.o

#test3D.o: test3D.f90 p3dfft++.o

test3D_c2c_inplace.o: test3D_c2c_inplace.f90 p3dfft_plus_plus.mod

test3D_r2c.o: test3D_r2c.f90 p3dfft_plus_plus.mod

test3D_c2c.o: test3D_c2c.f90 p3dfft_plus_plus.mod

test_deriv.o: test_deriv.f90 p3dfft_plus_plus.mod

test_deriv2.o: test_deriv2.f90 p3dfft_plus_plus.mod

test1D_cos_f_SOURCES = test1D_cos.f90 wrap.f90 p3dfft++.f90

test1D_cos_complex_f_SOURCES = test1D_cos_complex.f90 wrap.f90 p3dfft++.f90

test1D_r2c_f_SOURCES = test1D_r2c.f90 wrap.f90 p3dfft++.f90

#test2D_f_SOURCES = test2D.f90

#test3D_f_SOURCES = test3D.f90 wrap.f90

test3D_c2c_inplace_f_SOURCES = test3D_c2c_inplace.f90 wrap.f90 p3dfft++.f90

test3D_r2c_f_SOURCES = test3D_r2c.f90 wrap.f90 p3dfft++.f90

test3D_c2c_f_SOURCES = test3D_c2c.f90 wrap.f90 p3dfft++.f90

testIDIR_deriv_f_SOURCES = test_deriv.f90 wrap.f90 p3dfft++.f90

testIDIR_deriv2_f_SOURCES = test_deriv2.f90 wrap.f90 p3dfft++.f90

clean-local:
	-test -z "*.mod" || rm -f *.mod
