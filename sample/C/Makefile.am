AM_CPPFLAGS = -I$(top_builddir)/include -include $(top_builddir)/include/config.h#-03

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTWF_LIB) $(FFTW_LIB)

AM_LDFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -ldl

csampledir = $(datadir)/p3dfft++_samples/

csample_PROGRAMS = test3D_r2c_c test3D_c2c_c test_deriv_c test_trans1D_c

test3D_r2c_c_SOURCES = test3D_r2c.c

test3D_c2c_c_SOURCES = test3D_c2c.c

test_deriv_c_SOURCES = test_deriv.c

test_trans1D_c_SOURCES = test_trans1D.c
