AM_CXXFLAGS = -I$(top_builddir)/include -include $(top_builddir)/include/config.h -O3 -DGIT_VERSION=\"$(GIT_VERSION)\" -DGIT_DATE=\"$(GIT_DATE)\"

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTWF_LIB) $(FFTW_LIB)

cppsampledir = $(datadir)/p3dfft++_samples/

cppsample_PROGRAMS = test3D_r2c_cpp.ipo test_transplan_cpp test_cos_cpp test_cos_complex_cpp test_sin_cpp 

test3D_r2c_cpp_ipo_SOURCES = test3D_r2c.C

test_transplan_cpp_SOURCES = test_transplan.C

test_cos_cpp_SOURCES = test_cos.C

test_cos_complex_cpp_SOURCES = test_cos_complex.C 

test_sin_cpp_SOURCES = test_sin.C  
