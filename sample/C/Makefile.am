AM_CPPFLAGS = -I$(top_builddir)/include $(FFTW_INC) $(ESSL_INC) $(CUDA_INC) -DGIT_VERSION='"$(GIT_VERSION)"' -DGIT_DATE='"$(GIT_DATE)"'

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTW_LIB) $(FFTWF_LIB) $(ESSL_LIB) $(CUDA_LIB) -lstdc++

if ENABLE_CUDA
CCLD = nvcc
AM_CPPFLAGS += -DCUDA
AM_LDFLAGS = --std=c++11 -ccbin $(CXX_save)
else
CCLD = $(CXX)
endif

csampledir = $(datadir)/p3dfft++_samples/

csample_PROGRAMS = test1D_cos_c test1D_cos_complex_c test1D_r2c_c test3D_r2c_c test3D_c2c_c testIDIR_deriv_c test3D_c2c_inplace_c  test3D_r2c_memord_c test3D_r2c_single_c test2D+empty_c #testIDIR_deriv2_c

test3D_r2c_single_c_SOURCES = test3D_r2c_single.c

test1D_cos_c_SOURCES = test1D_cos.c

test1D_cos_complex_c_SOURCES = test1D_cos_complex.c

test1D_r2c_c_SOURCES = test1D_r2c.c

test3D_c2c_inplace_c_SOURCES = test3D_c2c_inplace.c

test3D_r2c_c_SOURCES = test3D_r2c.c

test3D_c2c_c_SOURCES = test3D_c2c.c

testIDIR_deriv_c_SOURCES = test_deriv.c

#testIDIR_deriv2_c_SOURCES = test_deriv2.c

test3D_r2c_memord_c_SOURCES = test3D_r2c_memord.c

test2D+empty_c_SOURCES = test2D+empty.c
