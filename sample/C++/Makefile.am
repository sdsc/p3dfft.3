AM_CXXFLAGS = -I$(top_builddir)/include $(FFTW_INC) $(ESSL_INC) $(CUDA_INC) -DGIT_VERSION='"$(GIT_VERSION)"' -DGIT_DATE='"$(GIT_DATE)"'

LDADD = -lm $(top_builddir)/build/libp3dfft.3.a $(FFTW_LIB) $(FFTWF_LIB) $(ESSL_LIB) $(CUDA_LIB) -lstdc++

if ENABLE_CUDA
CXX_LINKER = nvcc
AM_CXXFLAGS += -DCUDA
AM_LDFLAGS = --std=c++11 -ccbin $(CXX_save)
else
CXX_LINKER = $(CXX)
endif

CXXLD = $(CXX_LINKER)

cppsampledir = $(datadir)/p3dfft++_samples/

cppsample_PROGRAMS = test3D_c2c_cpp test3D_r2c_cpp test1D_transplan_cpp  test3D_c2c_inplace_cpp testIDIR_deriv_cpp  test3D_r2c_single_cpp #testIDIR_deriv2_cpp #test1D_cos_cpp test1D_cos_complex_cpp test1D_sin_cpp

testIDIR_deriv_cpp_SOURCES = test_deriv.C

testIDIR_deriv2_cpp_SOURCES = test_deriv2.C

test3D_c2c_inplace_cpp_SOURCES = test3D_c2c_inplace.C

test3D_c2c_cpp_SOURCES = test3D_c2c.C

test3D_r2c_cpp_SOURCES = test3D_r2c.C

test3D_r2c_single_cpp_SOURCES = test3D_r2c_single.C

test1D_transplan_cpp_SOURCES = test_transplan.C

test1D_cos_cpp_SOURCES = test1D_cos.C

test1D_cos_complex_cpp_SOURCES = test1D_cos_complex.C

test1D_sin_cpp_SOURCES = test1D_sin.C
