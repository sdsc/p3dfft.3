AM_CXXFLAGS = -I$(top_builddir)/include $(FFTW_INC) $(ESSL_INC) $(CUDA_INC) --std=c++11

ARFLAGS = -v -r -u

if ENABLE_CUDA
DEFS += -ccbin $(CXX_save) -x cu --compile -DCUDA
LINKER = nvcc
if ENABLE_CUBLAS
DEFS += -DCUBLAS
endif
if ENABLE_CUTENSOR
DEFS += -DCUTENSOR
endif
if ENABLE_CUTT
DEFS += -DCUTT
endif
else
LINKER = $(CXX_save)
endif

CXX = $(LINKER)


#plan.o: $(top_builddir)/include/p3dfft.h templ.C

init.o: $(top_builddir)/include/p3dfft.h

exec.o: $(top_builddir)/include/p3dfft.h

templ.o: $(top_builddir)/include/p3dfft.h

#wrap.o: $(top_builddir)/include/p3dfft.h

#deriv.o: $(top_builddir)/include/p3dfft.h

lib_LIBRARIES = libp3dfft.3.a

LDADD = $(FFTW_LIB) $(FFTWF_LIB) $(ESSL_LIB) $(CUDA_LIB) -lm -lstdc++

libp3dfft_3_a_SOURCES = init.C exec.C templ.C #wrap.C deriv.C

clean-local:
	-[ -z "*.mod" ] || rm -f *.mod

sources = $(libp3dfft_3_a_SOURCES)

.C.d:
	@set -e;
	rm -f $@;
	$(CXXCOMPILE) -MM $< > $@.temporary;
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.temporary > $@;
	rm -f $@.temporary;

-include $(sources:.C=.d)
