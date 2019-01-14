AM_CXXFLAGS = $(FFTW_INC) -I$(top_builddir)/include -include $(top_builddir)/include/config.h 

ARFLAGS = -v -r -u

plan.o: templ.C

lib_LIBRARIES = libp3dfft.3.a

LDADD = $(FFTWF_LIB) $(FFTW_LIB) -lm

libp3dfft_3_a_SOURCES = init.C plan.C exec.C templ.C wrap.C

clean-local:
	-[ -z "*.mod" ] || rm -f *.mod
