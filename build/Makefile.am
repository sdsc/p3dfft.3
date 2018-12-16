AM_CXXFLAGS = $(CPPFLAGS)

INCLUDES = $(FFTW_INC)

ARFLAGS = -v -r -u

plan.o: templ.C

lib_LIBRARIES = libp3dfft.3.a

#LDADD = -L$(FFTW_LIB)

libp3dfft_3_a_SOURCES = init.C plan.C exec.C templ.C wrap.C

libp3dfft_3_a_LIBADD = -L$(FFTW_LIB) -lfftw3 -lfftw3f -lm

clean-local:
	-[ -z "*.mod" ] || rm -f *.mod
