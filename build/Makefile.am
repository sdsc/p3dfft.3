AM_CXXFLAGS = $(FFTW_INC) -I$(top_builddir)/include $(ESSL_INC) $(MKL_INC) $(MKL_LIB)

ARFLAGS = -v -r -u

plan.o: $(top_builddir)/include/p3dfft.h templ.C

init.o: $(top_builddir)/include/p3dfft.h

exec.o: $(top_builddir)/include/p3dfft.h

templ.o: $(top_builddir)/include/p3dfft.h

wrap.o: $(top_builddir)/include/p3dfft.h

deriv.o: $(top_builddir)/include/p3dfft.h

lib_LIBRARIES = libp3dfft.3.a

LDADD = $(FFTWF_LIB) $(FFTW_LIB) $(ESSL_LIB) -lm

libp3dfft_3_a_SOURCES = init.C plan.C exec.C templ.C wrap.C deriv.C

clean-local:
	-[ -z "*.mod" ] || rm -f *.mod
