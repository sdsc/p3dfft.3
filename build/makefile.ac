AM_CXXFLAGS = $(CPPFLAGS)

INCLUDES = $(FFTW_INC)

LDADD = $(FFTW_LIB)

ARFLAGS = -v -r -u

lib_LIBRARIES = libp3dfft.3.a

libp3dfft_3_a_SOURCES = init.o plan.o exec.o templ.o wrap.o

plan.o: templ.C

all-local:
	-@[ -e "p3dfft.mod" ] && mv -f p3dfft.mod ../include

clean-local:
	-[ -z "*.mod"] || rm -f *.mod
