/* p3dfft3config.h.  Generated from p3dfft3config.h.in by configure.  */
/* p3dfft3config.h.in.  Generated from configure.ac by autoheader.  */

/* Define if you want to compile P3DFFT using Cray compiler */
/* #undef Cray */

/* Define if you want to use the ESSL library */
/* #undef ESSL */

/* Define if you want to use the FFTW library */
#define FFTW 1

/* Define if you want to use the estimate fftw planner flag */
/* #undef FFTW_FLAG_ESTIMATE */

/* Define if you want to use the measure fftw planner flag */
#define FFTW_FLAG_MEASURE 1

/* Define if you want to use the patient fftw planner flag */
/* #undef FFTW_FLAG_PATIENT */

/* Define if you want to compile P3DFFT using GNU compiler */
/* #undef GNU */

/* Define to 1 if you have the <fftw3.h> header file. */
/* #undef HAVE_FFTW3_H */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if the system has the type `_Bool'. */
#define HAVE__BOOL 1

/* Define if you want to compile P3DFFT using IBM compiler */
/* #undef IBM */

/* Define if you want to compile P3DFFT using Intel compiler */
#define INTEL 1

/* Define if you want to use the MKL library */
/* #undef MKL_BLAS */

/* Name of package */
#define PACKAGE "p3dfft--"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "dmitry@sdsc.edu"

/* Define to the full name of this package. */
#define PACKAGE_NAME "P3DFFT++"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "P3DFFT++ 3.0.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "p3dfft--"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "3.0.0"

/* Define if you want to compile P3DFFT using PGI compiler */
/* #undef PGI */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Define if you want to use timers */
/* #undef TIMERS */

/* Version number of package */
#define VERSION "3.0.0"

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
