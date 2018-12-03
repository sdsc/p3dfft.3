#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([P3DFFT++], [3.0.0], [dmitry@sdsc.edu])
AM_INIT_AUTOMAKE([-Wall])
AC_CONFIG_SRCDIR([build/init.C])
AC_CONFIG_HEADERS([config.h])

#TODO possibly make function for enabling certain compilers
AC_MSG_CHECKING([whether to use IBM compiler])
AC_ARG_ENABLE(ibm, [AC_HELP_STRING([--enable-ibm], [compile P3DFFT using IBM compiler])], enableibm=$enableval, enableibm=no)
AC_MSG_RESULT([$enableibm])
if test "$enableibm" = "yes"; then
	AC_DEFINE(IBM, 1, [Define if you want to compile P3DFFT using IBM compiler])
	compiler="ibm"
fi

AC_MSG_CHECKING([whether to use Intel compiler])
AC_ARG_ENABLE(intel, [AC_HELP_STRING([--enable-intel], [compile P3DFFT using Intel compiler])], enableintel=$enableval, enableintel=no)
AC_MSG_RESULT([$enableintel])
if test "$enableintel" = "yes"; then
	AC_DEFINE(INTEL, 1, [Define if you want to compile P3DFFT using Intel compiler])
	compiler="intel"
fi

AC_MSG_CHECKING([whether to use PGI compiler])
AC_ARG_ENABLE(pgi, [AC_HELP_STRING([--enable-pgi], [compile P3DFFT using PGI compiler])], enablepgi=$enableval, enablepgi=no)
AC_MSG_RESULT([$enablepgi])
if test "$enablepgi" = "yes"; then
	AC_DEFINE(PGI, 1, [Define if you want to compile P3DFFT using PGI compiler])
	compiler="pgi"
fi

AC_MSG_CHECKING([whether to use Cray compiler])
AC_ARG_ENABLE(cray, [AC_HELP_STRING([--enable-cray], [compile P3DFFT using IBM compiler])], enablecray=$enableval, enablecray=no)
AC_MSG_RESULT([$enablecray])
if test "$enablecray" = "yes"; then
	AC_DEFINE(Cray, 1, [Define if you want to compile P3DFFT using Cray compiler])
	compiler="cray"
fi

AC_MSG_CHECKING([whether to use GNU compiler])
AC_ARG_ENABLE(gnu, [AC_HELP_STRING([--enable-gnu], [compile P3DFFT using GNU compiler])], enablegnu=$enableval, enablegnu=no)
AC_MSG_RESULT([$enablegnu])
if test "$enablegnu" = "yes"; then
	AC_DEFINE(GNU, 1, [Define if you want to compile P3DFFT using GNU compiler])
	compiler="gnu"
fi

if test "$compiler" = ''; then
	AC_MSG_NOTICE([no compiler selected...defaulting to Intel])
	AC_DEFINE(INTEL, 1, [Define if you want to compile P3DFFT using Intel compiler])
	compiler="intel"
fi


AC_MSG_CHECKING([whether to use FFTW library])
AC_ARG_ENABLE(fftw, [AC_HELP_STRING([--enable-fftw], [for using the FFTW library])], fftwval=$enableval, fftwval=no)
if test "$fftwval" = "yes"; then
	AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
fi
AC_MSG_RESULT([$fftwval])

AC_MSG_CHECKING([whether to use ESSL library])
AC_ARG_ENABLE(essl, [AC_HELP_STRING([--enable-essl], [for using the ESSL library])], esslval=$enableval, esslval=no)
if test "$esslval" = "yes"; then
	AC_DEFINE(ESSL, 1, [Define if you want to use the ESSL library])
fi
AC_MSG_RESULT([$esslval])

# check to see if FFTW or ESSL has been selected
if test "$fftwval" = "no" && test "$esslval" = "no"; then
	AC_MSG_ERROR([You must choose to use either ESSL or FFTW!])
fi

# check to see if FFTW and ESSL have both been selected
if test "$fftwval" = "yes" && test "$esslval" = "yes"; then
	AC_MSG_ERROR([You must only choose either ESSL or FFTW!])
fi

AC_MSG_CHECKING([fftw location])
AC_ARG_WITH(fftw, [AC_HELP_STRING([--with-fftw], [path to FFTW directory])], withfftw=$withval, withfftw="")
if test "$withfftw" != "yes" && test "$withfftw" != ""; then
	AC_MSG_RESULT([$withfftw])
	if test -e [$withfftw/lib/libfftw3.a ]; then
		AC_SUBST(FFTW_INC, "-I$withfftw/include")
		AC_SUBST(FFTW_LIB, "$withfftw/lib/libfftw3.a")
	elif test -e [$withfftw/lib64/libfftw3.a ]; then
                AC_SUBST(FFTW_INC, "-I$withfftw/include")
                AC_SUBST(FFTW_LIB, "$withfftw/lib64/libfftw3.a")
	elif test -e [$withfftw/lib/x86_64-linux-gnu/libfftw3.a ]; then
                AC_SUBST(FFTW_INC, "-I$withfftw/include")
		AC_SUBST(FFTW_LIB, "$withfftw/lib/x86_64-linux-gnu/libfftw3.a")
	else
		AC_MSG_ERROR([libfftw3.a was not found in given location!])
	fi
else
	AC_MSG_RESULT([none])
fi

AC_MSG_CHECKING([fftw lib location])
AC_ARG_WITH(fftw-lib, [AC_HELP_STRING([--with-fftw-lib], [path to FFTW directory containing library file])], withfftwlib=$withval, withfftwlib="")
if test "$withfftwlib" != "" && test "$withfftw" != ""; then
	AC_MSG_ERROR([only choose --with-fftw or --with-fftw-lib!])
fi
if test "$withfftwlib" != "yes" && test "$withfftwlib" != "" && test "$fftwval" = "yes"; then
	AC_MSG_RESULT([$withfftwlib])
	if test -e [$withfftwlib/libfftw3.a ]
	then
		AC_SUBST(FFTW_LIB, "$withfftwlib/libfftw3.a")
	else
		AC_MSG_ERROR([libfftw3.a was not found in given location!])
	fi
else
	AC_MSG_RESULT([none])
fi

if test "$fftwval" = "yes" && test "$withfftwlib" = "" && test "$withfftw" = ""; then
	AC_MSG_ERROR([location to FFTW was not given!])
fi

AC_MSG_CHECKING([fftw inc location])
AC_ARG_WITH(fftw-inc, [AC_HELP_STRING([--with-fftw-inc], [path to FFTW include directory])], withfftwinc=$withval, withfftwinc="")
if (test "$withfftwinc" != "" && test "$withfftwlib" = "") || (test "$withfftwinc" = "" && test "$withfftwlib" != ""); then
	AC_MSG_ERROR([must use both --with-fftw-lib and --with-fftw-inc together!])
fi
if test "$withfftwinc" != "yes" && test "$withfftwinc" != "" && test "$fftwval" = "yes"; then
	AC_MSG_RESULT([$withfftwinc])
	AC_SUBST(FFTW_INC, "-I$withfftwinc")
else
	AC_MSG_RESULT([none])
fi

# Checks for programs.
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CPP
AM_PROG_AR

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lfftw3':
#AC_CHECK_LIB([fftw3], [main])
# FIXME: Replace `main' with a function in `-lfftw3f':
#AC_CHECK_LIB([fftw3f], [main])
# FIXME: Replace `main' with a function in `-lm':
#AC_CHECK_LIB([m], [main])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
		 build/Makefile])
		 #sample/Makefile
		 #sample/C++/Makefile
                 #sample/C/Makefile
                 #sample/FORTRAN/Makefile])
AC_OUTPUT
