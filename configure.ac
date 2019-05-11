#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([P3DFFT++], [3.0.0], [dmitry@sdsc.edu])
AM_INIT_AUTOMAKE([-Wall no-dependencies])
AC_CONFIG_SRCDIR([build/init.C])
AC_CONFIG_HEADERS([p3dfft3config.h])

#: ${CFLAGS="-g -O0"}

#: ${CXXFLAGS="-g -O0"}

#: ${FCFLAGS="-g -O0"}

DEFINE="#define"
UNDEFINE="#undef"

git_date=`git --no-pager show --date=short --format="%ad" --name-only 2> /dev/null | { read first _ ; echo $first ; } `
git_version=`git --no-pager describe --tags --always --dirty 2> /dev/null`
AC_SUBST(GIT_DATE, "$git_date")
AC_SUBST(GIT_VERSION, "$git_version")


AC_MSG_CHECKING([whether to use IBM compiler])
AC_ARG_ENABLE(ibm, [AC_HELP_STRING([--enable-ibm], [compile P3DFFT using IBM compiler])], enableibm=$enableval, enableibm=no)
AC_MSG_RESULT([$enableibm])
if test "$enableibm" = "yes"; then
	AC_DEFINE(IBM, 1, [Define if you want to compile P3DFFT using IBM compiler])
	if test "$compiler" == ''; then compiler="ibm"
	else AC_MSG_ERROR([Only choose one compiler!])
	fi
fi

AC_MSG_CHECKING([whether to use Intel compiler])
AC_ARG_ENABLE(intel, [AC_HELP_STRING([--enable-intel], [compile P3DFFT using Intel compiler])], enableintel=$enableval, enableintel=no)
AC_MSG_RESULT([$enableintel])
if test "$enableintel" = "yes"; then
	AC_DEFINE(INTEL, 1, [Define if you want to compile P3DFFT using Intel compiler])
	if test "$compiler" == ''; then compiler="intel"
	else AC_MSG_ERROR([Only choose one compiler!])
	fi
fi

AC_MSG_CHECKING([whether to use PGI compiler])
AC_ARG_ENABLE(pgi, [AC_HELP_STRING([--enable-pgi], [compile P3DFFT using PGI compiler])], enablepgi=$enableval, enablepgi=no)
AC_MSG_RESULT([$enablepgi])
if test "$enablepgi" = "yes"; then
	AC_DEFINE(PGI, 1, [Define if you want to compile P3DFFT using PGI compiler])
	if test "$compiler" == ''; then compiler="pgi"
	else AC_MSG_ERROR([Only choose one compiler!])
	fi
fi

AC_MSG_CHECKING([whether to use Cray compiler])
AC_ARG_ENABLE(cray, [AC_HELP_STRING([--enable-cray], [compile P3DFFT using IBM compiler])], enablecray=$enableval, enablecray=no)
AC_MSG_RESULT([$enablecray])
if test "$enablecray" = "yes"; then
	AC_DEFINE(Cray, 1, [Define if you want to compile P3DFFT using Cray compiler])
	if test "$compiler" == ''; then compiler="cray"
	else AC_MSG_ERROR([Only choose one compiler!])
	fi
fi

AC_MSG_CHECKING([whether to use GNU compiler])
AC_ARG_ENABLE(gnu, [AC_HELP_STRING([--enable-gnu], [compile P3DFFT using GNU compiler])], enablegnu=$enableval, enablegnu=no)
AC_MSG_RESULT([$enablegnu])
if test "$enablegnu" = "yes"; then
	AC_DEFINE(GNU, 1, [Define if you want to compile P3DFFT using GNU compiler])
	if test "$compiler" == ''; then compiler="gnu"
	else AC_MSG_ERROR([Only choose one compiler!])
	fi
fi

if test "$compiler" = ''; then
	AC_MSG_NOTICE([no compiler selected...defaulting to Intel])
	AC_DEFINE(INTEL, 1, [Define if you want to compile P3DFFT using Intel compiler])
	compiler="intel"
fi

if test "$CC" = ''; then
	CC="mpicc"
fi
if test "$CXX" = ''; then
	CXX="mpicxx"
fi
if test "$FC" = ''; then
	FC="mpif90"
fi


AC_MSG_CHECKING([whether to use TIMERS])
AC_ARG_ENABLE(timers, [AC_HELP_STRING([--enable-timers], [for using timers (DEBUGGING PURPOSES ONLY)])], timersval=$enableval, timersval=no)
if test "$timersval" = "yes"; then AC_DEFINE(TIMERS, 1, [Define if you want to use timers])
fi
AC_MSG_RESULT([$timersval])

fourlib=""
AC_MSG_CHECKING([whether to use FFTW library])
AC_ARG_ENABLE(fftw, [AC_HELP_STRING([--enable-fftw], [for using the FFTW library])], fftwval=$enableval, fftwval=no)
if test "$fftwval" = "yes"; then
	if "$fourlib" = ""; then fourlib="fftw"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
fi
AC_MSG_RESULT([$fftwval])

AC_MSG_CHECKING([whether to use ESSL library])
AC_ARG_ENABLE(essl, [AC_HELP_STRING([--enable-essl], [for using the ESSL library])], esslval=$enableval, esslval=no)
if test "$esslval" = "yes"; then
	if "$fourlib" = ""; then fourlib="essl"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
fi
AC_MSG_RESULT([$esslval])

AC_MSG_CHECKING([whether to use MKL library])
AC_ARG_ENABLE(mkl, [AC_HELP_STRING([--enable-mkl], [for using the MKL library])], mklval=$enableval, mklval=no)
if test "$mklval" = "yes"; then
	if "$fourlib" = ""; then fourlib="mkl"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
fi
AC_MSG_RESULT([$mklval])

##########DEPRECIATED###########################################################
AC_MSG_CHECKING([fftw location])
AC_ARG_WITH(fftw, [AC_HELP_STRING([--with-fftw], [path to FFTW directory])], withfftw=$withval, withfftw="")
if test "$withfftw" != "yes" && test "$withfftw" != "" && test "$withfftw" != "no"; then
	if test "$fourlib" = "fftw" || test "$fourlib" = ""; then fourlib="fftw"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
	AC_MSG_RESULT([$withfftw])
	if test -e [$withfftw/lib/libfftw3.a ] && test -e [$withfftw/lib/libfftw3f.a ]; then
 		AC_SUBST(FFTW_INC, "-I$withfftw/include")
 		AC_SUBST(FFTW_LIB, "$withfftw/lib/libfftw3.a")
 		AC_SUBST(FFTWF_LIB, "$withfftw/lib/libfftw3f.a")
 	elif test -e [$withfftw/lib64/libfftw3.a ] && test -e [$withfftw/lib64/libfftw3f.a ]; then
    AC_SUBST(FFTW_INC, "-I$withfftw/include")
    AC_SUBST(FFTW_LIB, "$withfftw/lib64/libfftw3.a")
    AC_SUBST(FFTWF_LIB, "$withfftw/lib64/libfftw3f.a")
 	elif test -e [$withfftw/lib/x86_64-linux-gnu/libfftw3.a ] && test -e [$withfftw/lib/x86_64-linux-gnu/libfftw3f.a ]; then
    AC_SUBST(FFTW_INC, "-I$withfftw/include")
 		AC_SUBST(FFTW_LIB, "$withfftw/lib/x86_64-linux-gnu/libfftw3.a")
 		AC_SUBST(FFTWF_LIB, "$withfftw/lib/x86_64-linux-gnu/libfftw3f.a")
 	else
 		AC_MSG_ERROR([libfftw3.a was not found in given location!])
 	fi
else
	AC_MSG_RESULT([none])
fi
################################################################################

AC_MSG_CHECKING([fftw lib location])
AC_ARG_WITH(fftw-lib, [AC_HELP_STRING([--with-fftw-lib], [path to FFTW directory containing library file])], withfftwlib=$withval, withfftwlib="")
if !( test "$withfftwlib" = "" || test "$withfftwlib" = "no" ) && !( test "$withfftw" = "" || test "$withfftw" = "no" ); then
	AC_MSG_ERROR([only choose --with-fftw or --with-fftw-lib!])
fi
if test "$withfftwlib" != "yes" && test "$withfftwlib" != ""; then
	if test "$fourlib" = "fftw" || test "$fourlib" = ""; then fourlib="fftw"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
	AC_MSG_RESULT([$withfftwlib])
	if test -e [$withfftwlib/libfftw3.a ] && test -e [$withfftwlib/libfftw3f.a ]; then
		AC_SUBST(FFTW_LIB, "$withfftwlib/libfftw3.a")
		AC_SUBST(FFTWF_LIB, "$withfftwlib/libfftw3f.a")
	else
		AC_MSG_ERROR([libfftw3.a was not found in given location!])
	fi
else
	AC_MSG_RESULT([none])
fi


AC_MSG_CHECKING([fftw inc location])
AC_ARG_WITH(fftw-inc, [AC_HELP_STRING([--with-fftw-inc], [path to FFTW include directory])], withfftwinc=$withval, withfftwinc="")
if (test "$withfftwinc" != "" && test "$withfftwlib" = "") || (test "$withfftwinc" = "" && test "$withfftwlib" != ""); then
	AC_MSG_ERROR([must use both --with-fftw-lib and --with-fftw-inc together!])
fi
if test "$withfftwinc" != "yes" && test "$withfftwinc" != ""; then
	AC_MSG_RESULT([$withfftwinc])
	AC_SUBST(FFTW_INC, "-I$withfftwinc")
else
	AC_MSG_RESULT([none])
fi


AC_MSG_CHECKING([ESSL directory])
AC_ARG_WITH(essl, [AC_HELP_STRING([--with-essl], [path to location of ESSL library])], withessl=$withval, withessl="")
if test "$withessl" != "yes" && test "$withessl" != "" && test "$withessl" != "no"; then
	if test "$fourlib" != "essl"; then AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(ESSL, 1, [Define if you want to use the ESSL library])
	AC_MSG_RESULT([$withessl])
	if [[ -d "$withessl/lib64" ]]; then
		AC_SUBST(ESSL_INC, "-I$withessl/include")
    AC_SUBST(ESSL_LIB, "-L$withessl/lib64 -lesslbg -L$IBM_MAIN_DIR/xlf_r/bg/14.1/lib64")
  else
    AC_MSG_ERROR([essl library was not found in given location!])
	fi
else
  AC_MSG_RESULT([none])
  if test "$esslval" = "yes"; then AC_MSG_ERROR([location to ESSL was not given!])
	fi
fi


AC_MSG_CHECKING([MKL directory])
AC_ARG_WITH(mkl, [AC_HELP_STRING([--with-mkl], [path to location of MKL library])], withmkl=$withval, withmkl="")
if test "$withmkl" != "yes" && test "$withmkl" != "" && test "$withmkl" != "no"; then
	if test "$fourlib" != "mkl"; then AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(MKL_BLAS, 1, [Define if you want to use the MKL library])
	AC_MSG_RESULT([$withmkl])
	if test -d $withmkl; then
		AC_SUBST(MKL_INC, "-I$withmkl/include")
		AC_SUBST(MKL_LIB, ["-Wl,--start-group $withmkl/lib/intel64/libmkl_intel_ilp64.a $withmkl/lib/intel64/libmkl_sequential.a $withmkl/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl"])
	else
		AC_MSG_ERROR([mkl library was not found in given location!])
	fi
else
	AC_MSG_RESULT([none])
	if test "$mklval" = "yes"; then AC_MSG_ERROR([location to MKL was not given!])
	fi
fi


# check to see if essl or ESSL or MKL has been selected
if test "$fourlib" = ""; then
	AC_MSG_ERROR([You must choose to use either ESSL or FFTW or MKL!])
fi


#default value for fftwflag
fftwflag=""
if test "$fourlib" = "fftw"; then
	AC_MSG_CHECKING([fftw flag])
	AC_ARG_ENABLE(fftwmeasure, [AC_HELP_STRING([--enable-fftwmeasure], [measure planner flag for FFTW library])], fftwmeasureval=$enableval, fftwmeasureval="")
	if test "$fftwmeasureval" = "yes"; then
		if test "$fftwflag" = ""; then fftwflag="measure"
		else AC_MSG_ERROR([Only choose one FFTW planner flag!])
		fi
		AC_DEFINE(FFTW_FLAG_MEASURE, 1, [Define if you want to use the measure fftw planner flag])
	fi

	AC_ARG_ENABLE(fftwestimate, [AC_HELP_STRING([--enable-fftwestimate], [estimate planner flag for FFTW library])], fftwestimateval=$enableval, fftwestimateval="")
	if test "$fftwestimateval" = "yes"; then
		if test "$fftwflag" = ""; then fftwflag="estimate"
		else AC_MSG_ERROR([Only choose one FFTW planner flag!])
		fi
		AC_DEFINE(FFTW_FLAG_ESTIMATE, 1, [Define if you want to use the estimate fftw planner flag])
	fi

	AC_ARG_ENABLE(fftwpatient, [AC_HELP_STRING([--enable-fftwpatient], [patient planner flag for FFTW library])], fftwpatientval=$enableval, fftwpatientval="")
	if test "$fftwpatientval" = "yes"; then
		if test "$fftwflag" = ""; then fftwflag="patient"
		else AC_MSG_ERROR([Only choose one FFTW planner flag!])
		fi
		AC_DEFINE(FFTW_FLAG_PATIENT, 1, [Define if you want to use the patient fftw planner flag])
	fi

	if test "$fftwflag" = ""; then
		fftwflag="measure"
		AC_DEFINE(FFTW_FLAG_MEASURE, 1, [Define if you want to use the measure fftw planner flag])
	fi
	AC_MSG_RESULT([$fftwflag])
fi


# Checks for programs.
AC_PROG_CXX
AC_PROG_FC
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_CPP
AM_PROG_AR
AC_PROG_RANLIB

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
		 build/Makefile
		 sample/Makefile
		 sample/C++/Makefile
                 sample/C/Makefile
                 sample/FORTRAN/Makefile])
AC_OUTPUT
mv p3dfft3config.h include/p3dfft3config.h
AC_MSG_NOTICE([configure script complete!])
