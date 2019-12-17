#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([P3DFFT++], [3.0.0], [dmitry@sdsc.edu])
AM_INIT_AUTOMAKE([-Wall no-dependencies])
AC_CONFIG_SRCDIR([build/init.C])
AC_CONFIG_HEADERS([p3dfft3config.h])
AM_MAINTAINER_MODE([disable])


#: ${CFLAGS="-g -O0"}

#CXXFLAGS="${CXXFLAGS} -std=c++11"

#: ${FCFLAGS="-g -O0"}


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
	if test "$fourlib" = ""; then fourlib="fftw"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
  AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
fi
AC_MSG_RESULT([$fftwval])

AC_MSG_CHECKING([whether to use ESSL library])
AC_ARG_ENABLE(essl, [AC_HELP_STRING([--enable-essl], [for using the ESSL library])], esslval=$enableval, esslval=no)
if test "$esslval" = "yes"; then
	if test "$fourlib" = ""; then fourlib="essl"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
  AC_DEFINE(ESSL, 1, [Define if you want to use the ESSL library])
fi
AC_MSG_RESULT([$esslval])

AC_MSG_CHECKING([whether to use MKL library])
AC_ARG_ENABLE(mkl, [AC_HELP_STRING([--enable-mkl], [for using the MKL library])], mklval=$enableval, mklval=no)
if test "$mklval" = "yes"; then
	if test "$fourlib" = ""; then fourlib="mkl"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(MKL_BLAS, 1, [Define if you want to use the MKL library])
fi
AC_MSG_RESULT([$mklval])

AS_IF([test "$fourlib" = "fftw"],
      [AC_MSG_CHECKING([FFTW directory])
       AC_ARG_WITH(fftw, [AC_HELP_STRING([--with-fftw], [path to location of FFTW directory])], withfftw=$withval, withfftw="")
       AS_IF([test "$withfftw" != ""],
             [AC_MSG_RESULT([$withfftw])

              LIBS_save=$LIBS

              AS_IF([test -d "${withfftw}/lib"],
                    [LDFLAGS_save=$LDFLAGS
                     LDFLAGS="-L${withfftw}/lib"
                     AC_SEARCH_LIBS([fftw_plan_dft], [fftw3],
                                    [AC_SUBST(FFTW_LIB, "-L${withfftw}/lib -lfftw3")],
                                    [AC_MSG_ERROR([libfftw3 not found!])])
                     AC_SEARCH_LIBS([fftwf_execute_r2r], [fftw3f],
                                    [AC_SUBST(FFTWF_LIB, "-lfftw3f")],
                                    [AC_MSG_ERROR([libfftw3f not found!])])
                    LDFLAGS=$LDFLAGS_save],
                    [AC_MSG_ERROR([location to FFTW lib directory is not found!])])

              AS_IF([test -d "${withfftw}/include"],
                    [AC_CHECK_HEADERS([${withfftw}/include/fftw3.h],
                                      [AC_SUBST(FFTW_INC, "-I${withfftw}/include")],
                                      [AC_MSG_ERROR([fftw3.h not found in "${withfftw}/include"])])],
                    [AC_MSG_ERROR([location to FFTW include directory is invalid!])])

              LIBS=$LIBS_save

              fftwflag=
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

              AC_MSG_RESULT([$fftwflag])],
             [AC_MSG_RESULT([none])
              if test "$fftwval" = "yes"; then AC_MSG_ERROR([location to FFTW library directory was not given!])
              fi])])

AC_MSG_CHECKING([ESSL directory])
AC_ARG_WITH(essl, [AC_HELP_STRING([--with-essl], [path to location of ESSL library])], withessl=$withval, withessl="")
if test "$withessl" != "yes" && test "$withessl" != "" && test "$withessl" != "no"; then
	if test "$fourlib" != "essl"; then AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi

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
