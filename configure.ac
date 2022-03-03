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
	CC="mpicxx"
fi
if test "$CXX" = ''; then
	CXX="mpicxx"
fi
if test "$FC" = ''; then
	FC="mpif90"
fi

AC_SUBST(CXX_save, "$CXX_save")
AC_SUBST(FC_save, "$FC_save")
fourlib=""

AC_MSG_CHECKING([whether to use Non-Blocking version])
AC_ARG_ENABLE(nonblocking, [AC_HELP_STRING([--enable-nonblocking], [for using nonblocking version])], nbval=$enableval, nbval=no)
if test "$nbval" = "yes"; then AC_DEFINE(NB, 1, [Define if you want to use nonblocking version])
fi
AC_MSG_RESULT([$nbval])

AS_IF([test "$nbval" = "yes"],
	    withnb="a2a"
      [AC_MSG_CHECKING([Nonblocking method])
       AC_ARG_WITH(nonblocking, [AC_HELP_STRING([--with-nonblocking], [nonblocking method])], withnb=$withval, withnb="")]
       AS_IF([test "$withnb" != ""],				      
             [AC_MSG_RESULT([$withnb])] 
	     if test "$withnb" = "p2p"; then AC_DEFINE(P2P, 1, [Point-to-Point NB communication])
	     fi
	     if test "$withnb" = "a2a"; then AC_DEFINE(A2A, 1, [All-to-all NB communication])
	     fi
	     )
	     )
	     

AC_MSG_CHECKING([whether to use TIMERS])
AC_ARG_ENABLE(timers, [AC_HELP_STRING([--enable-timers], [for using timers (DEBUGGING PURPOSES ONLY)])], timersval=$enableval, timersval=no)
if test "$timersval" = "yes"; then AC_DEFINE(TIMERS, 1, [Define if you want to use timers])
fi
AC_MSG_RESULT([$timersval])

AC_MSG_CHECKING([whether to use FFTW library])
AC_ARG_ENABLE(fftw, [AC_HELP_STRING([--enable-fftw], [for using the FFTW library])], fftwval=$enableval, fftwval=no)
if test "$fftwval" = "yes"; then
	if test "$fourlib" = ""; then fourlib="fftw"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
  AC_DEFINE(FFTW, 1, [Define if you want to use the FFTW library])
fi
AC_MSG_RESULT([$fftwval])

AC_MSG_CHECKING([whether to use CUDA library])
AC_ARG_ENABLE(cuda, [AC_HELP_STRING([--enable-cuda], [for using the CUDA library])], cudaval=$enableval, cudaval=no)
if test "$cudaval" = "yes"; then
	if test "$fourlib" = ""; then fourlib="cuda"
	else AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi
	AC_DEFINE(CUDA, 1, [Define if you want to use the CUDA library])
	AM_CONDITIONAL([ENABLE_CUDA], [true])
else
	AM_CONDITIONAL([ENABLE_CUDA], [false])
fi
AC_MSG_RESULT([$cudaval])

AC_MSG_CHECKING([whether to use cuBLAS library])
AC_ARG_ENABLE(cublas, [AC_HELP_STRING([--enable-cublas], [for using the cuBLAS library])], cublasval=$enableval, cublasval=no)
if test "$cublasval" = "yes"; then
	if test "$fourlib" = "cuda"; then
		AC_DEFINE(CUBLAS, 1, [Define if you want to use the cuBLAS library])
		AM_CONDITIONAL([ENABLE_CUBLAS], [true])
	else
		AC_MSG_ERROR([CUDA must be enabled in order to use cuBLAS!])
	fi
else
	AM_CONDITIONAL([ENABLE_CUBLAS], [false])
fi
AC_MSG_RESULT([$cublasval])

AC_MSG_CHECKING([whether to use cuTensor library])
AC_ARG_ENABLE(cutensor, [AC_HELP_STRING([--enable-cutensor], [for using the cuTensor library])], cutensorval=$enableval, cutensorval=no)
if test "$cutensorval" = "yes"; then
	if test "$fourlib" = "cuda"; then
		AC_DEFINE(CUTENSOR, 1, [Define if you want to use the cuTensor library])
		AM_CONDITIONAL([ENABLE_CUTENSOR], [true])
	else
		AC_MSG_ERROR([CUDA must be enabled in order to use cuTensor!])
	fi
else
	AM_CONDITIONAL([ENABLE_CUTENSOR], [false])
fi
AC_MSG_RESULT([$cutensorval])

AC_MSG_CHECKING([whether to use cuTT library])
AC_ARG_ENABLE(cutt, [AC_HELP_STRING([--enable-cutt], [for using the cuTT library])], cuttval=$enableval, cuttval=no)
if test "$cuttval" = "yes"; then
	if test "$fourlib" = "cuda"; then
		AC_DEFINE(CUTT, 1, [Define if you want to use the cuTT library])
		AM_CONDITIONAL([ENABLE_CUTT], [true])
	else
		AC_MSG_ERROR([CUDA must be enabled in order to use cuTT!])
	fi
else
	AM_CONDITIONAL([ENABLE_CUTT], [false])
fi
AC_MSG_RESULT([$cuttval])

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

              LIBS=$LIBS_save

              AS_IF([test -d "${withfftw}/include"],
                    [AC_CHECK_HEADERS([${withfftw}/include/fftw3.h],
                                      [AC_SUBST(FFTW_INC, "-I${withfftw}/include")],
                                      [AC_MSG_ERROR([fftw3.h not found in "${withfftw}/include"])])],
                    [AC_MSG_ERROR([location to FFTW include directory is invalid!])])

              CXX_save="$CXX"
	      FC_save="$FC"

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

AC_MSG_CHECKING([CUDA directory])
AC_ARG_WITH(cuda, [AC_HELP_STRING([--with-cuda], [path to location of CUDA library])], withcuda=$withval, withcuda="")
if test "withcuda" != "yes" && test "$withcuda" != "" && test "$withcuda" != "no"; then
	if test "$fourlib" != "cuda"; then AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi

	AC_MSG_RESULT([$withcuda])
	if [[ -d "$withcuda/samples/common/inc" ]]; then
		CUDA_LIB="-L$withcuda/lib64 -lcuda -lcufft"
		CUDA_INC="-I$withcuda/include -I$withcuda/samples/common/inc"
	else
		AC_MSG_ERROR([cuda library was not found in given location!])
	fi

	if test "$cublasval" = "yes"; then
		CUDA_LIB+=" -lcublas"
	fi

	if test "$cutensorval" = "yes"; then
		CUDA_LIB_save=$CUDA_LIB

		AC_ARG_WITH(cutensor, [AC_HELP_STRING([--with-cutensor], [path to location of cuTensor library])], withcutensor=$withval, withcutensor="")
		if test "withcutensor" != "yes" && test "$withcutensor" != "" && test "$withcutensor" != "no"; then
			CUDA_INC+=" -I$withcutensor/include"
#			CUDA_LIB=" -L$withcutensor/lib"
		fi
		AC_ARG_WITH(cutensor-lib, [AC_HELP_STRING([--with-cutensor-lib], [path to location of cuTensor lib directory])], withcutensorlib=$withval, withcutensorlib="")
		if test "$withcutensorlib" != ""; then
			CUDA_LIB=" -L$withcutensorlib"
		fi

		CUDA_LIB_save+=$CUDA_LIB
		CUDA_LIB=$CUDA_LIB_save
		CUDA_LIB+=" -lcutensor_static -lcudart"
	fi

	if test "$cuttval" = "yes"; then
		AC_ARG_WITH(cutt, [AC_HELP_STRING([--with-cutt], [path to location of cuTT library])], withcutt=$withval, withcutt="")
		if test "withcutt" != "yes" && test "$withcutt" != "" && test "$withcutt" != "no"; then
			CUDA_INC+=" -I$withcutt/include"
			CUDA_LIB+=" -L$withcutt/lib -lcutt"
		fi
	fi

	AC_SUBST(CUDA_LIB, "$CUDA_LIB")
	AC_SUBST(CUDA_INC, "$CUDA_INC")
	CXX_save="$CXX"
	FC_save="$FC"
else
	AC_MSG_RESULT([none])
	if test "$cudaval" = "yes"; then AC_MSG_ERROR([location to CUDA was not given!])
	fi
fi

AC_MSG_CHECKING([ESSL directory])
AC_ARG_WITH(essl, [AC_HELP_STRING([--with-essl], [path to location of ESSL library])], withessl=$withval, withessl="")
if test "$withessl" != "yes" && test "$withessl" != "" && test "$withessl" != "no"; then
	if test "$fourlib" != "essl"; then AC_MSG_ERROR([Only choose one Fourier transform library!])
	fi

	AC_MSG_RESULT([$withessl])
	if [[ -d "$withessl/lib64" ]]; then
		AC_SUBST(ESSL_INC, "-I$withessl/include")
    AC_SUBST(ESSL_LIB, "-L$withessl/lib64 -lessl") #-L$IBM_MAIN_DIR/xlf_r/bg/14.1/lib64")
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
	AC_MSG_ERROR([You must choose to use either ESSL or FFTW or MKL or CUDA!])
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

