Installation
============

To install P3DFFT, run the following:

1. run './configure'
2. run 'make'
3. run 'make install'

Common Configurations
^^^^^^^^^^^^^^^^^^^^^
**Comet**

To build with Intel and FFTW::

  ./configure --enable-fftw

**Bridges**

To build with Intel and FFTW::

  ./configure --enable-fftw CXX=mpiicpc CC=mpiicc FC=mpiifort

**Stampede2**

To build with Intel and FFTW::

  ./configure --enable-fftw --with-fftw-lib=$TACC_FFTW3_LIB --with-fftw-inc=$TACC_FFTW3_INC

NOTE: Make sure you have the correct modules loaded (i.e. make sure you have FFTW module loaded if using FFTW)

Configuration Options
^^^^^^^^^^^^^^^^^^^^^

Run "configure --help" for a complete list of options

Compilers
^^^^^^^^^
--enable-intel  specifies to use intel compiler [default if none selected]
--enable-pgi    specifies to use pgi compiler
--enable-cray   specifies to use cray compiler
--enable-ibm    specifies to use ibm compiler
--enable-gnu    specifies to use gnu compiler

Fourier Transform Library
^^^^^^^^^^^^^^^^^^^^^^^^^
**FFTW**

--enable-fftw               prepares p3dfft++ to be used with fftw library
--with-fftwlib=liblocation  Specifically sets the lib directory location
--with-fftwinc=inclocation  Specifically sets the inc directory location
--enable-measure            For search-once-for-the-fast algorithm (takes more time on p3dfft_setup()). [default if none specified]
--enable-patient            For search-once-for-the-fastest-algorithm (takes much more time on p3dfft_setup()).
--enable-estimate           If this argument is passed, the FFTW library will not use run-time tuning to select the fastest algorithm for computing FFTs.

**ESSL**

--enable-essl               prepares p3dfft++ to be used with essl library
--with-essl=baselocation    tells configure where to look for the essl library (must be directory above lib and inc directories)

**MKL**

--enable-mkl                prepares p3dfft++ to be used with mkl library
--with-mkl=baselocation      tells configure where to look for the mkl library (must be directory above lib and inc directories)

