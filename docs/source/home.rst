Introduction
============
This site provides tools for solution of numerical problems in multiscale phenomena in three dimensions (3D). The most common example of such problem is Fast Fourier Transform (FFT), which is an important algorithm for simulations in a wide range of fields, including studies of turbulence, climatology, astrophysics and material science. Other algorithms of importance include Chebyshev transforms and high-order finite difference compact schemes.

Parallel Three-Dimensional Fast Fourier Transforms, dubbed P3DFFT, as well as its extension P3DFFT++, is a library for large-scale computer simulations on parallel platforms.This project was initiated at San Diego Supercomputer Center (SDSC) at UC San Diego by its main author Dmitry Pekurovsky, Ph.D.

This library uses 2D, or pencil, decomposition. This overcomes an important limitation to scalability inherent in FFT libraries implementing 1D (or slab) decomposition: the number of processors/tasks used to run this problem in parallel can be as large as N2, where N is the linear problem size. This approach has shown good scalability up to 524,288 cores.

**P3DFFT**

P3DFFT is written in Fortran90 and is optimized for parallel performance. It uses Message Passing Interface (MPI) for interprocessor communication, and starting from v.2.7.5 there is a multithreading option for hybrid MPI/OpenMP implementation. C/C++ interface is available, as are detailed documentation and examples in both Fortran and C. A configure script is supplied for ease of installation. This package depends on a serial FFT library such as Fastest Fourier Transform in the West (`FFTW <http://www.fftw.org/>`_) or IBM's `ESSL <http://publibfp.boulder.ibm.com/epubs/pdf/am501405.pdf>`_. The library is provided under GPL 3.0 and is available from its `github page <https://github.com/sdsc/p3dfft>`_.

**P3DFFT++**

P3DFFT++ is the next generation of P3DFFT (versions starting with 3.0). It extends the interface of P3DFFT to allow a wider range of use scenarios. It provides the user with a choice in defining their own data layout formats beyond the predefined 2D pencil blocks. It is written in C++ with C and Fortran interfaces, and currently uses MPI. The library can be found at P3DFFT++ `github space <https://github.com/sdsc/p3dfft.3>`_. See P3DFFT++ Tutorial and P3DFFT++ reference pages in C++, C and Fortran.

The following table compares P3DFFT family 2.7.6 and 3.0 (P3DFFT++).

=======================================   ============== ========
Feature                                   P3DFFT 2.x     P3DFFT++
=======================================   ============== ========
real-to-complex and complex-to-real FFT   Yes            Yes   

complex FFT                               No             Yes           

sine and cosine transforms                In 1 dimension Yes        

Chebyshev transform                       In 1 dimension Yes               

pruned transforms                         Yes            No               

In-place and out-of-place                 Yes            Yes             

Multiple grids                            No             Yes            

Hybrid MPI/OpenMP                         Yes            No            
=======================================   ============== ========

**License of use**

This software is provided for free for educational and not-for-profit use under a UCSD license. License terms can be seen here. Users are requested to complete optional registration when downloading this software, and also acknowledge the use as below.  

**Citation information**

Please acknowledge/cite use of P3DFFT as follows: D. Pekurovsky, P3DFFT: a framework for parallel computations of Fourier transforms in three dimensions, SIAM Journal on Scientific Computing 2012, Vol. 34, No. 4, pp. C192-C209. This paper can be obtained  `here <http://arxiv.org/abs/1905.02803>`_.

To cite the software you can also use DOI for P3DFFT v. 2.7.9:

.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2634590.svg
        :target: https://doi.org/10.5281/zenodo.2634590

.. figure:: https://www.nsf.gov/images/logos/nsf1.jpg
        :align: center

        This project is supported by National Science Foundation grant OAC-1835885.
