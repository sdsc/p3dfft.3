# P3DFFT++
This is a repository for P3DFFT++ (a.k.a. P3DFFT v. 3), an open source numerical library for spectral transforms (Fast Fourier Transforms and related algorithms). These algorithms are commonly used in simulating multiscale phenomena in various domains of science and engineering, in particular solving partial differential equations and inverse propagation problems. It takes the essence of P3DFFT library (see http://www.p3dfft.net) further by creating an extensible, modular structure uniquely adaptable to a greater range of use cases. The users can specify in detail what kind of data layout they would like to use, both in terms of local memory ordering and the processor layout. 

Just like P3DFFT, P3DFFT++ is a parallel software package, using MPI as the primary method for interprocessor communication. It supports 1D, 2D and 3D (to come soon) domain decomposition schemes. As P3DFFT, P3DFFT++ also relies on lower-level libraries, for example FFTW to perform optimized 1D FFTs.

Unlike P3DFFT, which was written in Fortran90, P3DFFT++ is written in C++. Interfaces are provided for C and Fortran. To learn about using the code the user is encouraged to study example programs in C++, C and FORTRAN subdirectories. Consult documentation at the main project website http://www.p3dfft.net. 

Please e-mail Dmitry Pekurovsky (dmitry@sdsc.edu) for any questions or suggestions. Software contributions are welcome, assuming they follow the main ideas of the framework.

This work is supported by the U.S. National Science Foundation
