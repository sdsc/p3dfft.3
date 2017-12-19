# p3dfft.3
This is a repository for P3DFFT++ (a.k.a. P3DFFT v. 3), a library for simulating multiscale phenomena. It takes the essence of P3DFFT library (see http://www.p3dfft.net) further by creating an extensible, modular structure uniquely adaptable to a greater range of use cases. The users can specify in detail what kind of data layout they would like to use, both in terms of local memory ordering and the processor layout. 

Just like P3DFFT, P3DFFT++ is a distributed software package, using MPI as the primary method for interprocessor commubnication. 
It supports 1D, 2D and 3D (to come soon) domain decomposition schemes. As P3DFFT, P3DFFT++ also relies on lower-level libraries, for example FFTW to perform optimized 1D FFTs.

Unlike P3DFFT, which was written in Fortran90, P3DFFT++ is written in C++. Interfaces are provided for C and Fortran. To learn about using the code the user is encouraged to study example programs in C++, C and FORTRAN subdirectories. Please e-mail Dmitry Pekurovsky (dmitry@sdsc.edu) for any questions or suggestions. Software contributions are welcome, assuming they follow the main ideas of the framework.
