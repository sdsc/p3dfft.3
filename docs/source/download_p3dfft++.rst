Download P3DFFT++
=================
The latest release of P3DFFT++ can be downloaded `here <https://github.com/sdsc/p3dfft.3>`_.

To install, modify the *makefile* to set C++ MPI-enabled compiler appropriate to your system. Also edit the location of FFTW or MKL libraries. Several site-specific makefiles are provided as examples. *make lib* will build the library.

In addition, if you wish to build example programs in C++, C and/or Fortran, follow the same steps for *makefile* in each of the subdirectories. Then type *make samples* in the top directory to build all 3 language examples. *make all* or *make builds* both library and examples. 

.. raw:: html

        <iframe src="https://docs.google.com/forms/d/e/1FAIpQLScS7KDRdACJJNwZHi8khlVX1SYNAd61AjlPM23f_7bq8yVN2Q/viewform?embedded=true" width="640" height="1222" frameborder="0" marginheight="0" marginwidth="0">Loading…</iframe>
