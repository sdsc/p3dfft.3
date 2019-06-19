Common Errors
=============
I get an “undefined reference to MPI...” error when making the samples

  Make sure your compilers all link to the C++ MPI library. You can do this
  by adding the correct linking flag to your compiler flags as long as you are
  including the correct directory (example: FCFLAGS=”-lmpi_cxx”)
