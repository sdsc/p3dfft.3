P3DFFT++ Fortran examples
*************************

.. csv-table::
        :widths: auto

        "test3D_r2c.f90", "Test program for real-to-complex and complex-to-real 3D FFT. Can be used as a correctness test and a timing test."
        "test3D_c2c.f90", "Test program for complex-to-complex forward and backward 3D FFT. Can be used as a correctness test and a timing test."
        "test_trans1D.f90", "Test program for real-to-complex and complex-to-real 1D FFT. User can choose an arbitrary dimension of transform. This function performs a local transpose in addition to 1D FFT, as needed. The transpose is specified by input and output memory order, which maps logical onto storage array dimensions. No inter-processor transpose takes place here. 1D transforms can be used by themselves or as a stage in a multidimensional transform."
