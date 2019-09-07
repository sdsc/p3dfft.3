P3DFFT++ C++ examples
*********************
.. csv-table::
        :widths: auto

        ":download:`test3D_r2c.C <../../sample/C++/test3D_r2c.C>`", "Test program for real-to-complex and complex-to-real 3D FFT. Can be used as a correctness test and a timing test."
        ":download:`test3D_c2c.C <../../sample/C++/test3D_c2c.C>`", "Test program for complex-to-complex forward and backward 3D FFT. Can be used as a correctness test and a timing test."
        ":download:`test_transplan.C <../../sample/C++/test_transplan.C>`", "Test program for real-to-complex and complex-to-real 1D FFT. User can choose an arbitrary dimension of transform. This function performs a local transpose in addition to 1D FFT, as needed. The transpose is specified by input and output memory order, which maps logical onto storage array dimensions. No inter-processor transpose takes place here. 1D transforms can be used by themselves or as a stage in a multidimensional transform."
        ":download:`test1D_cos.C <../../sample/C++/test1D_cos.C>`, :download:`test1D_sin.C <../../sample/C++/test1D_sin.C>`, :download:`test1D_cos_complex.C <../../sample/C++/test1D_cos_complex.C>`", "Same as above, but use cosine and sine transforms, respectively; use cosine transform for complex-numbered inputs."
