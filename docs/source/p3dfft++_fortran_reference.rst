P3DFFT++ Fortran Reference
**************************

.. contents::

Setup and Grid Layout
=====================
In Fortran the grid structure is hidden and is operated on by integer handles.

p3dfft_setup
------------
.. code-block:: fortran

        subroutine p3dfft_setup

*Function*: called once in the beginning of use to initialize P3DFFT++.

p3dfft_cleanup
--------------
.. code-block:: fortran

        subroutine p3dfft_cleanup

*Function*: called once before exit and after the use to free up P3DFFT++ structures.

p3dfft_init_grid
----------------
.. code-block:: fortran

        function p3dfft_init_grid(ldims, glob_start, gdims, dim_conj_sym, pgrid, proc_order, mem_order, mpicomm)

        integer(C_INT) p3dfft_init_grid

*Function*: initialize a new grid

*Arguments*:

.. csv-table:: **IN**
        :widths: auto

        "integer gdims(3)", "Three global grid dimensions (logical order - X, Y, Z)"
        "integer dim_conj_sym", "Dimension of the array in which a little less than half of the elements are omitted due to conjugate symmetry. This argument should be non-negative only for complex-valued arrays resulting from real-to-complex FFT in the given dimension."
        "integer pgrid(3)", "Up to three dimensions of processor grid, decomposing the global grid array. Value =1 means the grid is not decomposed but is local in that logical dimension."
        "integer proc_order(3)", "A permutation of the 3 integers: 0, 1 and 2. Specifies the topology of processor grid on the interconnect. The dimension with lower value means the MPI tasks in that dimension are closer in ranks, e.g. value=0 means the ranks are adjacent (stride=1), value=1 means they are speard out with the stride equal to the pgrid value of the dimension with stride=1 etc"
        "integer mem_order(3)", "A permutation of the 3 integers: 0, 1 and 2. Specifies mapping of the logical dimension and memory storage dimensions for local memory for each MPI task. mem_order(i0) = 0 means that the i0's logical dimension is stored with stride=1 in memory. Similarly, mem_order(i1) =1 means that i1's logical dimension is stored with stride=ldims(i0) etc"
        "integer mpicomm", "The MPI communicator in which this grid lives"

|

.. csv-table:: **OUT**
        :widths: auto

        "integer ldims(3)", "Local dimensions of the grid for each MPI tasks, in order of logical dimensions numbering (XYZ). Essentially ldims = gdims / pgrid."
        "integer glob_start(3)", "Starting coordinates of the local portion of the grid within the global grid."

*Return value*: an integer handle of the initialized grid, to be used later by various routines accessing the grid.

p3dfft_free_grid_f
------------------
.. code-block:: fortran

        subroutine p3dfft_free_grid_f(grid)

*Function*: frees the grid handle

*Arguments*:

.. csv-table:: **IN**
        :widths: auto

        "integer(C_INT) grid", "The handle of the grid to be freed"

One-dimensional (1D) Transforms
===============================
The following predefined 1D transforms are available:

.. csv-table::
        :widths: auto

        "P3DFFT_EMPTY_TYPE", "Empty transform"
        "P3DFFT_R2CFFT_S, P3DFFT_R2CFFT_D", "Real-to-complex forward FFT (as defined in FFTW manual), in single and double precision respectively"
        "P3DFFT_C2RFFT_S, P3DFFT_C2RFFT_D", "Complex-to-real backward FFT (as defined in FFTW manual), in single and double precision respectively"
        "P3DFFT_CFFT_FORWARD_S, P3DFFT_CFFT_FORWARD_D", "Complex forward FFT (as defined in FFTW manual), in single and double precision respectively"
        "P3DFFT_CFFT_BACKWARD_S, P3DFFT_CFFT_BACKWARD_D", "Complex backward FFT (as defined in FFTW manual), in single and double precision respectively"
        "P3DFFT_DCT<x>_REAL_S, P3DFFT_DCT1_REAL_D", "Cosine transform for real-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DCT1, DCT2, DCT3 or DCT4"
        "P3DFFT_DST<x>_REAL_S, P3DFFT_DST1_REAL_D", "Sine transform for real-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DST1, DST2, DST3 or DST4"
        "P3DFFT_DCT<x>_COMPLEX_S, P3DFFT_DCT1_COMPLEX_D", "Cosine transform for complex-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DCT1, DCT2, DCT3 or DCT4"
        "P3DFFT_DST<x>_COMPLEX_S, P3DFFT_DST1_COMPLEX_D", "Sine transform for complex-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DST1, DST2, DST3 or DST4"

1D transform planning
---------------------
.. code-block:: fortran

        function p3dfft_plan_1Dtrans_f(gridIn, gridOut, type, dim, inplace)

        integer p3dfft_plan_1Dtrans

*Function*: defines and plans a 1D transform of a 3D array in a given dimension

*Arguments*:

.. csv-table:: **IN**
        :widths: auto

        "integer gridIn", "Initial grid handle"
        "integer gridOut", "Destination grid handle"
        "integer type", "1D transform type"
        "integer dim", "Dimension rank of the 3D array which should be transformed. valid values are 0, 1 or 2. Note that this is the logical dimension rank (0 for X, 1 for Y, 2 for Z), and may not be the same as the storage dimension, which depends on mem_order member of gridIn and gridOut. The transform dimension of the grid is assumed to be MPI task-local."
        "integer inplace", "Nonzero value if the transform is in-place."

1D transform execution
----------------------
.. code-block:: fortran

        subroutine p3dfft_exec_1Dtrans_single(mytrans,in,out)

        subroutine p3dfft_exec_1Dtrans_double(mytrans,in,out)

*Function*: Executes a 1D transform of a 3D array, in single or double precision

*Arguments*:

.. csv-table:: **IN**
        :widths: auto

        "mytrans", "The handle of a 1D transform predefined earlier with p3dfft_plan_1Dtrans."
        "in", "3D array to be transformed"
        "out", "Destination array (can be the same if inplace was nonzero when defining mytrans)"

.. note::

        1) If inplace was not defined the input and output arrays must be non-overlapping.

        2) This transform is done in the dimension specified in p3dfft_plan_1Dtrans, and this dimension should be local for both input and output arrays.

        3) This subroutine can be called multiple times with the same mytrans and same or different in/out.

Three-dimensional transforms
============================
3D transform planning
---------------------
.. code-block:: fortran

        function p3dfft_plan_3Dtrans_f(gridIn,gridOut,type,inplace, overwrite)

        integer p3dfft_plan_3Dtrans_f

*Function*: defines and plans a 3D transform

*Arguments*:

.. csv-table::
        :widths: auto

        "integer gridIn", "Initial grid handle"
        "integer gridOut", "Destination grid handle"
        "integer type(3)", "Three 1D transform types making up the desired 3D transform"
        "integer inplace", "If nonzero, the transform takes place in-place"
        "integer overwrite", "Nonzero if the input can be overwritten"

*Return value*: a handle of the 3D transform

.. note:: The final grid may or may not be the same as the initial grid. First, in real-to-complex and complex-to-real transforms the global grid dimensions change for example from (n0,n1,n2) to (n0/2+1,n1,n2), since most applications attempt to save memory by using the conjugate symmetry of the Fourier transform of real data. Secondly, the final grid may have different processor distribution and memory ordering, since for example many applications with convolution and those solving partial differential equations do not need the initial grid configuration in Fourier space. The flow of these applications is typically 1) transform from physical to Fourier space, 2) apply convolution or derivative calculation in Fourier space, and 3) inverse FFT to physical space. Since forward FFT's last step is 1D FFT in the third dimension, it is more efficient to leave this dimension local and stride-1, and since the first step of the inverse FFT is to start with the third dimension 1D FFT, this format naturally fits the algorithm and results in big savings of time due to elimination of several extra transposes.

3D transform execution
----------------------
.. code-block:: fortran

        subroutine p3dfft_exec_3Dtrans_single(mytrans,in,out)

        subroutine p3dfft_exec_3Dtrans_double(mytrans,in,out)

*Function*: Executes a predefined 3D transform in single or double precision

*Arguments*:

.. csv-table::
        :widths: auto

        "mytrans", "The handle of the predefined 3D transform"
        "in", "Input array"
        "out", "Output array"

.. note:: this subroutine can be called multiple times for the same mytrans and same or different in/out.Input and output arrays are local portions of the global 3D array, assumed to be stored contiguously in memory following the definition of the grids in planning stage.

3D Spectral Derivative
----------------------
.. code-block:: fortran

        p3dfft_exec_3Dtrans_single(mytrans, In, Out, idir)

        p3dfft_exec_3Dtrans_double(mytrans, In, Out, idir)

*Function*: execute 3D real-to-complex FFT, followed by spectral derivative calculation, i.e. multiplication by (ik), where i is the complex imaginary unit, and k is the wavenumber; in single or double precision, respectively

*Arguments*:

.. csv-table::
        :widths: auto

        "In and Out", "Input and output arrays, assumed to be the local portion of the 3D grid array stored contiguously in memory, consistent with definition of Grid in planning stage."
        "integer idir", "The dimension where derivative is to be taken in (this is logical dimension, NOT storage mapped). Valid values are 0 - 2."

.. note::

        1) Unless inplace was defined in the planning stage of mytrans, In and Out must be non-overlapping

        2) These functions can be used multiple times after the 3D transform has been defined and planned.
