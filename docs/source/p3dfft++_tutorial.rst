P3DFFT++ Tutorial
*****************

General considerations
======================
P3DFFT++ is written in C++ and contains wrappers providing easy interfaces with C and Fortran. 

For C++ users all P3DFFT++ objects are defined within the p3dfft namespace, in order to avoid confusion with user-defined objects. For example, to initialize P3DFFT++ it is necessary to call the function p3dfft::setup(), and to exit P3DFFT++ one should call p3dfft::cleanup() (alternatively, one can use namespace p3dfft and call setup() and cleanup()). From here on in this document we will omit the implicit p3dfft:: prefix from all C++ names. 

In C and Fortran these functions become p3dfft_setup and p3dfft_cleanup.  While C++ users can directly access P3DFFT objects such as grid class, C and Fortran users will access these through handles provided by corresponding wrappers (see more details below). 

Data types
==========
P3DFFT++ currently operates on four main data types:

1. float (single precision floating point)
2. double (double precision floating point)
3. mycomplex (single precision complex number) (equivalent to complex<float>)
4. complex_double (double precision complex number) (equivalent to complex<double>)

Data layout
===========
While P3DFFT had the assumption of predetermined 2D pencils in X and in Z dimensions as the primary data storage, P3DFFT++ relaxes this assumption to include more general formats, such as arbitrary shape and memory order 2D pencils as well as 3D blocks. Below is the technical description of how to specify the data layout formats. 

A basic P3DFFT++ descriptor is the "grid" construct. It defines all necessary information about decomposition of a grid among parallel tasks/processors. In C++ it is defined as a class, while in C and in Fortran it is defined through handles to a C++ object through inter-language wrappers. Below is the technical description of the definition for each language.

C++
---
The following is the main constructor call for the *grid* class:

.. code-block:: cpp

        grid(int gdims[3],int dim_conj_sym, int pgrid[3],int proc_order[3],int mem_order[3],MPI_Comm mpicomm);

*Arguments*:

.. csv-table::
        :widths: auto

        "*gdims*", "The three global dimensions of the grid to be decomposed. Here the Fortran-inspired convention is followed: the first of the three numbers specifies the dimension with the fastest changing index, i.e. the first logical dimension (X)."
        "*dim_conj_sym*", "The dimension of conjugate symmetry where we store N/2+1 of the data after Real-to-complex transform due to conjugate symmety;(-1 for none)"
        "*pgrid*", "The processor grid to be used in decomposition. For example, a 2D pencil with the first dimension local (X-pencil) would be described as having pgrid={1,P1,P2}, where P1 and P2 are the dimensions of 2D decomposition such that P1 x P2 = P, the total number of tasks. Of course a 2D grid could be defined as a Y-pencil (pgrid={P1,1,P2}) or a Z pencil (P1,P2,1). 1D decomposition (slabs) would be defined as (1,1,P), or (1,P,1) or (P,1,1), depending on the orientation of the slabs. 3D decomposition is also possible where each of the three values of pgrid is greater than 1."
        "*proc_order*", "The ordering of the processor grid dimensions with respect to the default layout of MPI tasks. For example, the simplest ordering where proc_order={0,1,2} and pgrid={1,P1,P2} corresponds to a grid with the second dimension decomposed among P1 MPI tasks adjacent to each other in the default MPI topology, such as tasks on the same node and/or neighboring nodes on the network, while the third dimension would be decomposed among P2 tasks on non-neighboring nodes (with stride equal to P1). On the other hand, proc_order={0,2,1} with the same pgrid would correspond to P2 tasks splitting the third dimension of the grid being adjacent in the MPI/network neighborhood, while the second dimension would be split over P1 tasks distant in the topology with stride equal to P2."
        "*mem_order*", "The relative ordering of the three dimensions in memory within the local portion of the grid. Here C-style indexing is used (indices start with 0). The simplest ordering {0,1,2} corresponds to the first logical dimension being stored with the fastest changing index (memory stride=1), followed by the second (stride=L0) and the third dimension (stride=L0*L1), where Li is the size of local grid in i's dimension for a given MPI task. This corresponds to a C array A[L2][L1][L0]. As another example, ordering {2,0,1} means that the second dimension (L1) is stored with the fastest-changing index (memory stride=1), the third dimension dimension (L2) with the medium stride =L1, and the first dimension is stored with the slowest index, stride=L1*L2. This would correspond to a C array A[L0][L2][L1]."
        "*mpicomm*", "The initial communicator for all subsequent library operations. It is recommended that the users define or duplicate a communicator for P3DFFT++ to be different from the one(s) used in their code, in order to avoid interference."

For example:

.. code-block:: cpp

        main() {
        
        ...
        
        int gdims[3],pgrid[3],proc_order[3], mem_order[3];
        
        MPI_Comm mpicomm;
        
        ...
        
        gdims= {128, 128, 128};
        
        pgrid={1,4,4}; //X-pencil
        
        proc_order = {0,1,2};
        
        mem_order={0,1,2};
        
        MPI_Comm_dup(MPI_COMM_WORLD, &mpicomm);
        
        grid mygrid(gdims, -1, pgrid, proc_order, mem_order, mpicomm);
        
        }

Upon construction the *grid* object defines several useful parameters, available by accessing the following public class members of *grid*:

.. csv-table::
        :widths: auto

        "*int ldims[3]*", "Dimensions of the local portion of the grid (ldims[0]=gdims[0]/pgrid[0] etc). Note: these dimensions are specified in the order of logical grid dimensions and may differ from memory storage order, which is defined by *mem_order*."
        "*int nd*", "Number of dimensions of the processor grid (1, 2 or 3)."
        "*int L[3]*", "0 to 3 local dimensions (i.e. not split)."
        "*int D[3]*", "0 to 3 split dimensions."
        "*int glob_start[3]*", "Coordinates of the lowest element of the local grid within the global array. This is useful for reconstructing the global grid from grid pieces for each MPI task."

and other useful information.  The grid class also provides a copy constructor. 

To release a grid object, simply delete it. 

C
^
For C users grid initialization is accomplished by a call to p3dfft_init_grid, returning a pointer to an object of type *Grid*. This type is a C structure containing a large part of the C++ class *grid*. Calling p3dfft_init_grid initializes the C++ *grid* object and also copies the information into a *Grid* object accessible from C, returning its pointer. For example:

.. code-block:: c

        int xdim;

        Grid *grid1;

        grid1 = p3dfft_init_grid(gdims, dim_conj_sym, pgrid, proc_order, mem_order, mpicomm);

        xdim = grid1->ldims[0]; /* Size of zero logical dimension of the local portion of the grid for a given processor */

        To release a grid object simply execute 

        p3dfft_free_grid(Grid *gr);

Fortran
-------
For Fortran users the grid object is represented as a handle of type *integer(C_INT)*. For example:

.. code-block:: fortran

        integer(C_INT) grid1

        integer ldims(3),glob_start(3),gdims(3),dim_conj_sym,pgrid(3),proc_order(3),mem_order(3),mpicomm

        grid1 = p3dfft_init_grid(ldims, glob_start, gdims, dim_conj_sym, pgrid, proc_order, mem_order, mpicomm)

This call initializes a C++ grid object as a global variable and assigns an integer ID, returned in this example as *grid1*. In addition this call also returns the dimensions of the local portion of the grid (*ldims*) and the position of this portion within the global array (*glob_start*).

Other elements of the C++ grid object can be accessed through respective functions, such as p3dfft\_grid_get_...

To release a grid object, simply call:

.. code-block:: fortran

        *p3dfft_free_grid_f(gr)*

where *gr* is the grid handle. 

P3DFFT++ Transforms
===================
P3DFFT++ aims to provide a versatile toolkit of algorithms/transforms in frequent use for solving multiscale problems. To give the user maximum flexibility there is a range of algorithms from top-level algorithms operating on the entire 3D array, to 1D algorithms which can function as building blocks the user can arrange to suit his/her needs. In addition, inter-processor exchanges/transposes are provided, so as to enable the user to rearrange the data from one orientation of  pencils to another, as well as other types of exchanges. In P3DFFT++ the one-dimensional transforms are assumed to be expensive in terms of memory bandwidth, and therefore such transforms are performed on local data (i.e. in the dimension that is not distributed across processor grid). Transforms in three dimensions consist of three transforms in one dimension, interspersed by inter-processor interchange as needed to rearrange the data.  The 3D transforms are  high-level functions saving the user work in arranging the 1D transforms and transposes, as well as often providing superior performance. **We recommend to use 3D transforms whenever they fit the user's algorithm.**

Although syntax for C++, C and Fortran is different, using P3DFFT++ follows the same logic. P3DFFT++ functions in a way similar to FFTW: first the user needs to plan a transform, using a planner function once per each transform type. The planner function initializes the transform, creates a plan and stores all information relevant to this transform inside P3DFFT++. The users gets a handle referring to this plan (the handle is a class in C++, and an integer variable in C or Fortran) that can be later used to execute this transform, which can be applied multiple times. The handles can be released after use.

In order to define and plan a transform (whether 1D or 3D, in C++, C or Fortran) one needs to first define initial and final grid objects. They contain all the necessary grid decomposition parameters. P3DFFT++ figures out the optimal way to transpose the data between these two grid configurations, assuming they are consistent (i.e. same grid size, number of tasks etc).

One-dimensional (1D) Transforms
===============================
1D transforms is the smaller building block for higher dimensional transforms in P3DFFT++. They include different flavors of Fast Fourier Transforms (FFTs), empty transform (provided for convenience, as in the case where a user might want to implement their own 1D transform, but is interested in memory reordering to arrange the transform dimension for stride-1 data access), and (in the future) other transforms that share the following property: they are memory bandwidth and latency intensive,  and are optimally done when the dimension the transform operates on is entirely within one MPI task's domain. 

1D transforms can be done with or without data exchange and/or memory reordering. In general, combining a transform with an exchange/reordering can be beneficial for performance due to cache reuse, compared to two separate calls to a transform and an exchange. 

The following predefined 1D transforms are available (in C++ the P3DFFT\_ prefix can be omitted if used within P3DFFT namespace).

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

C++
---
Below is an example of how a 1D transform can be called from C++. In this example, real-to-complex transform in double precision is planned and then performed. First a constructor for class transplan is called:

.. code-block:: cpp

        transplan<double,complex_double> trans_f(gridIn, gridOut, R2C_FFT_D, dim, false);

Here *gridIn* and *gridOut* are initial and final *grid* objects, describing, among other things, initial and final memory ordering of the grid storage array (ordering can be the same or different for input and output). *dim* is the dimension/rank to be transformed. Note that this is the logical dimension rank (0 for X, 1 for Y, 2 for Z), and may not be the same as the storage dimension, which depends on *mem_order* member of *gridIn* and *gridOut*. The transform dimension of the grid is assumed to be MPI task-local. The second last parameter is a bool variable telling P3DFFT++ whether this is an in-place or out-of-place transform. Note that in C++ the P3DFFT\_ prefix for transform types is optional. 

When a *transplan* constructor is called as above, P3DFFT++ stores the parameters of the 1D transform and if needed, plans its execution (i.e. as in FFTW planning) and stores the plan handle. This needs to be done once per transform type. In order to execute the transform, simply call exec member of the class, e.g.:

.. code-block:: cpp

        trans_f.exec((char *) In,(char *) Out);

Here *In* and *Out* are pointers to input and output arrays. In this case they are of type *double* and *complex_double*, however in this call they are cast as *char**, as required by P3DFFT++. They contain the local portion of the 3D input and output arrays, arranged as a contiguous sequence of numbers according to local grid dimensions and the memory order of *gridIn* and *gridOut* classes, respectively. If the transform is out-of-place, then these arrays must be non-overlapping. The execution can be performed many times with the same handle and same or different input and output arrays.This call will perform the 1D transform specified when the *transplan* object was constructed, along the dimension *dim*. Again, the logical dimension specified as *dim* in the planning stage must be MPI-local for both input and output arrays. Other utilities allow the user to transpose the grid arrays in MPI/processor space (*see MPIplan and transMPIplan*).

To release the transform handle simply delete the transplan class object. 

C
-
Here is an example of initializing and executing a 1D transform (again, a real-to-complex double precision FFT) in a C program.

.. code-block:: c

        Grid *gridIn, *gridOut;

        Plan3D trans_f;

        ...

        gridIn = p3dfft_init_grid(gdimsIn, pgridIn, proc_order, mem_orderIn, MPI_COMM_WORLD);
        gridOut = p3dfft_init_grid(gdimsOut, pgridOut, proc_order, mem_orderOut, MPI_COMM_WORLD);

        trans_f = p3dfft_plan_1Dtrans(gridIn, gridOut, P3DFFT_R2CFFT_D, dim, 0);

Here *gridIn* and *gridOut* are pointers to the C equivalent of P3DFFT++ *grid* object (initial and final), *trans_f* is the handle for the 1D transform after it has been initialized and planned, *dim* is the logical dimension of the transform (0, 1 or 2), and the last argument indicates that this is not an in-place transform (a non-zero argument would indicate in-place). This initialization/planning needs to be done once per transform type.

.. code-block:: c

        p3dfft_exec_1Dtrans_double(trans_f,IN,OUT);

This statement executes the 1D transformed planned and handled by *trans_f*. *IN* and *OUT* are pointers to one-dimensional input and output arrays containing the 3D grid stored contiguously in memory based on the local grid dimensions and storage order of *gridIn* and *gridOut*. The execution can be performed many times with the same handle and same or different input and output arrays. In case of out-of-place transform the input and output arrays must be non-overlapping. 

Fortran
-------
Here is an example of initializing and executing a 1D transform (again, a real-to-complex double precision FFT) in a Fortran program:

.. code-block:: fortran

        integer(C_INT) gridIn,gridOut
        integer trans_f

        gridIn = p3dfft_init_grid(ldimsIn, glob_startIn, gdimsIn, pgridIn, proc_order, mem_orderIn, MPI_COMM_WORLD)
        gridOut = p3dfft_init_grid(ldimsOut, glob_startOut, gdimsOut, pgridOut, proc_order, mem_orderOut, MPI_COMM_WORLD)
        trans_f = p3dfft_plan_1Dtrans_f(gridIn, gridOut, P3DFFT_R2CFFT_D, dim-1, 0)

These statement set up initial and final grids (gridIn and gridOut), initialize and plan the 1D real-to-complex double FFT and use trans_f as its handle. This needs to be done once per transform type. Note that we need to translate the transform dimension dim into C convention (so that X corresponds to 0, Y to 1 and Z to 2). The last argument is 0 for out-of-place and non-zero for in-place transform.

.. code-block:: fortran

        call p3dfft_1Dtrans_double(trans_f,Gin,Gout)

This statement executes the 1D transform planned before and handled by trans_f. Gin and Gout are 1D contiguous arrays of values (double precision and double complex) of the 3D grid array, according to the local grid dimensions and memory storage order of gridIn and gridOut, respectively. After the previous planning step is complete, the execution can be called many times with the same handle and same or different input and output arrays. If the transform was declared as out-of-place then Gin and Gout must be non-overlapping.

Three-dimensional Transforms
============================
As mentioned above, three-dimensional (3D) transforms consist of three one-dimensional transforms in sequence (one for each dimension), interspersed by inter-processor transposes. In order to specify a 3D transform, five main things are needed:

1. Initial *grid* (as described above, *grid* object defines all of the specifics of grid dimensions, memory ordering and distribution among processors).
2. Final *grid*.
3. The type of 3D transform.
4. Whether this is in-place transform
5. Whether this transform can overwrite input

The final grid may or may not be the same as the initial grid. First, in real-to-complex and complex-to-real transforms the global grid dimensions change for example from (n0,n1,n2) to (n0/2+1,n1,n2), since most applications attempt to save memory by using the conjugate symmetry of the Fourier transform of real data. Secondly, the final grid may have different processor distribution and memory ordering, since for example many applications with convolution and those solving partial differential equations do not need the initial grid configuration in Fourier space. The flow of these applications is typically 1) transform from physical to Fourier space, 2) apply convolution or derivative calculation in Fourier space, and 3) inverse FFT to physical space. Since forward FFT's last step is 1D FFT in the third dimension, it is more efficient to leave this dimension local and stride-1, and since the first step of the inverse FFT is to start with the third dimension 1D FFT, this format  naturally fits the algorithm and results in big savings of time due to elimination of several extra transposes. 

In order to define the 3D transform type one needs to know three 1D transform types comprising the 3D transform. Usage of 3D transforms is different depending on the language used and is described below.

C++
---
In C++ 3D transform type is interfaced through a class trans_type3D, which is constructed as in the following example:

.. code-block:: cpp

        trans_type3D name_type3D(int types1D[3]);

Here *types1D* is the array of three 1D transform types which define the 3D transform (empty transforms are permitted). Copy constructor is also provided for this class.

For example:

.. code-block:: cpp

        int type_rcc, type_ids[3];

        type_ids[0] = P3DFFT_R2CFFT_D;
        type_ids[1] = P3DFFT_CFFT_FORWARD_D;
        type_ids[2] = P3DFFT_CFFT_FORWARD_D;

        trans_type3D mytype3D(type_ids);

3D transforms are provided as the class template:

.. code-block:: cpp

        template<class TypeIn,class TypeOut> class transform3D;

Here *TypeIn* and *TypeOut* are initial and final data types. Most of the times these will be the same, however some transforms have different types on input and output, for example real-to-complex FFT. In all cases the floating point precision (single/double) of the initial and final types should match. 

The constructor of transform3D takes the following arguments:

.. code-block:: cpp

        transform3D<TypeIn,TypeOut>  my_transform_name(gridIn,gridOut,type,inplace,overwrite);

Here type is a 3D transform type (constructed as shown above), inplace is a bool variable indicating whether this is an in-place transform, and overwrites (also boolean) defines if the input can be rewritten (default is false). *gridIn* and *gridOut* are initial and final grid objects. Calling a *transform3D* constructor creates a detailed step-by-step plan for execution of the 3D transform and stores it in the *my_transform_name* object. 

Once a 3D transform has been defined and planned, execution of a 3D transform can be done by calling:

.. code-block:: cpp

        my_transform_name.exec(TypeIn *in,TypeOut *out);

Here *in* and *out* are initial and final data arrays of appropriate types. These are assumed to be one-dimensional contiguous arrays containing the three-dimensional grid for input and output, local to the memory of the given MPI task, and stored according to the dimensions and memory ordering specified in the *gridIn* and *gridOut* objects, respectively.  For example, if grid1.ldims={2,2,4} and grid1.mem_order={2,1,0}, then the in array will contain the following sequence: G000, G001, G002, G003, G010, G011, G012, G013, G100, G101, G102, G103, G110, G111, G112, G113. Again, we follow the Fortran convention that the fastest running index is the first, (i.e. G012 means the grid element at X=0, Y=1, Z=2).   

C
^
In C a unique datatype Type3D is used to define the 3D transform needed. *p3dfft_init_3Dtype* function is used to initialize a new 3D transform type, based on the three 1D transform types, as in the following example:

.. code-block:: c

        int type_rcc,  type_ids[3];

        type_ids[0] = P3DFFT_R2CFFT_D;
        type_ids[1] = P3DFFT_CFFT_FORWARD_D;
        type_ids[2] = P3DFFT_CFFT_FORWARD_D;

        type_rcc = p3dfft_init_3Dtype(type_ids);

In this example type_rcc will describe the real-to-complex (R2C) 3D transform (R2C in 1D followed by two complex 1D transforms).

To define and plan the 3D transform, use p3dfft_plan_3Dtrans function as follows:

.. code-block:: c

        int mytrans;

        mytrans = p3dfft_plan_3Dtrans(gridIn,gridOut,type,inplace,overwrite);

Here *gridIn* and *gridOut* are pointers to initial and final grid objects (of type *Grid*); *type* is the 3D transform type defined as above; *inplace* is an integer indicating an in-place transform if it's non-zero, out-of-place otherwise. Overwrite is an integer defining if the input can be overwritten (non-zero; default is zero). In this example *mytrans* contains the handle to the 3D transform that can be executed (many times) as follows:

.. code-block:: c

        p3dfft_exec_3Dtrans_double(mytrans,in,out);

Here *in* and *out* are pointers to input and output arrays, as before, assumed to be the local portion of the 3D grid array stored according to *gridIn* and *gridOut* descriptors. For single precision use *p3dfft_exec_3Dtrans_single*.

Fortran
-------
In Fortran, similar to C, to define a 3D transform the following routine is used:

.. code-block:: fortran

        mytrans = p3dfft_plan_3Dtrans_f(gridIn,gridOut,type,inplace, overwrite)

Here *gridIn* and *gridOut* are handles defining the initial and final grid configurations; *type* is the 3D transform type, defined as above; and *inplace* is the integer whose non-zero value indicates this is an in-place transform (or 0 for out-of-place). Non-zero overwrite indicates it is OK to overwrite input (default is no). Again, this planner routine is called once per transform. Execution can be called multiple times as follows:

.. code-block:: fortran

        call p3dfft_3Dtrans_double(mytrans,IN,OUT)

Here *IN* and *OUT* are the input and output arrays. For single precision use *p3dfft_3Dtrans_single_f*.
