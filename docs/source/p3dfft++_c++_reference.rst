P3DFFT++ C++ Reference
**********************

.. contents::

Introduction
============
For C++ users all P3DFFT++ objects are defined within the p3dfft namespace, in order to avoid confusion with user-defined objects. For example, to initialize P3DFFT++ it is necessary to call the function p3dfft::setup(), and to exit P3DFFT++ one should call p3dfft::cleanup() (alternatively, one can use namespace p3dfft and call setup() and cleanup()). From here on in this document we will omit the implicit p3dfft:: prefix from all C++ names.

Setup and Grid layout
=====================
The public portion of the grid class is below:

.. code-block:: cpp

        class grid {

        ...

        public :

        int taskid,numtasks;

        int nd; //number of dimensions the volume is split over

        int gdims[3]; //Global dimensions

        dim_conj_sym; // dimension of the array in which a little less than half of the elements are omitted due to conjugate symmetry. This argument should be non-negative only for complex-valued arrays resulting from real-to-complex FFT in the given dimension.

        int mem_order[3]; //Memory ordering inside the data volume
        int ldims[3]; //Local dimensions on THIS processor
        int pgrid[3]; //Processor grid
        int proc_order[3]; //Ordering of tasks in processor grid, e.g. (1,2,3) : first dimension - adjacent tasks,then second, then third dimension
        int P[3]; //Processor grid size (in inverse order of split dimensions, i.e. rows first, then columns etc
        int D[3]; //Ranks of Dimensions of physical grid split over rows and columns correspondingly
        int L[3]; //Rank of Local dimension (p=1)
        int grid_id[3]; //Position of this pencil/cube in the processor grid
        int grid_id_cart[3];
        int glob_start[3]; // Starting coords of this cube in the global grid
        MPI_Comm mpi_comm_glob; // Global MPi communicator we are starting from
        MPI_Comm mpi_comm_cart;
        MPI_Comm mpicomm[3]; //MPI communicators for each dimension
        // int (*st)[3],(*sz)[3],(*en)[3]; // Lowest, size and uppermost location in 3D, for each processor in subcommunicator
        int **st[3],**sz[3],**en[3]; // Lowest, size and uppermost location in 3D, for each processor in subcommunicator 

        bool is_set;
        grid(int gdims_[3],int pgrid_[3],int proc_order_[3],int mem_order[3],
        MPI_Comm mpicomm_);
        grid(const grid &rhs);
        grid() {};
        ~grid();
        void set_mo(int mo[3]) {for(int i=0;i<3;i++) mem_order[i] = mo[i];};

        ...
        };

*grid* constructor
------------------
.. code-block:: cpp

        grid(int gdims[3],int dim_conj_sym,int pgrid[3],int proc_order[3],int mem_order[3],MPI_Comm mpicomm);

*Function*:Initializes a new grid with specified parameters.

*Arguments*:

.. csv-table::
        :widths: auto

        "gdims", "Three global grid dimensions (logical order - X, Y, Z)"
        "dim_conj_sym", "Dimension of conjugate symmetry, non-negative only for complex arrays resulting from real-to-complex FFT in the given dimension. This is logical, not storage, dimension, with valid numbers 0 - 2, and -1 implying no conjugate symmetry."
        "pgrid", "Up to three dimensions of processor grid, decomposing the global grid array. Value =1 means the grid is not decomposed but is local in that logical dimension."
        "proc_order", "A permutation of the 3 integers: 0, 1 and 2. Specifies the topology of processor grid on the interconnect. The dimension with lower value means the MPI tasks in that dimension are closer in ranks, e.g. value=0 means the ranks are adjacent (stride=1), value=1 means they are speard out with the stride equal to the pgrid value of the dimension with stride=1 etc"
        "mem_order", "A permutation of the 3 integers: 0, 1 and 2. Specifies mapping of the logical dimension and memory storage dimensions for local memory for each MPI task. mem_order[i0] = 0 means that the i0's logical dimension is stored with stride=1 in memory. Similarly, mem_order[i1] =1 means that i1's logical dimension is stored with stride=ldims[i0] etc"
        "mpicomm", "The MPI communicator in which this grid lives"

P3DFFT++ Transforms
===================
P3DFFT++ functions in a way similar to FFTW: first the user needs to plan a transform, using a planner function once per each transform type. The planner function initializes the transform, creates a plan and stores all information relevant to this transform inside P3DFFT++. The users gets a handle referring to this plan (which is a class in C++) that can be later used to execute this transform, and can be applied multiple times. The handles can be released after use.

In order to define and plan a transform (whether 1D or 3D) one needs to first define initial and final grid objects. They contain all the necessary grid decomposition parameters. P3DFFT++ figures out the optimal way to transpose the data between these two grid configurations, assuming they are consistent (i.e. same grid size, number of tasks etc).

One-Dimensional (1D) Transforms
-------------------------------
The following predefined 1D transforms are available:

.. csv-table::
        :widths: auto

        "EMPTY_TYPE", "Empty transform"
        "R2CFFT_S, P3DFFT_R2CFFT_D", "Real-to-complex forward FFT (as defined in FFTW manual), in single and double precision respectively"
        "C2RFFT_S, P3DFFT_C2RFFT_D", "Complex-to-real backward FFT (as defined in FFTW manual), in single and double precision respectively"
        "CFFT_FORWARD_S, CFFT_FORWARD_D", "Complex forward FFT (as defined in FFTW manual), in single and double precision respectively"
        "CFFT_BACKWARD_S, CFFT_BACKWARD_D", "Complex backward FFT (as defined in FFTW manual), in single and double precision respectively"
        "DCT<x>_REAL_S, DCT1_REAL_D", "Cosine transform for real-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DCT1, DCT2, DCT3 or DCT4"
        "DST<x>_REAL_S, DST1_REAL_D", "Sine transform for real-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DST1, DST2, DST3 or DST4"
        "DCT<x>_COMPLEX_S, DCT1_COMPLEX_D", "Cosine transform for complex-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DCT1, DCT2, DCT3 or DCT4"
        "DST<x>_COMPLEX_S, DST1_COMPLEX_D", "Sine transform for complex-numbered data, in single and double precision, where <x> stands for the variant of the cosine transform, such as DST1, DST2, DST3 or DST4"

Custom transform types
^^^^^^^^^^^^^^^^^^^^^^
Custom 1D transforms can be defined by the user through trans_type1D class template.

.. code-block:: cpp

        template <class Type1,class Type2> class trans_type1D : public gen_trans_type{

        int ID;
        public :

        typedef long (*doplan_type)(const int *n,int howmany,Type1 *in,const int *inembed,int istride,int idist,Type2 *out,const int *onembed,int ostride,int odist,...);

        long (*doplan)(...);
        void (*exec)(...);

        trans_type1D(const char *name, long (*doplan_)(...),void (*exec)(...)=NULL,int isign=0);
        inline int getID() {return(ID);}
        trans_type1D(const trans_type1D &rhs); 
        ~trans_type1D();

        };

This class template is a derivative of gen_trans_type1D class, defined as follows:

.. code-block:: cpp

        class gen_trans_type {
        public :
        char *name;
        int isign; // forward (-1) or backward (+1), in case this is complex FFT
        bool is_set,is_empty;
        int dt1,dt2; //Datatype before and after
        int prec; // precision for a real value in bytes (4 or 8)
        gen_trans_type(const char *name_,int isign_=0);
        ~gen_trans_type();
        bool operator==(const gen_trans_type &) const;
        };

In order to define a custom transform type, the user needs to provide planning and execution functions (doplan and exec).  For example, in case of a complex FFT implemented through FFTW, the following is how the transform type is constructed:

.. code-block:: cpp

        char *name = "Complex-to-complex Fourier Transform, forward transform, double precision";
        int isign = FFTW_FORWARD;
        trans_type1D<complex_double,complex_double> *mytype = new trans_type1D<complex_double,complex_double>(name,(long (*)(...) ) fftw_plan_many_dft,(void (*)(...)) exec_c2c_d,isign);

where exec_c2c_d is defined as follows:

.. code-block:: cpp

        void exec_c2c_d(long plan,complex_double *in,complex_double *out)
        {
        fftw_execute_dft((fftw_plan) plan,(fftw_complex *) in,(fftw_complex *) out);
        }

Planning 1D transform 
^^^^^^^^^^^^^^^^^^^^^
1D transform in C++ is realized through transplan template class. TypeIn and TypeOut are the datatypes for input and output.

Two constructors are provided.

.. code-block:: cpp

        template <class TypeIn,class TypeOut> class transplan::transplan(const grid &gridIn,const grid &gridOut,const gen_trans_type *type,const int d, const bool inplace_);

        template <class TypeIn,class TypeOut> class transplan::transplan(const grid &gridIn,const grid &gridOut,const int type,const int d, const bool inplace_);

*Function*: define and plan a 1D transform of a 3D array

*Arguments*:

.. csv-table::
        :widths: auto

        "gridIn", "Initial grid descriptor"
        "gridOut", "Final grid descriptor"
        "type", "The type of the 1D transform (either as a predefined integer parameter, or as a class gen_trans_type."
        "d", "The dimension to be transformed. Note that this is the logical dimension rank (0 for X, 1 for Y, 2 for Z), and may not be the same as the storage dimension, which depends on mem_order member of gridIn and gridOut. The transform dimension of the grid is assumed to be MPI task-local."
        "inplace", "True for in-place transform, false for out-of-place."

Releasing 1D transform handle
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
To release a 1D transform handle, simply delete the corresponding transplan class.

Executing 1D transform
^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: cpp

        template <class TypeIn,class TypeOut> class transplan::exec(char *In, char *Out);

*Function*: Executes the pre-planned 1D transform of a 3D array

*Arguments*:

.. csv-table::
        :widths: auto

        "In and Out", "Pointers to input and output arrays, cast as pointers to char. They contain the local portion of the 3D input and output arrays, arranged as a contiguous sequence of numbers according to local grid dimensions and the memory order of initial and final grid objects respectively."

.. note:: If the transform is out-of-place, then these arrays must be non-overlapping. The execution can be performed many times with the same handle and same or different input and output arrays.

Three-dimensional Transforms
----------------------------
Three-dimensional (3D) transforms consist of three one-dimensional transforms in sequence (one for each dimension), interspersed by inter-processor transposes. In order to specify a 3D transform, three main things are needed:

1. Initial grid (as described above, grid object defines all of the specifics of grid dimensions, memory ordering and distribution among processors).
2. Final grid.
3. The type of 3D transform.

The final grid may or may not be the same as the initial grid. First, in real-to-complex and complex-to-real transforms the global grid dimensions change for example from (n0,n1,n2) to (n0/2+1,n1,n2), since most applications attempt to save memory by using the conjugate symmetry of the Fourier transform of real data. Secondly, the final grid may have different processor distribution and memory ordering, since for example many applications with convolution and those solving partial differential equations do not need the initial grid configuration in Fourier space. The flow of these applications is typically 1) transform from physical to Fourier space, 2) apply convolution or derivative calculation in Fourier space, and 3) inverse FFT to physical space. Since forward FFT's last step is 1D FFT in the third dimension, it is more efficient to leave this dimension local and stride-1, and since the first step of the inverse FFT is to start with the third dimension 1D FFT, this format naturally fits the algorithm and results in big savings of time due to elimination of several extra transposes.

In order to define the 3D transform type one needs to know three 1D transform types comprising the 3D transform. In C++ 3D transform type is interfaced through a class trans_type3D.

trans_type3D constructor
^^^^^^^^^^^^^^^^^^^^^^^^
Two constructors are provided for trans_type3D (in addition to a copy constructor):

.. code-block:: cpp

        trans_type3D::trans_type3D(const gen_trans_type *types_[3]); 
        trans_type3D::trans_type3D(const int types[3]);

Types is an array of 3 1D transform types, either as integer type IDs, or gen_trans_type classes.

trans_type3D class has the following public members:

.. code-block:: cpp

        char *name;
        int dtIn,dtOut; // Datatypes for input and output: 1 is real, 2 is complex
        int prec; // Datatype precision for a real value in bytes: 4 for single, 8 for double precision

        bool is_set;
        int types[3]; // 3 1D transform types

Transform3D constructor
^^^^^^^^^^^^^^^^^^^^^^^
In C++ 3D transforms are handled through class template transform3D, with input and output datatypes TypeIn and TypeOut. Often these will be the same, however some transforms have different types on input and output, for example real-to-complex FFT. In all cases the floating point precision (single/double) of the initial and final types should match.

.. code-block:: cpp

        template<class TypeIn,class TypeOut> class transform3D::transform3D( const grid &grid_in, const grid &grid_out, const trans_type3D *type, const bool inplace, const bool Overwrite);

*Function*: defines and plans a 3D transform

*Arguments*:

.. csv-table::
        :widths: auto

        "gridIn", "Initial grid configuration"
        "gridOut", "Final grid configuration"
        "type", "pointer to a 3D transform type class"
        "inplace", "true is this is an in-place transform; false if an out-of-place transform."
        "Overwrite (optional)", "Indicates whether input can be overwritten (true=yes, default is no)"

Transform3D Execution
^^^^^^^^^^^^^^^^^^^^^
.. code-block:: cpp

        template<class TypeIn,class TypeOut> class transform3D::exec(TypeIn *In,TypeOut *Out);

*Function*: executes a 3D transform

*Arguments*:

.. csv-table::
        :widths: auto

        "In and Out", "Pointers to input and output arrays. In case of in-place transform they can point to the same location. For out-of-place transforms the arrays must be non-overlapping."

Spectral Derivative for 3D array
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. code-block:: cpp

        template<class TypeIn,class TypeOut> class transform3D::exec_deriv(TypeIn *In,TypeOut *Out, int idir);

*Function*: execute 3D real-to-complex FFT, followed by spectral derivative calculation, i.e. multiplication by (ik), where i is the complex imaginary unit, and k is the wavenumber. This function is defined only for complex-valued output arrays (single or double precision), i.e. TypeOut must be either mycomplex or complex_double.

*Arguments*:

.. csv-table::
        :widths: auto

        "In and Out", "Pointers to input and output arrays, assumed to be the local portion of the 3D grid array stored contiguously in memory, consistent with definition of grids in planning stage."
        "idir", "The dimension where derivative is to be taken in (this is logical dimension, NOT storage mapped). Valid values are 0 - 2."

.. note::

        1) Unless inplace was defined in the planning stage of mytrans, In and Out must be non-overlapping

        2) This function can be used multiple times after the 3D transform has been defined and planned.
