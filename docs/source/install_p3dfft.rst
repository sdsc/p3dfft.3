.. _installing_p3dfft:

Installing P3DFFT
=================
The latest version of P3DFFT can be found `here <https://github.com/sdsc/p3dfft/releases/latest>`_. Once you have extracted the package, you must take the following steps to complete the setup:

1. Run the configure script
2. Run 'make'
3. Run 'make install'

How to compile P3DFFT
---------------------
P3DFFT uses a configure script to create Makefiles for compiling the library as well as several examples. This page will step you through the process of running the configure script properly.

Run "configure --help" for complete list of options. Recommended options: --enable-stride1. For Cray XT platforms also recommended --enable-useeven.

Currently the package supports four compiler suites: PGI, Intel, IBM and GNU. Some examples of compiling on several systems are given below. Users may need to customize as needed. If you wish to share more examples or to request or contribute in support for other compilers, please write to `dmitry@sdsc.edu <mailto:dmitry%40sdsc%2eedu>`_. If you give us an account on your system we will work with you to customize the installation.

.. csv-table::
        :header: "Argument", "Notes", "Description", "Example"
        :widths: auto
        :escape: '

        "--prefix=PREFIX", "Mandatory for users without access to /usr/local", "This argument will install p3dfft to PREFIX when you run make install. By default, configure will install to /usr/local", "--prefix=$HOME/local/"
        "--enable-gnu, --enable-ibm, --enable-intel, --enable-pgi, --enable-cray", "Mandatory", "These arguments will prepare p3dfft to be built by a specific compiler. You must only choose one option.", "--enable-pgi"
        "--enable-fftw, --enable-essl", "Mandatory", "These arguments will prepare p3dfft to be used with either the FFTW or ESSL library. You must only choose one option.", "--enable-fftw"
        "--with-fftw=FFTWLOCATION", "Mandatory if --enable-fftw is used", "This argument specifies the path location for the FFTW library; it is mandatory if you are planning to use p3dfft with the FFTW library.", "--enable-fftw --with-fftw=$FFTW_HOME"
        "--enable-openmp", "Mandatory if using multithreaded version", "This argument adds the appropriate compiler flags to enable OpenMP", "--enable-openmp"
        "--enable-openmpi", "Optional", "This argument uses the OpenMPI implementation of MPI", "--enable-openmpi"
        "--enable-oned", "Optional", "This argument is for 1D decomposition. The default is 2D decomposition but can be made to 1D by setting up a grid 1xn when running the code.", "--enable-oned"
        "--enable-estimate", "Optional, use only with --enable-fftw", "If this argument is passed, the FFTW library will not use run-time tuning to select the fastest algorithm for computing FFTs.", "--enable-estimate"
        "--enable-measure", "Optional, enabled by default, use only with --enable-fftw", "For search-once-for-the-fast algorithm (takes more time on p3dfft_setup()).", "--enable-measure"
        "--enable-patient", "Optional, use only with --enable-fftw", "For search-once-for-the-fastest-algorithm (takes much more time on p3dfft_setup()).", "--enable-patient"
        "--enable-dimsc", "Optional", "To assign processor rows and columns in the Cartesian processor grid according to C convention. The default is Fortran convention which is recommended. This option does not affect the order of storage of arrays in memory.", "--enable-dimsc"
        "--enable-useeven", "Optional, recommended for Cray XT", "This argument is for using MPI_Alltoall instead of MPI_Alltotallv. This will pad the send buffers with zeros to make them of equal size; not needed on most architecture but may lead to better results on Cray XT.", "--enable-useeven"
        "--enable-stride1", "Optional, recommended", "To enable stride-1 data structures on output (this may in some cases give some advantage in performance). You can define loop blocking factors NLBX and NBLY to experiment, otherwise they are set to default values.", "--enable-stride1"
        "--enable-nblx", "Optional", "To define loop blocking factor NBL_X", "--enable-nblx=32"
        "--enable-nbly1", "Optional", "To define loop blocking factor NBL_Y1", "--enable-nbly1=32"
        "--enable-nbly2", "Optional", "To define loop blocking factor NBL_Y2", "--enable-nbly2=32"
        "--enable-nblz", "Optional", "To define loop blocking factor NBL_Z", "--enable-nblz=32"
        "--enable-single", "Optional", "This argument will compile p3dfft in single-precision. By default, configure will setup p3dfft to be compiled in double-precision.", "--enable-single"
        "FC=<Fortran compiler>", "Strongly recommended", "Fortran compiler", "FC=mpif90"
        "FCFLAGS='"<Fortran compiler flags>'"", "Optional, recommended", "Fortran compiler flags", "FCFLAGS='"-O3'""
        "CC=<C compiler>", "Strongly Recommended", "C compiler", "CC=mpicc"
        "CFLAGS='"<C compiler flags>"", "Optional, recommended", "C compiler flags", "CFLAGS='"-O3'""
        "LDFLAGS='"<linker flags>"", "Optional", "Linker flags", ""

.. raw:: html

        <style>

        .tab {
                overflow: hidden;
                border: 1px solid #e1e4e5;
                background-color: #f1f1f1;
        }

        .tab button {
                background-color: inherit;
                float: left;
                border: none;
                outline: none;
                cursor: pointer;
                padding: 14px 16px;
                transition: 0.3s;
                font-size: 15px;
        }

        .tab button:hover {
                background-color: #ddd;
        }

        .tab button.active {
                background-color: #ccc;
        }

        .tabcontent {
                display: none;
                border-top: none;
        }
        </style>

Compiling on Comet (XSEDE/SDSC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Choose a MPI.

.. raw:: html

        <div class="tab">
                <button class="tablinks" onclick="openTab(event, 'Comet_Intel_MPI')">Intel MPI</button>
                <button class="tablinks" onclick="openTab(event, 'Comet_MVAPICH2')" id="defaultOpen">MVAPICH2</button>
                <button class="tablinks" onclick="openTab(event, 'Comet_Open_MPI')">Open MPI</button>
        </div>

        <div id="Comet_Intel_MPI" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw", "./configure --enable-intel --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc"
        "GNU", "gnu, fftw", "./configure --enable-gnu --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc"

.. raw:: html

        </div>

        <div id="Comet_MVAPICH2" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw", "./configure --enable-intel --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc CFLAGS=-lmpifort"
        "GNU", "gnu, fftw", "./configure --enable-gnu --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc CFLAGS=-lmpichf90"
        "PGI", "pgi, fftw", "./configure --enable-pgi --enable-fftw --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc"

.. raw:: html

        </div>

        <div id="Comet_Open_MPI" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw", "./configure --enable-intel --enable-fftw --enable-openmpi --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc"
        "GNU", "gnu, fftw", "./configure --enable-gnu --enable-fftw --enable-openmpi --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc CFLAGS=-lm"
        "PGI", "pgi, fftw", "./configure --enable-pgi --enable-fftw --enable-openmpi --with-fftw=$FFTWHOME FC=mpif90 CC=mpicc"

.. raw:: html

        </div><p></p>


Compiling on Stampede2 (XSEDE/TACC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Choose a MPI.

.. raw:: html

        <div class="tab">
                <button class="tablinks" onclick="openTab(event, 'Stampede2_Intel_MPI')" id="defaultOpen">Intel MPI</button>
                <button class="tablinks" onclick="openTab(event, 'Stampede2_MVAPICH2')">MVAPICH2</button>
        </div>

        <div id="Stampede2_Intel_MPI" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw3", "./configure --enable-intel --enable-fftw --with-fftw=$TACC_FFTW3_DIR FC=mpif90 CC=mpicc"

.. raw:: html

        </div>

        <div id="Stampede2_MVAPICH2" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel", "./configure --enable-intel --enable-fftw --with-fftw=/PATH/TO/FFTW/LIBRARY FC=mpif90 CC=mpicc CFLAGS=-lmpifort"

.. note::
        Stampede2's FFTW3 module is not compatible with its MVAPICH2 module. Users must install their own FFTW library.


.. raw:: html

        </div><p></p>

Compiling on Bridges (PSC)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Choose a MPI.

.. raw:: html

        <div class="tab">
                <button class="tablinks" onclick="openTab(event, 'Bridges_Intel_MPI')" id="defaultOpen">Intel MPI</button>
                <button class="tablinks" onclick="openTab(event, 'Bridges_MVAPICH2')">MVAPICH2</button>
                <button class="tablinks" onclick="openTab(event, 'Bridges_Open_MPI')">Open MPI</button>
        </div>

        <div id="Bridges_Intel_MPI" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw3", "./configure --enable-intel --enable-fftw --with-fftw=$FFTW3_LIB/.. FC=mpiifort CC=mpicc CFLAGS=-lm"

.. raw:: html

        </div>

        <div id="Bridges_MVAPICH2" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw3", "./configure --enable-intel --enable-fftw --with-fftw=$FFTW3_LIB/.. FC=mpif90 CC=mpicc CFLAGS=-lmpifort"
        "GNU", "fftw3", "./configure --enable-gnu --enable-fftw --with-fftw=$FFTW3_LIB/.. FC=mpif90 CC=mpicc CFLAGS=-lm"

.. raw:: html

        </div>

        <div id="Bridges_Open_MPI" class="tabcontent">

.. csv-table::
        :header: "Compiler", "Modules", "Arguments"
        :widths: auto

        "Intel", "intel, fftw3", "./configure --enable-intel --enable-fftw --enable-openmpi --with-fftw=$FFTW3_LIB/.. FC=mpif90 CC=mpicc"
        "GNU", "fftw3", "./configure --enable-gnu --enable-fftw --enable-openmpi --with-fftw=$FFTW3_LIB/.. FC=mpif90 CC=mpicc CFLAGS=-lm"
        "PGI", "pgi, fftw3", "./configure --enable-pgi --enable-fftw --enable-openmpi --with-fftw=$FFTW3_LIB/.. FC=mpif90 CC=mpicc"

.. raw:: html

        </div><p></p>

        <script>
        function openTab(evt, platform_mpi) {
                var i, tabcontent, tablinks;
                tabcontent = document.getElementsByClassName("tabcontent");
                for (i = 0; i < tabcontent.length; i++) {
                        tabcontent[i].style.display = "none";
                }
                tablinks = document.getElementsByClassName("tablinks");
                for (i = 0; i < tablinks.length; i++) {
                        tablinks[i].className = tablinks[i].className.replace(" active", "");
                }
                document.getElementById(platform_mpi).style.display = "block";
                evt.currentTarget.className += " active";
        }
        </script>


Compiling on Mira/Cetus/Vesta (ALCF)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. csv-table::
        :header: "Compiler", "Arguments"
        :widths: auto

        "IBM XL", "./configure --enable-ibm --enable-essl --with-essl=/soft/libraries/essl/current FC=mpixlf90_r CC=mpixlc_r"
        "GNU", "./configure --enable-gnu --enable-fftw --with-fftw=/soft/libraries/alcf/current/{xl,gcc}/FFTW3 FC=mpif90 CC=mpicc"
