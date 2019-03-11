/*struct Grid_struct_fort {
  int taskid;
  int numtasks;
  int nd;  //number of dimensions the volume is split over
  int gdims[3];  //Global dimensions
  int mem_order[3];  //Memory ordering inside the data volume
  int ldims[3];  //Local dimensions on THIS processor
  
  int pgrid[3];  //Processor grid
    int proc_order[3];   //Ordering of tasks in processor grid, e.g. (1,2,3) : first dimension - adjacent tasks,then second, then third dimension
  int glob_start[3]; // Starting coords of this cube in the global grid
  int mpi_comm_glob; // Global MPi communicator we are starting from
} ;

typedef struct Grid_struct_fort Grid_fort;
*/

  void p3dfft_init_3Dtype_f(int *,int[3]); //,char *);
  void p3dfft_plan_1Dtrans_f(int *,int *,int *,int *,int *,int *);
  void p3dfft_plan_3Dtrans_f(int *,int *,int *,Type3D *,int *);
void p3dfft_init_grid_f(int *,int *,int *,int *,int *,int *,int *,int *, int *);
void p3dfft_exec_1Dtrans_double_f(int *,double *,double *);
void p3dfft_exec_1Dtrans_single_f(int *,float *,float *);
void p3dfft_exec_3Dtrans_double_f(Plan3D *,double *,double *);
void p3dfft_exec_3Dtrans_single_f(Plan3D *,float *,float *);
void p3dfft_exec_3Dderiv_double_f(Plan3D *,double *,double *,int *);
void p3dfft_exec_3Dderiv_single_f(Plan3D *,float *,float *,int *);
void p3dfft_compute_deriv_single_f(float *,float *,int *,int *);
void p3dfft_compute_deriv_double_f(double *,double *,int *,int *);
