#include "TrilinosInterface.h"

extern "C" {

  void init_trilinos_(int& NumMyElements, int* MyGlobalElements, int& ForComm) {

    MPI_Comm Comm;
    // This function obtains a valid C handle to the Fortran MPI communicator
    Comm = MPI_Comm_f2c(ForComm);

    // Create the solver instance
    PICSolver = new JFNKSolver(NumMyElements, MyGlobalElements, Comm);
  }

  void end_trilinos_() {

    delete PICSolver;

    std::cout << "End trilinos not implemented yet!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, c_err_not_implemented);
  }

  void solve_gmres_(double* x, double *dir) {

    PICSolver->Solve(x, dir);

  }
}
