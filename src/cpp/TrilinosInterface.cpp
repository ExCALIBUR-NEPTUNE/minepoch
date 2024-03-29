#include "TrilinosInterface.h"

extern "C" {

  void init_trilinos_(int& NumMyElements, int* MyGlobalElements, int& ForComm, bool& verbose) {

    MPI_Comm Comm;
    // This function obtains a valid C handle to the Fortran MPI communicator
    Comm = MPI_Comm_f2c(ForComm);

    // Create the solver instance
    PICSolver = new JFNKSolver(NumMyElements, MyGlobalElements, Comm, verbose);
  }

  void end_trilinos_() {

    delete PICSolver;
  }

  void solve_gmres_(double* x, double *dir, double &eta) {

    PICSolver->Solve(x, dir, eta);

  }
}
