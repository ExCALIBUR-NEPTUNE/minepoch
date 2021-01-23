#include "TrilinosInterface.h"

extern "C" {

  void init_trilinos_() {
    std::cout << "Init trilinos not implemented yet!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, c_err_not_implemented);
  }

  void end_trilinos_() {
    std::cout << "End trilinos not implemented yet!" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, c_err_not_implemented);
  }
}
