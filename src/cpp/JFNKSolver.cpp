#include "JFNKSolver.h"

JFNKSolver::JFNKSolver(MPI_Comm &OldComm) {

  // Set-up the communicator and distributed objects
  // MPI_COMM_WORLD -is- an acceptable input here. See Epetra_MpiComm doxygen.
  this->Comm = new Epetra_MpiComm(OldComm);

  // Set-up Linear Problem
  AztecSolver = new AztecOO(Problem);

  // GMRES convergence criteria
  // AZ_r0 is the default option here
  // Other options AZ_rhs, AZ_Anorm, AZ_noscaled, AZ_sol, AZ_weighted
  // See user guide (pg 16): https://prod.sandia.gov/techlib-noauth/access-control.cgi/2004/043796.pdf
  AztecSolver->SetAztecOption(AZ_conv,AZ_r0);

  // Select solver method. GMRES.
  // AZ_gmres_condnum would provide additional info, including condition number estimate
  AztecSolver->SetAztecOption(AZ_solver,AZ_gmres);

  // Maximum dimension of Krylov subspace.
  // GMRES restarts every AZ_kspace iterations
  AztecSolver->SetAztecOption(AZ_kspace,40);

  // Set output level.
  // Possible options: AZ_none, AZ_summary, AZ_warning, AZ_last, AZ_all
  AztecSolver->SetAztecOption(AZ_output,AZ_summary);

  // Turn off built-in preconditioning method
  // Any preconditioning options should be set here.
  AztecSolver->SetAztecOption(AZ_precond,AZ_none);

  // Create Jacobian
  this->CreateJacobian();
}

JFNKSolver::~JFNKSolver() {
  delete AztecSolver;
}

void JFNKSolver::CreateJacobian() {
  std::cout << "JFNKSolver::CreateJacobian not implemented yet!" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, c_err_not_implemented);
}
