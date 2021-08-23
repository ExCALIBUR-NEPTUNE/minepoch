#include "JFNKSolver.h"

JFNKSolver::JFNKSolver(int NumMyElements_, int * MyGlobalElements, MPI_Comm &OldComm) :
  NumMyElements(NumMyElements_)
{

  // Set-up the communicator and distributed objects
  // MPI_COMM_WORLD -is- an acceptable input here. See Epetra_MpiComm doxygen.
  this->Comm = new Epetra_MpiComm(OldComm);

  // Set-up map
  this->Map = new Epetra_Map(-1, NumMyElements, MyGlobalElements, 0, *Comm);

  // Create Epetra vectors
  rhs = new Epetra_Vector(*Map);
  x = new Epetra_Vector(*Map);

  // Set LHS and RHS of problem
  Problem.SetLHS(x);
  Problem.SetRHS(rhs);

  // Create Jacobian
  this->CreateJacobian();

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
  // AztecSolver->SetAztecOption(AZ_output,AZ_summary);
  AztecSolver->SetAztecOption(AZ_output, AZ_none);

  // Turn off built-in preconditioning method
  // Any preconditioning options should be set here.
  AztecSolver->SetAztecOption(AZ_precond,AZ_none);
}

JFNKSolver::~JFNKSolver() {
  delete AztecSolver;
  delete rhs;
  delete x;
}

void JFNKSolver::CreateJacobian() {

  // Create the Jacobian interface
  Interface = Teuchos::rcp( new JFNKInterface() );

  // Create a NOX::Epetra::Vector (for the FiniteDifference class contructor)
  Teuchos::RCP<Epetra_Vector> Guess = Teuchos::rcp(new Epetra_Vector(*Map));
  NOX::Epetra::Vector noxGuess(Guess, NOX::Epetra::Vector::CreateView);

  // Create top level parameter list. Extra output sublist.
  Teuchos::RCP<Teuchos::ParameterList> nlParamsPtr = Teuchos::rcp(new Teuchos::ParameterList);
  Teuchos::ParameterList & nlParams = *(nlParamsPtr.get());
  Teuchos::ParameterList & printParams = nlParams.sublist("Printing");
  // TODO Fix MyPID
  printParams.set("MyPID",0);
  printParams.set("Output Precision",3);
  printParams.set("Output Processor",0);

  // Final argument is boolean controlling use of perturbation type
  JacFree = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams,Interface,noxGuess,false));
  Problem.SetOperator(&(*JacFree));
  // Used in calculation perturbation. 1e-6 is default value
  JacFree->setLambda(1e-6);
}

void JFNKSolver::Solve(double x0[], double p[], double eta) {

  double norm;

  Epetra_Vector guess(Copy, *Map, x0);

  Interface->computeF(guess,*rhs,JFNKInterface::Residual);

  // TODO toggle this output
  // rhs->Norm2(&norm);
  // std::cout << "||b|| =" << norm << std::endl;

  JacFree->computeJacobian(guess,*JacFree);

  // Set the initial guess to zero
  x->Scale(0.);

  AztecSolver->Iterate(300,eta);

  // Copy solution vector to p
  for(int i=0; i<NumMyElements; i++){
    p[i] = x[0][i];
  }
}
