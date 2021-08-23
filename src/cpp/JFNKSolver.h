#ifndef JFNKSOLVER_H
#define JFNKSOLVER_H

#include <mpi.h>
#include <iostream>

#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MpiComm.h"
#include <NOX_Epetra_MatrixFree.H>

#include "JFNKInterface.h"

// EPOCH3D errorcode
extern "C" int c_err_not_implemented;

class JFNKSolver {
 private:

  Epetra_LinearProblem Problem;
  AztecOO *AztecSolver;

  // Epetra objects
  Epetra_MpiComm* Comm;
  Epetra_Map* Map;

  // Epetra vectors
  Epetra_Vector* rhs;
  Epetra_Vector* x;

  // Jacobian interface
  Teuchos::RCP<JFNKInterface> Interface;

  // Matrix free operator
  Teuchos::RCP<NOX::Epetra::MatrixFree> JacFree;

  int NumMyElements;

public:
  JFNKSolver(int, int *, MPI_Comm &, bool);
  ~JFNKSolver();
  void CreateJacobian();
  void Solve(double *, double *, double);
};

#endif
