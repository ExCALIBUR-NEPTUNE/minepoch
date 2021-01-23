#ifndef JFNKSOLVER_H
#define JFNKSOLVER_H

#include <mpi.h>
#include <iostream>

#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_MpiComm.h"

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

 public:
  JFNKSolver(int &, int *, MPI_Comm &);
  ~JFNKSolver();
  void CreateJacobian();
};

#endif
