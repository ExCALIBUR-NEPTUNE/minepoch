#ifndef JFNKSOLVER_H
#define JFNKSOLVER_H

#include <mpi.h>
#include <iostream>

#include "AztecOO.h"
#include "Epetra_LinearProblem.h"

// EPOCH3D errorcode
extern "C" int c_err_not_implemented;

class JFNKSolver {
 private:

  Epetra_LinearProblem Problem;
  AztecOO *AztecSolver;

 public:
  JFNKSolver();
  ~JFNKSolver();
  void CreateJacobian();
};

#endif
