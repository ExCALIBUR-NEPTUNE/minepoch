#ifndef TRILINOSINTERFACE_H
#define TRILINOSINTERFACE_H

#include <mpi.h>
#include <iostream>

#include "JFNKSolver.h"

// EPOCH3D errorcode
extern "C" int c_err_not_implemented;

// Global variables

JFNKSolver *PICSolver;

// Fortran/C/C++ interface routines.
// Called from Fortran code

extern "C" {

  void init_trilinos_(int &, int *, int &, bool &);
  void end_trilinos_();

}

#endif
