#ifndef TRILINOSINTERFACE_H
#define TRILINOSINTERFACE_H

#include <mpi.h>
#include <iostream>

// EPOCH3D errorcode
extern "C" int c_err_not_implemented;

// Fortran/C/C++ interface routines.
// Called from Fortran code

extern "C" {

  void init_trilinos_();
  void end_trilinos_();

}

#endif
