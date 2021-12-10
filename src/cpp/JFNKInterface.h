#ifndef JFNKINTERFACE_H
#define JFNKINTERFACE_H

#include <mpi.h>

#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra.H"

// EPOCH3D errorcode
extern "C" int c_err_not_implemented;

class JFNKInterface :public NOX::Epetra::Interface::Required {
public:
  JFNKInterface();
  ~JFNKInterface();
  bool computeF(const Epetra_Vector&, Epetra_Vector&, FillType);
};

#endif
