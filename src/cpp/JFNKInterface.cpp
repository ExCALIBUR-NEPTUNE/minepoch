#include "JFNKInterface.h"

extern "C" {
  //  External Fortran routine to calculate residual
  void jfnk_computef(double *,double *);
}

JFNKInterface::JFNKInterface(){
}

JFNKInterface::~JFNKInterface(){
}

bool JFNKInterface::computeF(const Epetra_Vector& x,Epetra_Vector& f
  ,FillType flag){

  // Extract a view for x and f
  double * fptr, * xptr;

  // Use to inform Fortran modules what needs to be calculated
  // TODO: Use named constants defined in Fortran
  int iflag;

  x.ExtractView(&xptr);
  f.ExtractView(&fptr);

  if(flag==NOX::Epetra::Interface::Required::Residual){
    iflag = 0;
  }
  else if(flag==NOX::Epetra::Interface::Required::FD_Res){
    iflag = 1;
  }
  else {
    // TODO. Check.
    iflag = 1;
  }

  // Compute rhs by calling the external Fortran routine. Should an iflag be passed here?
  jfnk_computef(xptr,fptr);

  return true;
}
