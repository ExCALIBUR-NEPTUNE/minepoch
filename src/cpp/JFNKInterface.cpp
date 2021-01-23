#include "JFNKInterface.h"

//extern "C" {
  // External Fortran routine to calculate residual
  // void computef(double *,double *,int *);
//}

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

  std::cout << "JFNKInterface::computeF not implemented yet!" << std::endl;
  MPI_Abort(MPI_COMM_WORLD, c_err_not_implemented);
  // Compute rhs by calling the external Fortran routine
  // computef(xptr,fptr,&iflag);

  return true;
}
