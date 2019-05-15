#include "fd_grad.h"
#include "fd_partial_derivative.h"
#include "igl/cat.h"

void fd_grad(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  Eigen::SparseMatrix<double> & G)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  int rows = (nx - 1) * ny * nz + nx * (ny - 1) * nz + nx * ny * (nz - 1);
  Eigen::SparseMatrix<double> Dx, Dy, Dz, tmp;
  fd_partial_derivative(nx, ny, nz, h, 0, Dx);
  fd_partial_derivative(nx, ny, nz, h, 1, Dy);
  fd_partial_derivative(nx, ny, nz, h, 2, Dz);
  
  G.resize(rows, nx * ny * nz);
  igl::cat(1, Dx, Dy, tmp);
  igl::cat(1, tmp, Dz, G);
}
