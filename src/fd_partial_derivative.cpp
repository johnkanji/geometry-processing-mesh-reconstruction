#include "fd_partial_derivative.h"
#include <iostream>

void fd_partial_derivative(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const int dir,
  Eigen::SparseMatrix<double> & D)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  using T = Eigen::Triplet<double>;
  std::vector<T> tripletList;
  tripletList.reserve(nx * ny * nz * 2);

  int ndx = dir == 0 ? nx - 1 : nx;
  int ndy = dir == 1 ? ny - 1 : ny;
  int ndz = dir == 2 ? nz - 1 : nz;

  for (int i = 0; i < ndx; i++) {
    for (int j = 0; j < ndy; j++) {
      for (int k = 0; k < ndz; k++) {
        int i_offset = dir == 0 ? i + 1 : i;
        int j_offset = dir == 1 ? j + 1 : j;
        int k_offset = dir == 2 ? k + 1 : k;
        tripletList.push_back(T(i + ndx * j + ndx*ndy * k, i + nx * j + nx*ny * k, -1./h));
        tripletList.push_back(T(i + ndx * j + ndx*ndy * k, i_offset + nx * j_offset + nx*ny * k_offset, 1./h));
      }
    }
  }
  D.resize(ndx * ndy * ndz, nx * ny * nz);
  D.setFromTriplets(tripletList.begin(), tripletList.end());
}
