#include "fd_interpolate.h"
#include "igl/floor.h"
#include <vector>

void fd_interpolate(
  const int nx,
  const int ny,
  const int nz,
  const double h,
  const Eigen::RowVector3d & corner,
  const Eigen::MatrixXd & P,
  Eigen::SparseMatrix<double> & W)
{
  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////

  Eigen::MatrixXd P_grid = (P.rowwise() - corner) / h;
  Eigen::MatrixXd P_grid_corners;
  igl::floor(P_grid, P_grid_corners);

  using T = Eigen::Triplet<double>;
  std::vector<T> tripletList;
  tripletList.reserve(P.rows() * 8);

  int i, j, k;
  double wx, wy, wz;
  for (int l = 0; l < P.rows(); l++) {
    i = P_grid_corners.row(l)(0);
    j = P_grid_corners.row(l)(1);
    k = P_grid_corners.row(l)(2);

    wx = P_grid.row(l)(0) - i;
    wy = P_grid.row(l)(1) - j;
    wz = P_grid.row(l)(2) - k;

    tripletList.push_back(T(l, i + nx * j + nx * ny * k,                   (1. - wx) * (1. - wy) * (1. - wz)));
    tripletList.push_back(T(l, (i + 1) + nx * j + nx * ny * k,             wx * (1. - wy) * (1. - wz)));
    tripletList.push_back(T(l, i + nx * (j + 1) + nx * ny * k,             (1. - wx) * wy * (1. - wz)));
    tripletList.push_back(T(l, (i + 1) + nx * (j + 1) + nx * ny * k,       wx * wy * (1. - wz)));
    tripletList.push_back(T(l, i + nx * j + nx * ny * (k + 1),             (1. - wx) * (1. - wy) * wz));
    tripletList.push_back(T(l, (i + 1) + nx * j + nx * ny * (k + 1),       wx * (1. - wy) * wz));
    tripletList.push_back(T(l, i + nx * (j + 1) + nx * ny * (k + 1),       (1. - wx) * wy * wz));
    tripletList.push_back(T(l, (i + 1) + nx * (j + 1) + nx * ny * (k + 1), wx * wy * wz));
  }
  W.resize(P.rows(), nx * ny * nz);
  W.setFromTriplets(tripletList.begin(), tripletList.end());
}
