#include "poisson_surface_reconstruction.h"
#include "fd_interpolate.h"
#include "fd_grad.h"
#include <igl/copyleft/marching_cubes.h>
#include <iostream>
#include <algorithm>

#define print_size(s, M) std::cout << s << ": " << M.rows() << "x" << M.cols() << std::endl;

void poisson_surface_reconstruction(
    const Eigen::MatrixXd & P,
    const Eigen::MatrixXd & N,
    Eigen::MatrixXd & V,
    Eigen::MatrixXi & F)
{
  ////////////////////////////////////////////////////////////////////////////
  // Construct FD grid, CONGRATULATIONS! You get this for free!
  ////////////////////////////////////////////////////////////////////////////
  // number of input points
  const int n = P.rows();
  // Grid dimensions
  int nx, ny, nz;
  // Maximum extent (side length of bounding box) of points
  double max_extent =
    (P.colwise().maxCoeff()-P.colwise().minCoeff()).maxCoeff();
  // padding: number of cells beyond bounding box of input points
  const double pad = 8;
  // choose grid spacing (h) so that shortest side gets 30+2*pad samples
  double h  = max_extent/double(30+2*pad);
  // Place bottom-left-front corner of grid at minimum of points minus padding
  Eigen::RowVector3d corner = P.colwise().minCoeff().array()-pad*h;
  // Grid dimensions should be at least 3 
  nx = std::max((P.col(0).maxCoeff()-P.col(0).minCoeff()+(2.*pad)*h)/h,3.);
  ny = std::max((P.col(1).maxCoeff()-P.col(1).minCoeff()+(2.*pad)*h)/h,3.);
  nz = std::max((P.col(2).maxCoeff()-P.col(2).minCoeff()+(2.*pad)*h)/h,3.);
  // Compute positions of grid nodes
  Eigen::MatrixXd x(nx*ny*nz, 3);
  for(int i = 0; i < nx; i++) 
  {
    for(int j = 0; j < ny; j++)
    {
      for(int k = 0; k < nz; k++)
      {
         // Convert subscript to index
         const auto ind = i + nx*(j + k * ny);
         x.row(ind) = corner + h*Eigen::RowVector3d(i,j,k);
      }
    }
  }
  Eigen::VectorXd g = Eigen::VectorXd::Zero(nx*ny*nz);

  ////////////////////////////////////////////////////////////////////////////
  // Add your code here
  ////////////////////////////////////////////////////////////////////////////
  Eigen::SparseMatrix<double> W, Wx, Wy, Wz, G;
  Eigen::VectorXd v, vx, vy, vz;
  Eigen::RowVector3d dx(h/2., 0., 0.);
  Eigen::RowVector3d dy(0., h/2., 0.);
  Eigen::RowVector3d dz(0., 0., h/2.);
  fd_interpolate(nx - 1, ny, nz, h, corner + dx, P.rowwise() + dx, Wx);
  fd_interpolate(nx, ny - 1, nz, h, corner + dy, P.rowwise() + dy, Wy);
  fd_interpolate(nx, ny, nz - 1, h, corner + dz, P.rowwise() + dz, Wz);
  fd_interpolate(nx, ny, nz, h, corner, P, W);
  vx = Wx.transpose() * N.col(0);
  vy = Wy.transpose() * N.col(1);
  vz = Wz.transpose() * N.col(2);
  v.resize(vx.size() + vy.size() + vz.size());
  v << vx, vy, vz;

  fd_grad(nx, ny, nz, h, G);
  auto A = G.transpose() * G;
  auto b = G.transpose() * v;

  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
  solver.compute(A);
  g = solver.solve(b);
  std::cout << "#iterations:     " << solver.iterations() << std::endl;
  std::cout << "estimated error: " << solver.error()      << std::endl;

  double sigma = 1./W.rows() * Eigen::RowVectorXd::Ones(g.rows()) * W * g;
  g = g - Eigen::VectorXd::Ones(g.size()) * sigma;


  ////////////////////////////////////////////////////////////////////////////
  // Run black box algorithm to compute mesh from implicit function: this
  // function always extracts g=0, so "pre-shift" your g values by -sigma
  ////////////////////////////////////////////////////////////////////////////
  igl::copyleft::marching_cubes(g, x, nx, ny, nz, V, F);
}
