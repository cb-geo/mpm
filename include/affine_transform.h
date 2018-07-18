#ifndef MPM_AFFINE_TRANSFORM_H_
#define MPM_AFFINE_TRANSFORM_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/LU"

namespace mpm {

// The algorithm to compute the affine approximation to the point on the
// unit cell does the following steps:
// <ul>
// <li> find the least square dim-dimensional plane approximating the cell
// vertices, i.e. we find an affine map A x_hat + b from the reference cell
// to the real space.
// <li> Solve the equation A x_hat + b = p for x_hat
// </ul>
//
// We look for a matrix X such that X * M = Y where M is
// a (Dim+1) x nvertices matrix (M[Dim][nvertices]) and
// Y is a dim x nvertices (Y[dim][nvertices]).
// The i-th column of M is unit_vertex[i] and the last row all 1's.
// The i-th column of Y is real_vertex[i]. If we split X=[A|b],
// the least square approx is A x_hat+b,
// Classically X = Y * (M^t (M M^t)^{-1}).
// Let K = M^t * (M M^t)^{-1} = [KA Kb] this can be precomputed,
// and that is exactly what we do. Finally A = Y*KA and b = Y*Kb.

//! \brief Affine transform real to unit cell
//! \tparam Tdim Dimension
//! \tparam Tnfunc Number of functions (vertices)
template <int Tdim, int Tnfunc>
struct TransformR2UAffine {
  static const Eigen::Matrix<double, Tnfunc, Tdim> KA;
  static const Eigen::Matrix<double, Tnfunc, 1> Kb;
};
}  // namespace mpm

#endif  // MPM_AFFINE_TRANSFORM_H_
