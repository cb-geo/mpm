#ifndef MPM_MULTIMATERIAL_H_
#define MPM_MULTIMATERIAL_H_

#include <vector>

#include "Eigen/Dense"
#include "Eigen/LU"

namespace mpm {
// \brief Multimaterial parameters on each node
// \tparam Tdim Dimension
// \tparam Tnmat Number of materials in the model
// \tparam Tnnodes Number of nodes in the mesh
template <int Tdim, int Tnmat, int Tnnodes>
struct multimaterial_properties {
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  Eigen::Matrix<double, Tnmat, Tnnodes> mass;
  Eigen::Matrix<VectorDim, Tnmat, Tnnodes> velocity;
  Eigen::Matrix<VectorDim, Tnmat, Tnnodes> momentum;
};

}  // namespace mpm

#endif  // MPM_MULTIMATERIAL_H_
