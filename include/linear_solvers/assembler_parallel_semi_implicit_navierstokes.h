#ifndef MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_
#define MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_eigen_semi_implicit_navierstokes.h"
#include "spdlog/spdlog.h"

#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerParallelSemiImplicitNavierStokes
    : public AssemblerEigenSemiImplicitNavierStokes<Tdim> {
 public:
  //! Constructor
  AssemblerParallelSemiImplicitNavierStokes();

  //! Create a pair between nodes and index in Matrix / Vector
  bool assign_global_node_indices(unsigned nactive_node,
                                  unsigned nglobal_active_node) override;

  unsigned global_active_dof() override { return global_active_dof_; };

  std::vector<int> rank_global_mapper() override {
    return rank_global_mapper_;
  };

 protected:
  //! Logger
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::console_;
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Global node indices
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::global_node_indices_;
  //! Number of total active_dof in all rank
  unsigned global_active_dof_;
  //! Rank to Global mapper
  std::vector<int> rank_global_mapper_;
};
}  // namespace mpm

#include "assembler_parallel_semi_implicit_navierstokes.tcc"
#endif  // MPM_ASSEMBLER_PARALLEL_SEMI_IMPLICIT_NAVIERSTOKES_H_
