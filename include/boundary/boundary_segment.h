#ifndef MPM_BOUNDARY_SEGMENT_H_
#define MPM_BOUNDARY_SEGMENT_H_

#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <vector>

#include "Eigen/Dense"

#include "cell.h"
#include "geometry.h"
#include "logger.h"
#include "particle.h"
#include "particle_base.h"

namespace mpm {

//! Boundar Segment class
//! \brief Brief class that stores the information about boundary segments
//! \tparam Tdim Dimension
template <unsigned Tdim>
class BoundarySegment {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Define DOFs
  static const unsigned Tdof = (Tdim == 1) ? 1 : 3 * (Tdim - 1);

  //! Constructor with id, pointers to boundary points
  //! \param[in] id Global boundary segment id
  //! \param[in] pointptrs Pointers to boundary point
  BoundarySegment(
      mpm::Index id,
      const std::vector<std::shared_ptr<ParticleBase<Tdim>>>& pointptrs);

  //! Default destructor
  ~BoundarySegment() = default;

  //! Delete copy constructor
  BoundarySegment(const BoundarySegment<Tdim>&) = delete;

  //! Delete assignement operator
  BoundarySegment& operator=(const BoundarySegment<Tdim>&) = delete;

  //! Return id of the cell
  mpm::Index id() const { return id_; }

  //! Return point coordinates
  VectorDim coordinates(const unsigned index) const {
    return points_.at(index)->coordinates();
  }

 private:
  // id
  mpm::Index id_{std::numeric_limits<Index>::max()};
  // Points which defines the boundary segment
  std::vector<std::shared_ptr<ParticleBase<Tdim>>> points_;
  // cells
  std::set<mpm::Index> cells_;
  //! Logger
  std::unique_ptr<spdlog::logger> console_;

};  // BoundarySegment class
}  // namespace mpm

#include "boundary_segment.tcc"

#endif  // MPM_BOUNDARY_SEGMENT_H_
