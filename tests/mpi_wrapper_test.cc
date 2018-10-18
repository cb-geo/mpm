#include <numeric>
#include <vector>

#include "Eigen/Dense"
#include "catch.hpp"
// MPI
#ifdef USE_MPI
#include "mpi.h"
#endif

#include "mpi_wrapper.h"
#include "node.h"

//! \brief Check MPI wrapper
TEST_CASE("MPI wrapper functions are checked", "[mpi]") {
  const double Tolerance = 1.E-7;

  //! Check chunked scalars
  SECTION("Check chunked scalars") {
    // Get all particle ids
    std::vector<mpm::Index> all_particles_ids(10);
    std::iota(all_particles_ids.begin(), all_particles_ids.end(), 0);

    // Get local particles ids chunks
    std::vector<mpm::Index> particles_ids;
    REQUIRE(particles_ids.size() == 0);
    mpm::chunk_scalar_quantities(all_particles_ids, particles_ids);
    REQUIRE(particles_ids.size() == all_particles_ids.size());

    // Check values
    for (unsigned i = 0; i < all_particles_ids.size(); ++i)
      REQUIRE(particles_ids.at(i) == all_particles_ids.at(i));
  }

  //! Check chunked vectors
  SECTION("Check chunked vectors") {
    // Get all particle ids
    std::vector<Eigen::Matrix<double, 3, 1>> all_particles;
    for (unsigned i = 0; i < 10; ++i)
      all_particles.emplace_back(Eigen::Matrix<double, 3, 1>::Constant(i));

    // Get local particles ids chunks
    std::vector<Eigen::Matrix<double, 3, 1>> particles;
    REQUIRE(particles.size() == 0);
    mpm::chunk_vector_quantities(all_particles, particles);
    REQUIRE(particles.size() == all_particles.size());

    // Check values
    for (unsigned i = 0; i < all_particles.size(); ++i)
      for (unsigned j = 0; j < all_particles.at(i).size(); ++j)
        REQUIRE(particles.at(i)[j] ==
                Approx(all_particles.at(i)[j]).epsilon(Tolerance));
  }
}
