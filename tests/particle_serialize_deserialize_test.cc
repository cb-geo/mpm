#include <chrono>
#include <limits>

#include "catch.hpp"

#include "cell.h"
#include "data_types.h"
#include "element.h"
#include "function_base.h"
#include "hdf5_particle.h"
#include "hexahedron_element.h"
#include "linear_function.h"
#include "material.h"
#include "node.h"
#include "particle.h"
#include "quadrilateral_element.h"

//! \brief Check particle class for serialization and deserialization
TEST_CASE("Particle is checked for serialization and deserialization",
          "[particle][3D][serialize]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 3;
  // Number of phases
  const unsigned Nphases = 1;
  // Phase
  const unsigned phase = 0;

  // Check initialise particle from HDF5 file
  SECTION("Check initialise particle HDF5") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    // Coordinates
    Eigen::Matrix<double, Dim, 1> pcoords;
    pcoords.setZero();

    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::Particle<Dim>>(id, pcoords);

    mpm::HDF5Particle h5_particle;
    h5_particle.id = 13;
    h5_particle.mass = 501.5;

    Eigen::Vector3d coords;
    coords << 1., 2., 0.;
    h5_particle.coord_x = coords[0];
    h5_particle.coord_y = coords[1];
    h5_particle.coord_z = coords[2];

    Eigen::Vector3d displacement;
    displacement << 0.01, 0.02, 0.0;
    h5_particle.displacement_x = displacement[0];
    h5_particle.displacement_y = displacement[1];
    h5_particle.displacement_z = displacement[2];

    Eigen::Vector3d lsize;
    lsize << 0.25, 0.5, 0.;
    h5_particle.nsize_x = lsize[0];
    h5_particle.nsize_y = lsize[1];
    h5_particle.nsize_z = lsize[2];

    Eigen::Vector3d velocity;
    velocity << 1.5, 2.5, 0.0;
    h5_particle.velocity_x = velocity[0];
    h5_particle.velocity_y = velocity[1];
    h5_particle.velocity_z = velocity[2];

    Eigen::Matrix<double, 6, 1> stress;
    stress << 11.5, -12.5, 13.5, 14.5, -15.5, 16.5;
    h5_particle.stress_xx = stress[0];
    h5_particle.stress_yy = stress[1];
    h5_particle.stress_zz = stress[2];
    h5_particle.tau_xy = stress[3];
    h5_particle.tau_yz = stress[4];
    h5_particle.tau_xz = stress[5];

    Eigen::Matrix<double, 6, 1> strain;
    strain << 0.115, -0.125, 0.135, 0.145, -0.155, 0.165;
    h5_particle.strain_xx = strain[0];
    h5_particle.strain_yy = strain[1];
    h5_particle.strain_zz = strain[2];
    h5_particle.gamma_xy = strain[3];
    h5_particle.gamma_yz = strain[4];
    h5_particle.gamma_xz = strain[5];

    h5_particle.epsilon_v = strain.head(Dim).sum();

    h5_particle.status = true;

    h5_particle.cell_id = 1;

    h5_particle.volume = 2.;

    h5_particle.material_id = 1;

    // Reinitialise particle from HDF5 data
    REQUIRE(particle->initialise_particle(h5_particle) == true);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    unsigned mid = 1;

    auto material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(mid), jmaterial);
    std::vector<std::shared_ptr<mpm::Material<Dim>>> materials;
    materials.emplace_back(material);

    // Serialize particle
    auto buffer = particle->serialize();
    REQUIRE(buffer.size() > 0);

    // Deserialize particle
    std::shared_ptr<mpm::ParticleBase<Dim>> rparticle =
        std::make_shared<mpm::Particle<Dim>>(id, pcoords);

    REQUIRE_NOTHROW(rparticle->deserialize(buffer, materials));

    // Check particle id
    REQUIRE(particle->id() == particle->id());
    // Check particle mass
    REQUIRE(particle->mass() == rparticle->mass());
    // Check particle volume
    REQUIRE(particle->volume() == rparticle->volume());
    // Check particle status
    REQUIRE(particle->status() == rparticle->status());

    // Check for coordinates
    auto coordinates = rparticle->coordinates();
    REQUIRE(coordinates.size() == Dim);
    for (unsigned i = 0; i < coordinates.size(); ++i)
      REQUIRE(coordinates(i) == Approx(coords(i)).epsilon(Tolerance));

    // Check for displacement
    auto pdisplacement = rparticle->displacement();
    REQUIRE(pdisplacement.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pdisplacement(i) == Approx(displacement(i)).epsilon(Tolerance));

    // Check for size
    auto size = rparticle->natural_size();
    REQUIRE(size.size() == Dim);
    for (unsigned i = 0; i < size.size(); ++i)
      REQUIRE(size(i) == Approx(lsize(i)).epsilon(Tolerance));

    // Check velocity
    auto pvelocity = rparticle->velocity();
    REQUIRE(pvelocity.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pvelocity(i) == Approx(velocity(i)).epsilon(Tolerance));

    // Check stress
    auto pstress = rparticle->stress();
    REQUIRE(pstress.size() == stress.size());
    for (unsigned i = 0; i < stress.size(); ++i)
      REQUIRE(pstress(i) == Approx(stress(i)).epsilon(Tolerance));

    // Check strain
    auto pstrain = rparticle->strain();
    REQUIRE(pstrain.size() == strain.size());
    for (unsigned i = 0; i < strain.size(); ++i)
      REQUIRE(pstrain(i) == Approx(strain(i)).epsilon(Tolerance));

    // Check particle volumetric strain centroid
    REQUIRE(particle->volumetric_strain_centroid() ==
            rparticle->volumetric_strain_centroid());

    // Check cell id
    REQUIRE(particle->cell_id() == rparticle->cell_id());

    // Check material id
    REQUIRE(particle->material_id() == rparticle->material_id());

    SECTION("Performance benchmarks") {
      // Number of iterations
      unsigned niterations = 1000;

      // Serialization benchmarks
      auto serialize_start = std::chrono::steady_clock::now();
      for (unsigned i = 0; i < niterations; ++i) {
        // Serialize particle
        auto buffer = particle->serialize();
        // Deserialize particle
        std::shared_ptr<mpm::ParticleBase<Dim>> rparticle =
            std::make_shared<mpm::Particle<Dim>>(id, pcoords);

        REQUIRE_NOTHROW(rparticle->deserialize(buffer, materials));
      }
      auto serialize_end = std::chrono::steady_clock::now();
    }
  }
}
