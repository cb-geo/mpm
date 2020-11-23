#include <limits>

#include "catch.hpp"

#include "data_types.h"
#include "material.h"
#include "particle.h"
#include "particle_twophase.h"
#include "pod_particle_twophase.h"

//! \brief Check particle class for serialization and deserialization
TEST_CASE("Twophase particle is checked for serialization and deserialization",
          "[particle][3D][serialize][2Phase]") {
  // Dimension
  const unsigned Dim = 3;
  // Dimension
  const unsigned Dof = 3;
  // Number of phases
  const unsigned Nphases = 2;
  // Phase
  const unsigned phase = 0;

  // Check initialise particle from POD file
  SECTION("Check initialise particle POD") {
    mpm::Index id = 0;
    const double Tolerance = 1.E-7;
    // Coordinates
    Eigen::Matrix<double, Dim, 1> pcoords;
    pcoords.setZero();

    std::shared_ptr<mpm::ParticleBase<Dim>> particle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, pcoords);

    // Assign material
    unsigned solid_mid = 1;
    unsigned liquid_mid = 2;
    // Initialise material
    Json jsolid_material;
    Json jliquid_material;
    jsolid_material["density"] = 1000.;
    jsolid_material["youngs_modulus"] = 1.0E+7;
    jsolid_material["poisson_ratio"] = 0.3;
    jsolid_material["porosity"] = 0.3;
    jsolid_material["k_x"] = 0.001;
    jsolid_material["k_y"] = 0.001;
    jsolid_material["k_z"] = 0.001;
    jliquid_material["density"] = 1000.;
    jliquid_material["bulk_modulus"] = 2.0E9;
    jliquid_material["dynamic_viscosity"] = 8.90E-4;

    auto solid_material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "LinearElastic3D", std::move(solid_mid), jsolid_material);
    auto liquid_material =
        Factory<mpm::Material<Dim>, unsigned, const Json&>::instance()->create(
            "Newtonian3D", std::move(liquid_mid), jliquid_material);
    std::vector<std::shared_ptr<mpm::Material<Dim>>> materials;
    materials.emplace_back(solid_material);
    materials.emplace_back(liquid_material);

    mpm::PODParticleTwoPhase h5_particle;
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

    h5_particle.nstate_vars = 0;

    for (unsigned i = 0; i < h5_particle.nstate_vars; ++i)
      h5_particle.svars[i] = 0.;

    h5_particle.liquid_mass = 100.1;

    Eigen::Vector3d liquid_velocity;
    liquid_velocity << 5.5, 2.1, 4.2;
    h5_particle.liquid_velocity_x = liquid_velocity[0];
    h5_particle.liquid_velocity_y = liquid_velocity[1];
    h5_particle.liquid_velocity_z = liquid_velocity[2];

    h5_particle.porosity = 0.33;

    h5_particle.liquid_saturation = 1.;

    h5_particle.liquid_material_id = 2;

    h5_particle.nliquid_state_vars = 1;

    for (unsigned i = 0; i < h5_particle.nliquid_state_vars; ++i)
      h5_particle.liquid_svars[i] = 0.;

    // Reinitialise particle from POD data
    REQUIRE(particle->initialise_particle(h5_particle, materials) == true);

    // Serialize particle
    auto buffer = particle->serialize();
    REQUIRE(buffer.size() > 0);

    // Deserialize particle
    std::shared_ptr<mpm::ParticleBase<Dim>> rparticle =
        std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, pcoords);

    REQUIRE_NOTHROW(rparticle->deserialize(buffer, materials));

    // Check particle id
    REQUIRE(particle->id() == particle->id());
    // Check particle mass
    REQUIRE(particle->mass() == rparticle->mass());
    // Check particle volume
    REQUIRE(particle->volume() == rparticle->volume());
    // Check particle mass density
    REQUIRE(particle->mass_density() == rparticle->mass_density());
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

    // Check liquid mass
    REQUIRE(particle->liquid_mass() == rparticle->liquid_mass());

    // Check liquid velocity
    auto pliquid_velocity = rparticle->liquid_velocity();
    REQUIRE(pliquid_velocity.size() == Dim);
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(pliquid_velocity(i) ==
              Approx(liquid_velocity(i)).epsilon(Tolerance));

    // Check porosity
    REQUIRE(particle->porosity() == rparticle->porosity());

    // Check liquid material id
    REQUIRE(particle->material_id(mpm::ParticlePhase::Liquid) ==
            rparticle->material_id(mpm::ParticlePhase::Liquid));

    // Check state variables
    for (unsigned phase = 0; phase < materials.size(); phase++) {
      REQUIRE(particle->state_variables(phase).size() ==
              rparticle->state_variables(phase).size());
      auto state_variables = materials[phase]->state_variables();
      for (const auto& state_var : state_variables) {
        REQUIRE(particle->state_variable(state_var, phase) ==
                rparticle->state_variable(state_var, phase));
      }
    }

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
            std::make_shared<mpm::TwoPhaseParticle<Dim>>(id, pcoords);

        REQUIRE_NOTHROW(rparticle->deserialize(buffer, materials));
      }
      auto serialize_end = std::chrono::steady_clock::now();
    }
  }
}
