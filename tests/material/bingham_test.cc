#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "material/material.h"

//! \brief Check Bingham class
TEST_CASE("Bingham is checked in 3D", "[material][bingham][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;

  //! Check for id = 0
  SECTION("Bingham id is zero") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);
  }

  SECTION("Bingham id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Read material properties
  SECTION("Bingham check properties") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    jmaterial["tau0"] = 771.8;
    jmaterial["mu"] = 0.0451;
    jmaterial["strain_cutoff"] = 0.2;

    // Check material status before assigning material property
    REQUIRE(material->status() == false);
    
    // Get material properties
    REQUIRE(material->property("density") ==
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    // Check for property that does not exist
    REQUIRE(material->property("noproperty") ==
            Approx(std::numeric_limits<double>::max()).epsilon(Tolerance));

    material->properties(jmaterial);

    // Check material status after assigning material property
    REQUIRE(material->status() == true);

    // Get material properties
    REQUIRE(material->property("density") == Approx(jmaterial["density"]).epsilon(Tolerance));
    
  }

  SECTION("Bingham check stresses") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "Bingham3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;
    jmaterial["tau0"] = 771.8;
    jmaterial["mu"] = 0.0451;
    jmaterial["strain_cutoff"] = 0.2;

    material->properties(jmaterial);

    // Check material status after assigning material property
    REQUIRE(material->property_handle() == true);

    // Initialise stress
    mpm::Material<Dim>::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise strain
    mpm::Material<Dim>::Vector6d strain;
    strain.setZero();
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0005000;
    strain(3) = 0.0000000;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Initialise strain rate
    mpm::Material<Dim>::Vector6d strain_rate;
    strain_rate.setZero();
    strain_rate(0) = 0.0000100;
    strain_rate(1) = 0.0000050;
    strain_rate(2) = 0.0000050;
    strain_rate(3) = 0.0000000;
    strain_rate(4) = 0.0000000;
    strain_rate(5) = 0.0000000;

    // Add particle
    mpm::Index pid = 0;
    Eigen::Matrix<double, Dim, 1> coords;
    coords << 0.75, 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

    // Reset stress
    stress.setZero();
    // Compute updated stress
    stress = material->compute_stress(stress, strain, particle.get());

    // Check stressees
    REQUIRE(stress(0) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(16666.66666667).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));



    // // Check stress with some shear strain_rate but does not reach threshold
    // // Initialise strain
    // strain(0) = 0.0010000;
    // strain(1) = 0.0005000;
    // strain(2) = 0.0005000;
    // strain(3) = 0.0000100;
    // strain(4) = 0.0000200;
    // strain(5) = 0.0000300;
    
    // // Initialise strain rate
    // strain_rate(0) = 0.0000100;
    // strain_rate(1) = 0.0000050;
    // strain_rate(2) = 0.0000050;
    // strain_rate(3) = 0.0000001;
    // strain_rate(4) = 0.0000002;
    // strain_rate(5) = 0.0000003;

    // // Reset stress
    // stress.setZero();

    // // Compute updated stress
    // stress = material->compute_stress(stress, strain, particle.get());

    // // Check stressees
    // REQUIRE(stress(0) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));



    // // Check stress with some shear strain_rate that does reach threshold
    // // Initialise strain
    // strain(0) = 0.0010000;
    // strain(1) = 0.0005000;
    // strain(2) = 0.0005000;
    // strain(3) = 0.0000100;
    // strain(4) = 0.0000200;
    // strain(5) = 0.0000300;
    
    // // Initialise strain rate
    // strain_rate(0) = 1.000000;
    // strain_rate(1) = 0.500000;
    // strain_rate(2) = 0.500000;
    // strain_rate(3) = 0.010000;
    // strain_rate(4) = 0.020000;
    // strain_rate(5) = 0.030000;

    // // Reset stress
    // stress.setZero();

    // // Compute updated stress
    // stress = material->compute_stress(stress, strain, particle.get());

    // // Check stressees
    // REQUIRE(stress(0) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(1) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(2) == Approx(16666.66666666667).epsilon(Tolerance));
    // REQUIRE(stress(3) == Approx(8.908724740775902).epsilon(Tolerance));
    // REQUIRE(stress(4) == Approx(17.81744948155181).epsilon(Tolerance));
    // REQUIRE(stress(5) == Approx(26.72617422232771).epsilon(Tolerance));

  }
}
