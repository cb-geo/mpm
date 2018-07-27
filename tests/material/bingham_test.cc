#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "material/material.h"

//! \brief Check Bingham class
TEST_CASE("Bingham is checked", "[material][bingham]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  //! Check for id = 0
  SECTION("Bingham id is zero") {
    unsigned id = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "Bingham", std::move(id));
    REQUIRE(material->id() == 0);
  }

  SECTION("Bingham id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "Bingham", std::move(id));
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Read material properties
  SECTION("Bingham check properties") {
    unsigned id = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "Bingham", std::move(id));
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
    
    // Calculate modulus values
    const double K = 8333333.333333333;
    const double G = 3846153.846153846;

  }

  SECTION("Bingham check stresses") {
    unsigned id = 0;
    auto material = Factory<mpm::Material, unsigned>::instance()->create(
        "Bingham", std::move(id));
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

    // Initialise stress
    mpm::Material::Vector6d stress;
    stress.setZero();
    REQUIRE(stress(0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.).epsilon(Tolerance));

    // Initialise dstrain
    mpm::Material::Vector6d dstrain;
    dstrain.setZero();
    dstrain(0) = 0.0010000;
    dstrain(1) = 0.0005000;
    dstrain(2) = 0.0005000;
    dstrain(3) = 0.0000000;
    dstrain(4) = 0.0000000;
    dstrain(5) = 0.0000000;

    // Initialise strain rate
    mpm::Material::Vector6d strain_rate;
    strain_rate.setZero();
    strain_rate(0) = 0.0000100;
    strain_rate(1) = 0.0000050;
    strain_rate(2) = 0.0000050;
    strain_rate(3) = 0.0000000;
    strain_rate(4) = 0.0000000;
    strain_rate(5) = 0.0000000;

    // Compute updated stress
    material->compute_stress(stress, dstrain, strain_rate);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise dstrain
    dstrain(0) = 0.0010000;
    dstrain(1) = 0.0005000;
    dstrain(2) = 0.0005000;
    dstrain(3) = 0.0000100;
    dstrain(4) = 0.0000200;
    dstrain(5) = 0.0000300;
    
    // Initialise strain rate
    strain_rate(0) = 0.0000100;
    strain_rate(1) = 0.0000050;
    strain_rate(2) = 0.0000050;
    strain_rate(3) = 0.0000001;
    strain_rate(4) = 0.0000002;
    strain_rate(5) = 0.0000003;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    material->compute_stress(stress, dstrain, strain_rate);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(7.69230769230769e+01).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(1.15384615384615e+02).epsilon(Tolerance));
  }
}
