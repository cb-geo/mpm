#include <limits>

#include "Eigen/Dense"
#include "catch.hpp"
#include "json.hpp"

#include "cell.h"
#include "hex_shapefn.h"
#include "material/material.h"
#include "node.h"
#include "particle.h"
#include "quad_shapefn.h"
#include "shapefn.h"

//! Check linearelastic class in 2D
TEST_CASE("LinearElastic is checked in 2D", "[material][linear_elastic][2D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 2;
  
  //! Check for id = 0
  SECTION("LinearElastic id is zero") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic2D", std::move(id));
    REQUIRE(material->id() == 0);
  }

  SECTION("LinearElastic id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic2D", std::move(id));
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Read material properties
  SECTION("LinearElastic check stiffness matrix") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic2D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

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
    const double a1 = 13461538.461566667;
    const double a2 = 5769230.769166667;

    mpm::Material<Dim>::Matrix6x6 de = material->elastic_tensor();
    REQUIRE(de(0, 0) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(0, 1) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(0, 2) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(0, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(0, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(0, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(1, 0) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(1, 1) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(1, 2) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(1, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(1, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(1, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(2, 0) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(2, 1) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(2, 2) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(2, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(2, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(2, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(3, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 3) == Approx(G).epsilon(Tolerance));
    REQUIRE(de(3, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(4, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 4) == Approx(G).epsilon(Tolerance));
    REQUIRE(de(4, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(5, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 5) == Approx(G).epsilon(Tolerance));
  }

  SECTION("LinearElastic check stresses") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic2D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    material->properties(jmaterial);

    mpm::Material<Dim>::Matrix6x6 de = material->elastic_tensor();

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
    strain(2) = 0.0000000;
    strain(3) = 0.0000000;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Compute updated stress
    stress = material->compute_stress(stress, strain);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise strain
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0000000;
    strain(3) = 0.0000100;
    strain(4) = 0.0000000;
    strain(5) = 0.0000000;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress = material->compute_stress(stress, strain);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));

    // Add particle
    mpm::Index pid = 0;
    Eigen::Matrix<double, Dim, 1> coords;
    coords << 0.75, 0.75;
    auto particle = std::make_shared<mpm::Particle<Dim, 1>>(pid, coords);

    // Reset stress
    stress.setZero();
    // Compute updated stress
    stress = material->compute_stress(stress, strain, particle.get());

    // Check stressees
    REQUIRE(stress(0) == Approx(1.63461538461538e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.25000000000000e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(0.86538461538462e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.00000000000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.00000000000000e+00).epsilon(Tolerance));
  }
}


//! Check linearelastic class in 3D
TEST_CASE("LinearElastic is checked in 3D", "[material][linear_elastic][3D]") {
  // Tolerance
  const double Tolerance = 1.E-7;

  const unsigned Dim = 3;
  
  //! Check for id = 0
  SECTION("LinearElastic id is zero") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic3D", std::move(id));
    REQUIRE(material->id() == 0);
  }

  SECTION("LinearElastic id is positive") {
    //! Check for id is a positive value
    unsigned id = std::numeric_limits<unsigned>::max();
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic3D", std::move(id));
    REQUIRE(material->id() == std::numeric_limits<unsigned>::max());
  }

  //! Read material properties
  SECTION("LinearElastic check stiffness matrix") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

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
    const double a1 = 13461538.461566667;
    const double a2 = 5769230.769166667;

    mpm::Material<Dim>::Matrix6x6 de = material->elastic_tensor();
    REQUIRE(de(0, 0) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(0, 1) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(0, 2) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(0, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(0, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(0, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(1, 0) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(1, 1) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(1, 2) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(1, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(1, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(1, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(2, 0) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(2, 1) == Approx(a2).epsilon(Tolerance));
    REQUIRE(de(2, 2) == Approx(a1).epsilon(Tolerance));
    REQUIRE(de(2, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(2, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(2, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(3, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 3) == Approx(G).epsilon(Tolerance));
    REQUIRE(de(3, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(3, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(4, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(4, 4) == Approx(G).epsilon(Tolerance));
    REQUIRE(de(4, 5) == Approx(0.).epsilon(Tolerance));

    REQUIRE(de(5, 0) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 1) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 2) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 3) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 4) == Approx(0.).epsilon(Tolerance));
    REQUIRE(de(5, 5) == Approx(G).epsilon(Tolerance));
  }

  SECTION("LinearElastic check stresses") {
    unsigned id = 0;
    auto material = Factory<mpm::Material<Dim>, unsigned>::instance()->create(
        "LinearElastic3D", std::move(id));
    REQUIRE(material->id() == 0);

    // Initialise material
    Json jmaterial;
    jmaterial["density"] = 1000.;
    jmaterial["youngs_modulus"] = 1.0E+7;
    jmaterial["poisson_ratio"] = 0.3;

    material->properties(jmaterial);

    mpm::Material<Dim>::Matrix6x6 de = material->elastic_tensor();

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

    // Compute updated stress
    stress = material->compute_stress(stress, strain);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(0.000000e+00).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(0.000000e+00).epsilon(Tolerance));

    // Initialise strain
    strain(0) = 0.0010000;
    strain(1) = 0.0005000;
    strain(2) = 0.0005000;
    strain(3) = 0.0000100;
    strain(4) = 0.0000200;
    strain(5) = 0.0000300;

    // Reset stress
    stress.setZero();

    // Compute updated stress
    stress = material->compute_stress(stress, strain);

    // Check stressees
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(7.69230769230769e+01).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(1.15384615384615e+02).epsilon(Tolerance));

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
    REQUIRE(stress(0) == Approx(1.92307692307333e+04).epsilon(Tolerance));
    REQUIRE(stress(1) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(2) == Approx(1.53846153845333e+04).epsilon(Tolerance));
    REQUIRE(stress(3) == Approx(3.84615384615385e+01).epsilon(Tolerance));
    REQUIRE(stress(4) == Approx(7.69230769230769e+01).epsilon(Tolerance));
    REQUIRE(stress(5) == Approx(1.15384615384615e+02).epsilon(Tolerance));
  }
}
