#include <cmath>
#include <limits>
#include <memory>

#include "Eigen/Dense"
#include "catch.hpp"

#include "function_base.h"
#include "geometry.h"
#include "node.h"

// Check node class for 1D case
TEST_CASE("Implicit Linear Node is checked for 1D case",
          "[node][1D][Implicit]") {
  const unsigned Dim = 1;
  const unsigned Dof = 1;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;
  Eigen::Matrix<double, 1, 1> coords;
  coords.setZero();

  mpm::Index id = 0;
  const double Tolerance = 1.E-7;
  std::shared_ptr<mpm::NodeBase<Dim>> node =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

  double mass = 0.;

  SECTION("Check update of velocity and acceleration using Newmark scheme") {
    // Time step
    const double dt = 0.1;
    // Parameters of Newmark scheme
    const double newmark_beta = 0.25;
    const double newmark_gamma = 0.50;

    // Initialization
    node->initialise_implicit();

    // Nodal displacement
    // clang-format off
    Eigen::VectorXd nodal_displacement;
    nodal_displacement.resize(Dim);
    nodal_displacement << 0.1;
    // clang-format on
    // Assign nodal active_id
    node->assign_active_id(0);
    // Assign nodal displacement
    node->update_displacement_increment(nodal_displacement, Nphase, 1);
    // Check nodal displacement
    REQUIRE(node->displacement(Nphase)(0) ==
            Approx(nodal_displacement(0)).epsilon(Tolerance));

    // Nodal velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    for (unsigned i = 0; i < Dim; ++i)
      node->assign_velocity_constraint(i, velocity(i));
    node->apply_velocity_constraints();

    // Update nodal velocity and acceleration using Newmark scheme
    node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                               newmark_gamma, dt);

    // Values of nodal velocity
    Eigen::Matrix<double, 1, 1> nodal_velocity;
    // clang-format off
    nodal_velocity << 2.;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
      REQUIRE(node->velocity(Nphase)(i) ==
              Approx(nodal_velocity(0, i)).epsilon(Tolerance));

    // Values of nodal acceleration
    Eigen::Matrix<double, 1, 1> nodal_acceleration;
    // clang-format off
    nodal_acceleration << 40.;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) ==
              Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

    // Exception check when mass is zero
    mass = 0.;
    // Update mass to 0.
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
    REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
  }

  SECTION("Check momentum, inertia, velocity and acceleration") {
    // Initialization
    node->initialise_implicit();

    // Check momentum
    Eigen::Matrix<double, Dim, 1> momentum;
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

    // Check initial momentum
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update momentum to 20
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check inertia
    Eigen::Matrix<double, Dim, 1> inertia;
    for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

    // Check initial inertia
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update inertia to 20
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check zero velocity
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check zero acceleration
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check mass
    double mass = 0.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    // Compute and check velocity and acceleration this should throw zero mass
    node->compute_velocity_acceleration();

    mass = 100.;
    // Update mass to 100.
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

    // Compute and check velocity and acceleration
    node->compute_velocity_acceleration();
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
  }
}

// \brief Check node class for 2D case
TEST_CASE("Implicit Linear Node is checked for 2D case",
          "[node][2D][Implicit]") {
  const unsigned Dim = 2;
  const unsigned Dof = 2;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;
  Eigen::Vector2d coords;
  coords.setZero();

  mpm::Index id = 0;
  const double Tolerance = 1.E-7;
  std::shared_ptr<mpm::NodeBase<Dim>> node =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

  double mass = 0.;

  SECTION("Check update of velocity and acceleration using Newmark scheme") {
    // Time step
    const double dt = 0.1;
    // Parameters of Newmark scheme
    const double newmark_beta = 0.25;
    const double newmark_gamma = 0.50;

    // Initialization
    node->initialise_implicit();

    // Nodal displacement
    // clang-format off
    Eigen::VectorXd nodal_displacement;
    nodal_displacement.resize(Dim);
    nodal_displacement << 0.,
                          0.1;
    // clang-format on
    // Assign nodal active_id
    node->assign_active_id(0);
    // Assign nodal displacement
    node->update_displacement_increment(nodal_displacement, Nphase, 1);
    // Check nodal displacement
    for (unsigned i = 0; i < nodal_displacement.size(); ++i)
      REQUIRE(node->displacement(Nphase)(i) ==
              Approx(nodal_displacement(i)).epsilon(Tolerance));

    // Nodal velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    for (unsigned i = 0; i < Dim; ++i)
      node->assign_velocity_constraint(i, velocity(i));
    node->apply_velocity_constraints();

    // Update nodal velocity and acceleration using Newmark scheme
    node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                               newmark_gamma, dt);

    // Values of nodal velocity
    Eigen::Matrix<double, 2, 1> nodal_velocity;
    // clang-format off
    nodal_velocity << 0.,
                      1.;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
      REQUIRE(node->velocity(Nphase)(i) ==
              Approx(nodal_velocity(0, i)).epsilon(Tolerance));

    // Values of nodal acceleration
    Eigen::Matrix<double, 2, 1> nodal_acceleration;
    // clang-format off
    nodal_acceleration << 0.,
                          0.;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) ==
              Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

    // Exception check when mass is zero
    mass = 0.;
    // Update mass to 0.
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
    REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
  }

  SECTION("Check momentum, inertia, velocity and acceleration") {
    // Time step
    const double dt = 0.1;
    // Initialization
    node->initialise_implicit();

    // Check momentum
    Eigen::Matrix<double, Dim, 1> momentum;
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

    // Check initial momentum
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update momentum to 20
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check inertia
    Eigen::Matrix<double, Dim, 1> inertia;
    for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

    // Check initial inertia
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update inertia to 20
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check zero velocity
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check zero acceleration
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check mass
    double mass = 0.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    // Compute and check velocity this should throw zero mass
    node->compute_velocity_acceleration();

    mass = 100.;
    // Update mass to 100.
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

    // Compute and check velocity and acceleration
    node->compute_velocity_acceleration();
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
  }
}

// \brief Check node class for 3D case
TEST_CASE("Implicit Linear Node is checked for 3D case",
          "[node][3D][Implicit]") {
  const unsigned Dim = 3;
  const unsigned Dof = 3;
  const unsigned Nphases = 1;
  const unsigned Nphase = 0;

  Eigen::Vector3d coords;
  coords.setZero();

  mpm::Index id = 0;
  const double Tolerance = 1.E-7;
  std::shared_ptr<mpm::NodeBase<Dim>> node =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coords);

  double mass = 0.;

  SECTION("Check update of velocity and acceleration using Newmark scheme") {
    // Time step
    const double dt = 0.1;
    // Parameters of Newmark scheme
    const double newmark_beta = 0.25;
    const double newmark_gamma = 0.50;

    // Initialization
    node->initialise_implicit();

    // Nodal displacement
    // clang-format off
    Eigen::VectorXd nodal_displacement;
    nodal_displacement.resize(Dim);
    nodal_displacement << 0.,
                          0.1,
                          0.2;
    // clang-format on
    // Assign nodal active_id
    node->assign_active_id(0);
    // Assign nodal displacement
    node->update_displacement_increment(nodal_displacement, Nphase, 1);
    // Check nodal displacement
    for (unsigned i = 0; i < nodal_displacement.size(); ++i)
      REQUIRE(node->displacement(Nphase)(i) ==
              Approx(nodal_displacement(i)).epsilon(Tolerance));

    // Nodal velocity
    Eigen::VectorXd velocity;
    velocity.resize(Dim);
    for (unsigned i = 0; i < velocity.size(); ++i) velocity(i) = i;
    for (unsigned i = 0; i < Dim; ++i)
      node->assign_velocity_constraint(i, velocity(i));
    node->apply_velocity_constraints();

    // Update nodal velocity and acceleration using Newmark scheme
    node->update_velocity_acceleration_newmark(Nphase, newmark_beta,
                                               newmark_gamma, dt);

    // Values of nodal velocity
    Eigen::Matrix<double, 3, 1> nodal_velocity;
    // clang-format off
    nodal_velocity << 0.,
                      1.,
                      2.;
    // clang-format on
    // Check nodal velocity
    for (unsigned i = 0; i < nodal_velocity.cols(); ++i)
      REQUIRE(node->velocity(Nphase)(i) ==
              Approx(nodal_velocity(0, i)).epsilon(Tolerance));

    // Values of nodal acceleration
    Eigen::Matrix<double, 3, 1> nodal_acceleration;
    // clang-format off
    nodal_acceleration << 0.,
                          0.,
                          0.;
    // clang-format on
    // Check nodal acceleration
    for (unsigned i = 0; i < nodal_acceleration.cols(); ++i)
      REQUIRE(node->acceleration(Nphase)(i) ==
              Approx(nodal_acceleration(0, i)).epsilon(Tolerance));

    // Exception check when mass is zero
    mass = 0.;
    // Update mass to 0.
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(mass).epsilon(Tolerance));
    REQUIRE(node->compute_acceleration_velocity(Nphase, dt) == false);
  }

  SECTION("Check momentum, inertia, velocity and acceleration") {
    // Time step
    const double dt = 0.1;
    // Initialization
    node->initialise_implicit();

    // Check momentum
    Eigen::Matrix<double, Dim, 1> momentum;
    for (unsigned i = 0; i < momentum.size(); ++i) momentum(i) = 10.;

    // Check initial momentum
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update momentum to 20
    REQUIRE_NOTHROW(node->update_momentum(true, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign momentum to 10
    REQUIRE_NOTHROW(node->update_momentum(false, Nphase, momentum));
    for (unsigned i = 0; i < momentum.size(); ++i)
      REQUIRE(node->momentum(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check inertia
    Eigen::Matrix<double, Dim, 1> inertia;
    for (unsigned i = 0; i < inertia.size(); ++i) inertia(i) = 10.;

    // Check initial inertia
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Check update inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check update inertia to 20
    REQUIRE_NOTHROW(node->update_inertia(true, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(20.).epsilon(Tolerance));

    // Check assign inertia to 10
    REQUIRE_NOTHROW(node->update_inertia(false, Nphase, inertia));
    for (unsigned i = 0; i < inertia.size(); ++i)
      REQUIRE(node->inertia(Nphase)(i) == Approx(10.).epsilon(Tolerance));

    // Check mass
    double mass = 0.;
    REQUIRE_NOTHROW(node->update_mass(false, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(0.0).epsilon(Tolerance));
    // Compute and check velocity this should throw zero mass
    node->compute_velocity();

    mass = 100.;
    // Update mass to 100.5
    REQUIRE_NOTHROW(node->update_mass(true, Nphase, mass));
    REQUIRE(node->mass(Nphase) == Approx(100.).epsilon(Tolerance));

    // Check zero velocity
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.).epsilon(Tolerance));
    // Check zero acceleration
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.).epsilon(Tolerance));

    // Compute and check velocity and acceleration
    node->compute_velocity_acceleration();
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->velocity(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
    for (unsigned i = 0; i < Dim; ++i)
      REQUIRE(node->acceleration(Nphase)(i) == Approx(0.1).epsilon(Tolerance));
  }
}
