#include <Eigen/Dense>
#include <memory>

#include "catch.hpp"

#include "nodal_properties.h"

//! \brief Check NodalProperties struct and its functions
TEST_CASE("NodalProperties is checked", "[nodal_properties]") {
  // Number of nodes
  const unsigned nnodes = 3;
  // Number of materials
  const unsigned nmaterials = 2;
  // Tolerance
  const double tolerance = 1.E-7;

  // Check scalar property creation
  SECTION("Create scalar properties") {
    // Declare nodal properties
    mpm::NodalProperties nodal_properties;

    // Define properties to be created
    std::string property1 = "masses";
    std::string property2 = "areas";

    // Define dimension (1D for scalar)
    const unsigned dim = 1;

    // Check property creation
    REQUIRE(
        nodal_properties.create_property(property1, nnodes * dim, nmaterials));
    REQUIRE(
        nodal_properties.create_property(property2, nnodes * dim, nmaterials));

    // Check size of matrix of property data
    REQUIRE(nodal_properties.properties_.at(property1).rows() == 3);
    REQUIRE(nodal_properties.properties_.at(property1).cols() == 2);
    REQUIRE(nodal_properties.properties_.at(property1).size() == 6);
    REQUIRE(nodal_properties.properties_.at(property2).rows() == 3);
    REQUIRE(nodal_properties.properties_.at(property2).cols() == 2);
    REQUIRE(nodal_properties.properties_.at(property2).size() == 6);

    // Check if values in nodal properties are equal to 0.0
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
        }
      }
    }

    // clang-format off
    Eigen::Matrix<double,3,2> data1;
    data1 << 0.0, 3.3,
             0.5, 2.2,
             1.0, 1.1;
    Eigen::Matrix<double,3,2> data2;
    data2 << 2.3, 0.0,
             1.8, 0.2,
             0.1, 4.1;
    // clang-format on

    // Assign values in the property data matrix and check assignment
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        REQUIRE_NOTHROW(nodal_properties.assign_property(
            property1, i, j, data1.block<dim, 1>(i * dim, j), dim));
        REQUIRE_NOTHROW(nodal_properties.assign_property(
            property2, i, j, data2.block<dim, 1>(i * dim, j), dim));
      }
    }

    // Check values in matrix of property data
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(data1(i * dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(data2(i * dim + k, j)).epsilon(tolerance));
        }
      }
    }

    // clang-format off
    data1 << 4.2, 3.1,
             0.5, 1.7,
             6.6, 8.2;
    data2 << 4.1, 3.2,
             6.8, 5.1,
             9.1, 2.0;
    // clang-format on

    // Update values in the property data matrix and check update
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        REQUIRE_NOTHROW(nodal_properties.update_property(
            property1, i, j, data1.block<dim, 1>(i * dim, j), dim));
        REQUIRE_NOTHROW(nodal_properties.update_property(
            property2, i, j, data2.block<dim, 1>(i * dim, j), dim));
      }
    }

    // clang-format off
    data1 << 4.2, 6.4,
             1.0, 3.9,
             7.6, 9.3;
    data2 << 6.4, 3.2,
             8.6, 5.3,
             9.2, 6.1;
    // clang-format on

    // Check values in matrix of property data
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(data1(i * dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(data2(i * dim + k, j)).epsilon(tolerance));
        }
      }
    }

    // initialise all nodal properties (set all values to zero)
    nodal_properties.initialise_properties();

    // Check if all values for all properties are zero
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
        }
      }
    }
  }

  // Check vector property creation
  SECTION("Create vector properties") {
    // Declare nodal properties
    mpm::NodalProperties nodal_properties;

    // Define properties to be created
    std::string property1 = "velocities";
    std::string property2 = "momenta";

    // Define dimension (2D)
    const unsigned dim = 2;

    // Check property creation
    REQUIRE(
        nodal_properties.create_property(property1, nnodes * dim, nmaterials));
    REQUIRE(
        nodal_properties.create_property(property2, nnodes * dim, nmaterials));

    // Check size of matrix of property data
    REQUIRE(nodal_properties.properties_.at(property1).rows() == 6);
    REQUIRE(nodal_properties.properties_.at(property1).cols() == 2);
    REQUIRE(nodal_properties.properties_.at(property1).size() == 12);
    REQUIRE(nodal_properties.properties_.at(property2).rows() == 6);
    REQUIRE(nodal_properties.properties_.at(property2).cols() == 2);
    REQUIRE(nodal_properties.properties_.at(property2).size() == 12);

    // Check if values in nodal properties are equal to 0.0
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
        }
      }
    }

    // clang-format off
    Eigen::Matrix<double,6,2> data1;
    data1 <<  3.5, 1.1,
              4.2, 3.4, // end of first node
              0.0, 0.5,
             -2.1, 6.1, // end of second node
             -1.2, 1.1,
              2.3, 0.7; // end of third node
    Eigen::Matrix<double,6,2> data2;
    data2 <<  0.1, 3.2,
             -1.9, 1.3, // end first node
              0.0, 2.3,
              0.3, 2.2, // end of second node
              0.0, 1.3,
              2.1, 5.5; // end of third node
    // clang-format on

    // Assign values in the property data matrix and check assignment
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        REQUIRE_NOTHROW(nodal_properties.assign_property(
            property1, i, j, data1.block<dim, 1>(i * dim, j), dim));
        REQUIRE_NOTHROW(nodal_properties.assign_property(
            property2, i, j, data2.block<dim, 1>(i * dim, j), dim));
      }
    }

    // Check values in matrix of property data
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(data1(i * dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(data2(i * dim + k, j)).epsilon(tolerance));
        }
      }
    }

    // clang-format off
    data1 << 4.2, 3.1,
             1.0, 1.7, // end of first node
             2.2, 8.2,
             3.5, 2.5, // end of second node
             9.7, 1.1,
             6.4, 0.7; // end of third node
    data2 << 4.1, 3.2,
             7.7, 5.1, // end first node
             2.4, 2.0,
             2.5, 4.3, // end of second node
             9.8, 8.3,
             6.1, 0.1; // end of third node
    // clang-format on

    // Update values in the property data matrix and check update
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        REQUIRE_NOTHROW(nodal_properties.update_property(
            property1, i, j, data1.block<dim, 1>(i * dim, j), dim));
        REQUIRE_NOTHROW(nodal_properties.update_property(
            property2, i, j, data2.block<dim, 1>(i * dim, j), dim));
      }
    }

    // clang-format off
    data1 << 7.7, 4.2,
             5.2, 5.1, // end of first node
             2.2, 8.7,
             1.4, 8.6, // end of second node
             8.5, 2.2,
             8.7, 1.4; // end of third node
    data2 << 4.2, 6.4,
             5.8, 6.4, // end first node
             2.4, 4.3,
             2.8, 6.5, // end of second node
             9.8, 9.6,
             8.2, 5.6; // end of third node
    // clang-format on

    // Check values in matrix of property data
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(data1(i * dim + k, j)).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(data2(i * dim + k, j)).epsilon(tolerance));
        }
      }
    }

    // initialise all nodal properties (set all values to zero)
    nodal_properties.initialise_properties();

    // Check if all values for all properties are zero
    for (int i = 0; i < nnodes; ++i) {
      for (int j = 0; j < nmaterials; ++j) {
        for (int k = 0; k < dim; ++k) {
          REQUIRE(nodal_properties.property(property1, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
          REQUIRE(nodal_properties.property(property2, i, j, dim)(k, 0) ==
                  Approx(0.0).epsilon(tolerance));
        }
      }
    }
  }
}
