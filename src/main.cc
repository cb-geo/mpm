#include <array>
#include <iostream>
#include <memory>

#include "node.h"

#include "Eigen/Dense"

int main(int argc, char** argv) {
  long long id = 0;
  const unsigned Dim = 3;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();

  auto node = std::make_shared<mpm::Node<Dim>>(id, coord);
  node->info();
}
