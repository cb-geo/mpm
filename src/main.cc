#include <array>
#include <iostream>
#include <memory>

#include "node.h"
#include "node_base.h"
#include "explicitonephase_mpm.h"


#include "Eigen/Dense"

int main(int argc, char** argv) {
  unsigned long long id = 0;
  const unsigned Dim = 3;
  const unsigned Dof = 6;
  const unsigned Nphases = 1;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();

  std::shared_ptr<mpm::NodeBase<Dim>> node =
      std::make_shared<mpm::Node<Dim, Dof, Nphases>>(id, coord);
  std::cout << "Node id: " << node->id() << '\n';

  std::shared_ptr<mpm::MPM<Dim>> explicit_solver =
    std::make_shared<mpm::MPMExplicitSinglePhase<Dim>>();
}
