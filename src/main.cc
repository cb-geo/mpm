#include <array>
#include <iostream>
#include <memory>

#include "node.h"
#include "node_handler.h"

#include "Eigen/Dense"

int main(int argc, char** argv) {
  unsigned long long id = 0;
  const unsigned Dim = 3;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();

  auto node = std::make_shared<mpm::Node<Dim>>(id, coord);
  std::cout << "Node id: " << node->id() << '\n';

  auto nodehandler = std::make_shared<mpm::NodeHandler<Dim>>();
  nodehandler->insert_node(node);

  for (auto itr = nodehandler->nodes_begin(); itr != nodehandler->nodes_end(); ++itr)
    std::cout << ((*itr).second)->id() << '\n';
}
