#include <array>
#include <iostream>
#include <memory>

#include "container.h"
#include "handler.h"
#include "node.h"

#include "Eigen/Dense"

int main(int argc, char** argv) {
  unsigned long long id = 0;
  const unsigned Dim = 3;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();

  auto node = std::make_shared<mpm::Node<Dim>>(id, coord);
  std::cout << "Node id: " << node->id() << '\n';

  auto nodehandler = std::make_shared<mpm::Handler<mpm::Node<Dim>>>();
  nodehandler->insert(node);

  for (auto itr = nodehandler->begin(); itr != nodehandler->end(); ++itr)
    std::cout << ((*itr).second)->id() << '\n';

  auto nodecontainer = std::make_shared<mpm::Container<mpm::Node<Dim>>>();
  nodecontainer->insert(node);
  nodecontainer->insert(node);

}
