#include <array>
#include <iostream>
#include <memory>

#include "container.h"
#include "handler.h"
#include "node_base.h"

#include "Eigen/Dense"

int main(int argc, char** argv) {
  unsigned long long id = 0;
  const unsigned Dim = 3;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();

  auto nodebase = std::make_shared<mpm::NodeBase<Dim>>(id, coord);
  std::cout << "NodeBase id: " << nodebase->id() << '\n';

  auto nodebasehandler = std::make_shared<mpm::Handler<mpm::NodeBase<Dim>>>();
  nodebasehandler->insert(nodebase);

  for (auto itr = nodebasehandler->begin(); itr != nodebasehandler->end();
       ++itr)
    std::cout << ((*itr).second)->id() << '\n';

  auto nodebasecontainer =
      std::make_shared<mpm::Container<mpm::NodeBase<Dim>>>();
  nodebasecontainer->insert(nodebase);
  nodebasecontainer->insert(nodebase);
}
