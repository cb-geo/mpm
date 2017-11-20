#include "node_base.h"


#include <array>
#include <iostream>
#include <memory>

int main(int argc, char** argv) {
  long long id = 0;
  const unsigned Dim = 3;
  std::array<double, Dim> coord = {{0.}};

  auto node = std::make_shared<mpm::NodeBase<Dim>>(id, coord);
  node->info();
}
