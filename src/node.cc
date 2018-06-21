#include "node.h"
#include "factory.h"
#include "node_base.h"

static Register<mpm::NodeBase<2>, mpm::Node<2, 2, 1>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    node2d("Node2D");
