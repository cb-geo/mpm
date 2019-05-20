#include "node.h"
#include "factory.h"
#include "node_base.h"

// Node2D (2 DoF, 1 Phase)
static Register<mpm::NodeBase<2>, mpm::Node<2, 2, 1>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    node2d("N2D");

// Node3D (3 DoF, 1 Phase)
static Register<mpm::NodeBase<3>, mpm::Node<3, 3, 1>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    node3d("N3D");

// Node2D (2 DoF, 2 Phase)
static Register<mpm::NodeBase<2>, mpm::Node<2, 2, 2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    node2d2p("N2D2P");

// Node3D (3 DoF, 2 Phase)
static Register<mpm::NodeBase<3>, mpm::Node<3, 3, 2>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    node3d2p("N3D2P");
