#include "affine_transform.h"

// 1D
// M=[0 1; 1 1];
// K1 = transpose(M) * inverse (M*transpose(M));
template <>
const Eigen::Matrix<double, 2, 1> mpm::TransformR2UAffine<1, 2>::KA =
    (Eigen::Matrix<double, 2, 1>() << -1.000000f, 1.000000f).finished();

template <>
const Eigen::Matrix<double, 2, 1> mpm::TransformR2UAffine<1, 2>::Kb =
    (Eigen::Matrix<double, 2, 1>() << 1.000000f, 0.000000f).finished();

// 2D
// M=[0 1 0 1;0 0 1 1;1 1 1 1];
// K2 = transpose(M) * inverse (M*transpose(M));
// clang-format off
template <>
const Eigen::Matrix<double, 4, 2> mpm::TransformR2UAffine<2, 4>::KA =
    (Eigen::Matrix<double, 4, 2>() <<
     -0.500000f, -0.500000f,
      0.500000f, -0.500000f,
     -0.500000f,  0.500000f,
      0.500000f,  0.500000f).finished();
// clang-format on

// M=[0 1 0 1 0 1 0 1;0 0 1 1 0 0 1 1; 0 0 0 0 1 1 1 1; 1 1 1 1 1 1 1 1];
// K3 = transpose(M) * inverse(M*transpose(M))
template <>
const Eigen::Matrix<double, 4, 1> mpm::TransformR2UAffine<2, 4>::Kb =
    (Eigen::Matrix<double, 4, 1>() << 0.750000f, 0.250000f, 0.250000f,
     -0.250000f)
        .finished();

// 3D
// clang-format off
template <>
const Eigen::Matrix<double, 8, 3> mpm::TransformR2UAffine<3, 8>::KA =
    (Eigen::Matrix<double, 8, 3>() <<
     -0.250000f, -0.250000f, -0.250000f,
      0.250000f, -0.250000f, -0.250000f,
     -0.250000f,  0.250000f, -0.250000f,
      0.250000f,  0.250000f, -0.250000f,
     -0.250000f, -0.250000f,  0.250000f,
      0.250000f, -0.250000f,  0.250000f,
     -0.250000f,  0.250000f,  0.250000f,
      0.250000f,  0.250000f,  0.250000f).finished();
// clang-format on

template <>
const Eigen::Matrix<double, 8, 1> mpm::TransformR2UAffine<3, 8>::Kb =
    (Eigen::Matrix<double, 8, 1>() << 0.500000f, 0.250000f, 0.250000f,
     0.000000f, 0.250000f, 0.000000f, 0.000000f, -0.250000f)
        .finished();
