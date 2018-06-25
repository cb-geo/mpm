#ifndef MPM_LOGGER_H_
#define MPM_LOGGER_H_

// Speed log
#include "spdlog/spdlog.h"

//! MPM namespace
namespace mpm {

// Create a logger for reading mesh
auto io_logger = spdlog::stdout_color_st("IO");

// Create a logger for reading mesh
auto read_mesh_logger = spdlog::stdout_color_st("ReadMesh");

// Create a logger for reading ascii mesh
auto read_mesh_ascii_logger = spdlog::stdout_color_st("ReadMeshAscii");

// Create a logger for MPM
auto mpm_logger = spdlog::stdout_color_st("MPM");

// Create a logger for MPM
auto mpm_explicit_logger = spdlog::stdout_color_st("MPMExplicit");

}  // namespace mpm

#endif  // MPM_LOGGER_H_
