#include "logger.h"

// Create a logger for IO
const std::shared_ptr<spdlog::logger> mpm::Logger::io_logger =
    spdlog::stdout_color_st("IO");

// Create a logger for reading mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::read_mesh =
    spdlog::stdout_color_st("ReadMesh");

// Create a logger for reading ascii mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::read_mesh_ascii =
    spdlog::stdout_color_st("ReadMeshAscii");

// Create a logger for MPM
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_logger =
    spdlog::stdout_color_st("MPM");

// Create a logger for MPM Base
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_base_logger =
    spdlog::stdout_color_st("MPMBase");

// Create a logger for MPM Explicit
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_logger =
    spdlog::stdout_color_st("MPMExplicit");

// Create a logger for MPM Explicit USF
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_usf_logger =
    spdlog::stdout_color_st("MPMExplicitUSF");

// Create a logger for MPM Explicit USL
const std::shared_ptr<spdlog::logger> mpm::Logger::mpm_explicit_usl_logger =
    spdlog::stdout_color_st("MPMExplicitUSL");
