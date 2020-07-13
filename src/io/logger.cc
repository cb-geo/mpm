#include "logger.h"

// Create a logger for IO
const std::shared_ptr<spdlog::logger> mpm::Logger::io_logger =
    spdlog::stdout_color_st("IO");

// Create a logger for reading mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::io_mesh_logger =
    spdlog::stdout_color_st("IOMesh");

// Create a logger for reading ascii mesh
const std::shared_ptr<spdlog::logger> mpm::Logger::io_mesh_ascii_logger =
    spdlog::stdout_color_st("IOMeshAscii");

// Create a logger for point generator
const std::shared_ptr<spdlog::logger> mpm::Logger::point_generator_logger =
    spdlog::stdout_color_st("PointGenerator");

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

// Create a logger for MPM Explicit TwoPhase
const std::shared_ptr<spdlog::logger>
    mpm::Logger::mpm_explicit_twophase_logger =
        spdlog::stdout_color_st("MPMExplicitTwoPhase");