#ifdef USE_PARTIO
#include "partio_writer.h"

// Write particles
bool mpm::partio::write_particles(
    const std::string& filename,
    const std::vector<mpm::PODParticle>& particles) {
  bool status = false;

  if (!particles.empty()) {
    Partio::ParticlesDataMutable* parts = Partio::create();
    Partio::ParticleAttribute postition, velocity, mass, index;

    index = parts->addAttribute("index", Partio::INT, 1);
    mass = parts->addAttribute("m", Partio::VECTOR, 1);
    postition = parts->addAttribute("position", Partio::VECTOR, 3);
    velocity = parts->addAttribute("v", Partio::VECTOR, 3);

    for (const auto& particle : particles) {
      // Add particle and get index
      int idx = parts->addParticle();
      // Index
      int* index_p = parts->dataWrite<int>(index, idx);
      index_p[0] = particle.id;

      // Write mass
      float* mass_p = parts->dataWrite<float>(mass, idx);
      mass_p[0] = particle.mass;

      // Write velocity and position
      float* velocity_p = parts->dataWrite<float>(velocity, idx);
      velocity_p[0] = particle.velocity_x;
      velocity_p[1] = particle.velocity_y;
      velocity_p[2] = particle.velocity_z;

      float* position_p = parts->dataWrite<float>(postition, idx);
      position_p[0] = particle.coord_x;
      position_p[1] = particle.coord_y;
      position_p[2] = particle.coord_z;
    }
    Partio::write(filename.c_str(), *parts);
    parts->release();

    status = true;
  }
  return status;
}
#endif
