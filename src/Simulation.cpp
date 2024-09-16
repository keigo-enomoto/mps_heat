#include <cmath>
#include <iostream>

#include "Simulation.hpp"
#include "Parameters.hpp"

// Function to apply periodic boundary conditions to displacement
void Simulation::applyPBC(std::vector<double>& displacement, const CellList& celllist) {
    if (displacement[0] < -celllist.Lx * 0.5) displacement[0] += celllist.Lx;
    else if (displacement[0] > celllist.Lx * 0.5) displacement[0] -= celllist.Lx;

    if (displacement[1] < -celllist.Ly * 0.5) displacement[1] += celllist.Ly;
    else if (displacement[1] > celllist.Ly * 0.5) displacement[1] -= celllist.Ly;

    if (displacement[2] < -celllist.Lz * 0.5) displacement[2] += celllist.Lz;
    else if (displacement[2] > celllist.Lz * 0.5) displacement[2] -= celllist.Lz;
}

// Function to calculate viscosity
void Simulation::calViscosity(Parameters* param, Gradient_constant* constant, CellList* celllist) {
    double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z;
    double distance, w, a;
    a = (2.0 * parameters->DIM) / (constant->N0_forLaplacian * constant->Lambda);

    #pragma omp parallel for collapse(3)
    for (int b_ix = 0; b_ix < celllist->grid_size_x; ++b_ix) {
        for (int b_iy = 0; b_iy < celllist->grid_size_y; ++b_iy) {
            for (int b_iz = 0; b_iz < celllist->grid_size_z; ++b_iz) {
                for (int i = celllist->BucketFirst[b_ix][b_iy][b_iz]; i != -1; i = celllist->Nextof[i]) {
                    Particle& particle_i = particles[i];
                    viscosityTerm_x = viscosityTerm_y = viscosityTerm_z = 0.0;

                    for (int j = 0; j < particles.size(); ++j) {
                        if (j == i) continue;
                        Particle& particle_j = particles[j];
                        distance = particle_i.calculateDistancePBC(particle_j, *celllist);
                        w = Particle::calculateWeight(distance, constant->Re_forLaplacian);
                        viscosityTerm_x += (particle_j.vel[0] - particle_i.vel[0]) * w;
                        viscosityTerm_y += (particle_j.vel[1] - particle_i.vel[1]) * w;
                        viscosityTerm_z += (particle_j.vel[2] - particle_i.vel[2]) * w;
                    }

                    particle_i.acl[0] += a * viscosityTerm_x;
                    particle_i.acl[1] += a * viscosityTerm_y;
                    particle_i.acl[2] += a * viscosityTerm_z;
                }
            }
        }
    }
}

// Function to handle collisions between particles
void Simulation::collision(Parameters* param, CellList* celllist) {
    const double e = 0.2;  // Coefficient of restitution
    const double dt = parameters->DT;
    const double collisionDistance2 = 0.25 * (parameters->ParticleDistance * parameters->ParticleDistance);

    #pragma omp parallel for
    for (int i = 0; i < particles.size(); ++i) {
        Particle& particle_i = particles[i];
        for (int j = 0; j < particles.size(); ++j) {
            if (i == j) continue;
            Particle& particle_j = particles[j];

            std::vector<double> displacement(3);
            for (int d = 0; d < 3; ++d) {
                displacement[d] = particle_j.pos[d] - particle_i.pos[d];
            }
            applyPBC(displacement, *celllist);

            double distance2 = displacement[0] * displacement[0] +
                               displacement[1] * displacement[1] +
                               displacement[2] * displacement[2];
            if (distance2 < collisionDistance2) {
                double mi = particle_i.rho;
                double mj = particle_j.rho;
                double forceDT = (1.0 + e) * ((displacement[0] * (particle_i.vel[0] - particle_j.vel[0])) +
                                              (displacement[1] * (particle_i.vel[1] - particle_j.vel[1])) +
                                              (displacement[2] * (particle_i.vel[2] - particle_j.vel[2]))) /
                                 (dt * (mi + mj));

                // Update velocities
                for (int d = 0; d < 3; ++d) {
                    particle_i.vel[d] -= forceDT * (mj / (mi + mj)) * displacement[d];
                    particle_j.vel[d] += forceDT * (mi / (mi + mj)) * displacement[d];
                }
            }
        }
    }
}
