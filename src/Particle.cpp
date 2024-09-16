#include <cmath>  // For sqrt

#include "Particle.hpp"
#include "CellList.hpp"

Particle::Particle(const std::vector<double>& pos, const std::vector<double>& vel, double m)
    : pos(pos), vel(vel), rho(m), acl(3, 0.0) {}

void Particle::update(double dt) {
    for (int i = 0; i < 3; ++i) {
        vel[i] += acl[i] * dt;
        pos[i] += vel[i] * dt;
    }
}

void Particle::applyForce(const std::vector<double>& force) {
    for (int i = 0; i < 3; ++i) {
        acl[i] += force[i] / rho;
    }
}

void Particle::resetAcceleration() {
    acl.assign(3, 0.0);
}

// Calculate the distance between this particle and another, applying PBC
double Particle::calculateDistancePBC(const Particle& other, const CellList& celllist) const {
    std::vector<double> displacement(3);
    for (int i = 0; i < 3; ++i) {
        displacement[i] = other.pos[i] - pos[i];
    }

    // Apply periodic boundary conditions
    if (displacement[0] < -celllist.Lx * 0.5) displacement[0] += celllist.Lx;
    else if (displacement[0] > celllist.Lx * 0.5) displacement[0] -= celllist.Lx;

    if (displacement[1] < -celllist.Ly * 0.5) displacement[1] += celllist.Ly;
    else if (displacement[1] > celllist.Ly * 0.5) displacement[1] -= celllist.Ly;

    if (displacement[2] < -celllist.Lz * 0.5) displacement[2] += celllist.Lz;
    else if (displacement[2] > celllist.Lz * 0.5) displacement[2] -= celllist.Lz;

    return sqrt(displacement[0] * displacement[0] +
                displacement[1] * displacement[1] +
                displacement[2] * displacement[2]);
}

// Calculate the weight function based on distance
double Particle::calculateWeight(double distance, double re) {
    if (distance >= re) return 0.0;
    return re / distance - 1.0;
}
