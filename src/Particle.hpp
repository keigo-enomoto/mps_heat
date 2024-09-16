#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include "CellList.hpp"

enum{

 GHOST = -1,
 FLUID =  0,
 SOLID = 1,
 WALL  =  2,
 DUMMY_WALL = 3

};


/**************************************************************
  *                  Particles                                *
***************************************************************/

class Particle {
public:
    std::vector<double> pos;   // 3D position vector     [m]
    std::vector<double> vel;   // 3D velocity vector     [m/s]
    std::vector<double> acl;   // 3D acceleration vector [m/s^2]
    double rho;                // mass density           [kg/m^3]
    double press;              // pressure               [Pa]
    double n;                  // number density



    // Constructor
    Particle(const std::vector<double>& pos, const std::vector<double>& vel, double mass);

    // Methods
    void update(double dt);
    void applyForce(const std::vector<double>& force);
    void resetAcceleration();

    // Calculate distance between this particle and another, applying PBC
    double calculateDistancePBC(const Particle& other, const CellList& celllist) const;

    // Calculate the weight function based on distance
    static double calculateWeight(double distance, double re);
};

#endif // PARTICLE_H
