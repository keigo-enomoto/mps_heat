#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "Particle.hpp"
#include "Parameters.hpp"
#include "CellList.hpp"

/**************************************************************
  *         constats for Gradient and Laplacian model         *
***************************************************************/

// Parameters でもいいし、Kernelクラスを別で定義してもいいかも
// typedef struct Gradient_constant
// {
//     double Re_forNumberDensity ;
//     double Re_forGradient ;
//     double Re_forLaplacian ;
//     double N0_forNumberDensity ;
//     double N0_forGradient ;
//     double N0_forLaplacian ;
//     double Lambda ;

// } Gradient_constant;

/**************************************************************
  *                      Simulation Class                     *
***************************************************************/

class Simulation {
public:
    std::vector<Particle> particles;
    double dt;
    double gravity;

    // Constructor
    Simulation(double timeStep, double gravityAccel);

    // Methods
    void initializeParticles();
    void applyForces();
    void updateParticles();
    void run(double endTime);

    // // Apply periodic boundary conditions to displacement
    // void applyPBC(std::vector<double>& displacement, const CellList& celllist);

    // // Methods for viscosity calculation and collision handling
    // void calViscosity(Parameters* param, Gradient_constant* constant, CellList* celllist);
    // void collision(Parameters* param, CellList* celllist);
};

#endif // SIMULATION_H
