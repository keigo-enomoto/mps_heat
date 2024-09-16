#ifndef OUTPUT_H
#define OUTPUT_H

#include <vector>
#include <string>
#include "Particle.hpp"
#include "Parameters.hpp"

class Output {
public:
    Output();

    // Methods to output data
    // void outputParameters(const simulation_Parameters* param, const Gradient_constant* constant, const linklist_constant* linklist, int NumberOfParticles);
    // void outputAnalysis(int n, double ParticleDistance, int iTimeStep, const simulation_Parameters* param, int n_fluid, int n_solid, double Time, int ny, clock_t start);
    // void outputParticleData(const std::string& filename, const std::vector<Particle>& particles);
};

#endif // OUTPUT_H
