#ifndef PRESSURE_SOLVER_H
#define PRESSURE_SOLVER_H

#include <vector>
#include <Eigen/Sparse>

#include "Particle.hpp"

class PressureSolver {
public:
    Eigen::SparseMatrix<double> A;  // Sparse matrix for system
    Eigen::VectorXd b;              // RHS vector
    Eigen::VectorXd pressure;       // Solution vector

    PressureSolver(int size);

    // Methods
    void setupSystem(const std::vector<Particle>& particles);
    void solve();
};

#endif // PRESSURE_SOLVER_H
