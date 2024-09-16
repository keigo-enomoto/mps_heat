#include "PressureSolver.hpp"
#include <Eigen/IterativeLinearSolvers>
#include <iostream>

PressureSolver::PressureSolver(int size)
    : A(size, size), b(size), pressure(size) {}

void PressureSolver::setupSystem(const std::vector<Particle>& particles) {
    // Setup sparse matrix A and vector b using particle positions, etc.
}

void PressureSolver::solve() {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    pressure = solver.compute(A).solve(b);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Pressure solver failed!" << std::endl;
    }
}
