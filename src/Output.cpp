#include <fstream>
#include <iostream>
#include <cstdio>
#include <ctime>

#include "Output.hpp"

// constructer
Output::Output() {}

// Output simulation parameters to file
void Output::outputParameters(const simulation_Parameters* param, const Gradient_constant* constant, const linklist_constant* linklist, int NumberOfParticles) {
    FILE* fp;
    char fileName[1024];
    sprintf(fileName, "%s/input_parameters.txt", parameters->output_dir);
    fp = fopen(fileName, "w");

    // Write simulation parameters
    fprintf(fp, "*******input parameter *********\n");
    fprintf(fp, "DIM : %d\n", parameters->DIM);
    fprintf(fp, "ParticleDistance : %lf\n", parameters->ParticleDistance);
    fprintf(fp, "Rho_fluid  : %lf\n", parameters->Rho[FLUID]);
    fprintf(fp, "Rho_solid  : %lf\n", parameters->Rho[SOLID]);
    fprintf(fp, "Rho_wall  : %lf\n", parameters->Rho[WALL]);
    fprintf(fp, "Re_inv fluid  : %lf\n", parameters->Re_inv[FLUID]);
    fprintf(fp, "Re_inv solid  : %lf\n", parameters->Re_inv[SOLID]);
    fprintf(fp, "Re_inv wall : %lf\n", parameters->Re_inv[WALL]);
    fprintf(fp, "strain_rate : %lf\n", parameters->strain_rate);
    fprintf(fp, "wall_speed : %lf\n", parameters->wall_speed);
    fprintf(fp, "Gravity : %lf\n", parameters->Gravity);
    fprintf(fp, "DT : %lf\n", parameters->DT);
    fprintf(fp, "FINISH_TIME : %lf\n", parameters->FINISH_TIME);
    fprintf(fp, "OUTPUT_INTERVAL : %d\n", parameters->OUTPUT_INTERVAL);
    fprintf(fp, "output_dir : %s\n", parameters->output_dir);
    fprintf(fp, "Init_place : %s\n", parameters->Init_place);
    fprintf(fp, "ks : %lf\n", parameters->ks);
    fprintf(fp, "kb : %lf\n", parameters->kb);
    fprintf(fp, "VTU_flag : %d\n", parameters->VTU_flag);
    fprintf(fp, "\n");

    // Write constants
    fprintf(fp, "*******constant parameter *********\n");
    fprintf(fp, "Re_forNumberDensity %lf\n", constant->Re_forNumberDensity);
    fprintf(fp, "Re_forGradient : %lf\n", constant->Re_forGradient);
    fprintf(fp, "Re_forLaplacian : %lf\n", constant->Re_forLaplacian);
    fprintf(fp, "N0_forNumberDensity : %lf\n", constant->N0_forNumberDensity);
    fprintf(fp, "N0_forGradient : %lf\n", constant->N0_forGradient);
    fprintf(fp, "N0_forLaplacian : %lf\n", constant->N0_forLaplacian);
    fprintf(fp, "Lambda : %lf\n", constant->Lambda);
    fprintf(fp, "\n");

    // Write linklist constants
    fprintf(fp, "******** linklist parameters ***********\n");
    fprintf(fp, "DB  : %lf\n", linklist->DB);
    fprintf(fp, "nBx : %d\n", linklist->nBx);
    fprintf(fp, "nBy : %d\n", linklist->nBy);
    fprintf(fp, "nBz : %d\n", linklist->nBz);
    fprintf(fp, "nBxy : %d\n", linklist->nBxy);
    fprintf(fp, "nBxyz : %d\n", linklist->nBxyz);
    fprintf(fp, "LX : %lf\n", linklist->Lx);
    fprintf(fp, "LY : %lf\n", linklist->Ly);
    fprintf(fp, "LZ : %lf\n", linklist->Lz);
    fprintf(fp, "MIN_X : %lf\n", linklist->MIN_X);
    fprintf(fp, "MIN_Y : %lf\n", linklist->MIN_Y);
    fprintf(fp, "MIN_Z : %lf\n", linklist->MIN_Z);
    fprintf(fp, "MAX_X : %lf\n", linklist->MAX_X);
    fprintf(fp, "MAX_Y : %lf\n", linklist->MAX_Y);
    fprintf(fp, "MAX_Z : %lf\n", linklist->MAX_Z);

    fclose(fp);
}

// Output analysis results to file
void Output::outputAnalysis(int n, double ParticleDistance, int iTimeStep, const simulation_Parameters* param, int n_fluid, int n_solid, double Time, int ny, clock_t start) {
    FILE* fp;
    char outputfile[1024];
    double width = (ny - 5) * ParticleDistance;
    clock_t end = clock();

    sprintf(outputfile, "%s/output_analysis.txt", parameters->output_dir);
    fp = fopen(outputfile, "w");

    fprintf(fp, "NumberOfParticle : %d\n", n);
    fprintf(fp, "IterationTime : %d\n", iTimeStep);
    fprintf(fp, "n_fluid : %d\n", n_fluid);
    fprintf(fp, "n_solid : %d\n", n_solid);
    fprintf(fp, "width_y : %lf\n", width);
    fprintf(fp, "Past_Time : %lf\n", Time);
    fprintf(fp, "Past_Time(real) : %lf\n", (double)(end - start) / CLOCKS_PER_SEC);

    fclose(fp);
}

// Output particle data to a file (e.g., CSV or other formats)
void Output::outputParticleData(const std::string& filename, const std::vector<Particle>& particles) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    // Output particle positions and velocities
    file << "id,position_x,position_y,position_z,velocity_x,velocity_y,velocity_z\n";
    for (size_t i = 0; i < particles.size(); ++i) {
        file << i << ","
             << particles[i].position[0] << "," << particles[i].position[1] << "," << particles[i].position[2] << ","
             << particles[i].velocity[0] << "," << particles[i].velocity[1] << "," << particles[i].velocity[2] << "\n";
    }

    file.close();
}
