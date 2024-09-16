#include <iostream>
#include "Parameters.hpp"

int main(int argc, char** argv) {
    // Check if an input file was provided as the first argument
    if (argc < 2) {
        std::cerr << "Error: No input file specified." << std::endl;
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    // Get the input file name from the first argument
    std::string input_filename = argv[1];

    // Create a Parameters object and load the parameters from the file
    Parameters params(input_filename);

    // Output some of the loaded parameters to verify
    std::cout << "Loaded parameters from: " << input_filename << std::endl;
    std::cout << "ParticleDistance: " << params.ParticleDistance << std::endl;
    std::cout << "DIM: " << params.DIM << std::endl;
    std::cout << "strain_rate: " << params.strain_rate << std::endl;
    std::cout << "wall_speed: " << params.wall_speed << std::endl;
    std::cout << "Gravity: " << params.Gravity << std::endl;
    std::cout << "VTU_flag: " << params.VTU_flag << std::endl;

    // Program execution can continue using the loaded parameters
    // ...

    return 0;
}
