#pragma once
#include <string>
#include "InputParams.hpp"  // Include the InputParams class

class Parameters {
public:
    double ParticleDistance = 1.0;
    int DIM = 3;
    double strain_rate = 0.0;
    double wall_speed = 0.0;
    double Gravity = 9.81;
    double DT = 0.01;
    double FINISH_TIME = 10.0;
    int OUTPUT_INTERVAL = 100;
    std::string Init_place = ".";
    std::string output_dir = ".";
    int VTU_flag = 0;        // if 1, output vtu file, if 2, output only fiber and wall
    int stress_flag = 0;     // if 1, output stress
    int flow_flag = 0;       // if 1, output flow profile

    // Constructor that loads parameters using InputParams
    Parameters(const std::string& filename) {
        param::InputParams input_param(filename);

        // Use the get method from InputParams to initialize parameters
        ParticleDistance = input_param.get<double>("ParticleDistance", ParticleDistance);
        DIM = input_param.get<int>("DIM", DIM);
        strain_rate = input_param.get<double>("strain_rate", strain_rate);
        wall_speed = input_param.get<double>("wall_speed", wall_speed);
        Gravity = input_param.get<double>("Gravity", Gravity);
        DT = input_param.get<double>("DT", DT);
        FINISH_TIME = input_param.get<double>("FINISH_TIME", FINISH_TIME);
        OUTPUT_INTERVAL = input_param.get<int>("OUTPUT_INTERVAL", OUTPUT_INTERVAL);
        Init_place = input_param.get<std::string>("Init_place", Init_place);
        output_dir = input_param.get<std::string>("output_dir", output_dir);
        VTU_flag = input_param.get<int>("VTU_flag", VTU_flag);
        stress_flag = input_param.get<int>("stress_flag", stress_flag);
        flow_flag = input_param.get<int>("flow_flag", flow_flag);
    }
};
