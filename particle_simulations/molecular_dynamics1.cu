#include <iostream>
#include <map>
#include <vector>
#include "configparser.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configfile>" << std::endl;
        return 1;
    }

    std::map<std::string, std::string> config;
    std::vector<std::vector<float>> particlePositions;
    std::vector<std::vector<float>> particleVelocities;

    try {
        parseConfigFile(argv[1], config, particlePositions, particleVelocities);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    // Process the parsed configuration
    int numTimeSteps = std::stoi(config["num_time_step"]);
    int timeStepLength = std::stoi(config["time_step_lenth"]);
    float sigma = std::stof(config["sigma"]);
    float epsilon = std::stof(config["epsilon"]);
    int particleNum = std::stoi(config["particle_num"]);
    float particleMass = std::stof(config["particle_mass"]);

    // Output the parsed data for verification
    std::cout << "Simulation Parameters:" << std::endl;
    std::cout << "Number of Time Steps: " << numTimeSteps << std::endl;
    std::cout << "Time Step Length: " << timeStepLength << std::endl;
    std::cout << "Sigma: " << sigma << std::endl;
    std::cout << "Epsilon: " << epsilon << std::endl;
    std::cout << "Number of Particles: " << particleNum << std::endl;
    std::cout << "Particle Mass: " << particleMass << std::endl;

    std::cout << "\nParticle Positions:" << std::endl;
    for (const auto& pos : particlePositions) {
        std::cout << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")" << std::endl;
    }

    std::cout << "\nParticle Velocities:" << std::endl;
    for (const auto& vel : particleVelocities) {
        std::cout << "(" << vel[0] << ", " << vel[1] << ", " << vel[2] << ")" << std::endl;
    }

    return 0;
}