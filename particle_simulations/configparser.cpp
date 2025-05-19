#include "configparser.h"
#include <fstream>
#include <sstream>
#include <regex>
#include <iostream>

Config ConfigParser::parse(const std::string& filename) {
    Config config;
    std::ifstream infile(filename);
    std::string line;

    std::regex number_regex(R"([+-]?\d*\.?\d+)");
    std::regex vector_regex(R"(\(\s*([+-]?\d*\.?\d+),\s*([+-]?\d*\.?\d+),\s*([+-]?\d*\.?\d+)\))");

    std::vector<Vector3> positions, velocities;
    std::vector<float> masses;

    while (std::getline(infile, line)) {
        if (line.find("num_time_step") != std::string::npos) {
            config.num_time_step = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("time_step_lenth") != std::string::npos) {
            config.time_step_length = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("sigma") != std::string::npos) {
            config.sigma = std::stof(line.substr(line.find("=") + 1));
        } else if (line.find("epsilon") != std::string::npos) {
            config.epsilon = std::stof(line.substr(line.find("=") + 1));
        } else if (line.find("particle_num") != std::string::npos) {
            config.particle_num = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("particle_positions") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, vector_regex)) {
                    positions.push_back({
                        std::stof(match[1]),
                        std::stof(match[2]),
                        std::stof(match[3])
                    });
                }
            }
        } else if (line.find("particle_velocities") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, vector_regex)) {
                    velocities.push_back({
                        std::stof(match[1]),
                        std::stof(match[2]),
                        std::stof(match[3])
                    });
                }
            }
        } else if (line.find("particle_masses") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, number_regex)) {
                    masses.push_back(std::stof(match[0]));
                }
            }
        }
    }

    // Zusammensetzen der Partikel
    size_t n = std::min({positions.size(), velocities.size(), masses.size()});
    for (size_t i = 0; i < n; ++i) {
        Particle p;
        p.position = positions[i];
        p.velocity = velocities[i];
        p.mass = masses[i];
        p.force = {0, 0, 0};         
        p.acceleration = {0, 0, 0};  
        config.particles.push_back(p);
    }

    return config;
}

// Ausgabeoperatoren

std::ostream& operator<<(std::ostream& os, const Vector3& vec) {
    os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
    return os;
}

std::ostream& operator<<(std::ostream& os, const Particle& particle) {
    os << "  Position: " << particle.position << "\n"
       << "  Velocity: " << particle.velocity << "\n"
       << "  Force: " << particle.force << "\n"
       << "  Acceleration: " << particle.acceleration << "\n"
       << "  Mass: " << particle.mass;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Config& config) {
    os << "Time Steps: " << config.num_time_step << "\n";
    os << "Time Step Length: " << config.time_step_length << "\n";
    os << "Sigma: " << config.sigma << "\n";
    os << "Epsilon: " << config.epsilon << "\n";
    os << "Particle Count: " << config.particle_num << "\n";

    for (size_t i = 0; i < config.particles.size(); ++i) {
        os << "Particle " << i + 1 << ":\n" << config.particles[i] << "\n";
    }

    return os;
}
