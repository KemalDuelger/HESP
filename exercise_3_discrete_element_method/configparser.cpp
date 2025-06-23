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
    std::vector<double> masses, radiuses;

    while (std::getline(infile, line)) {
        if (line.find("num_time_step") != std::string::npos) {
            config.num_time_step = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("time_step_length") != std::string::npos) {
            config.time_step_length = std::stod(line.substr(line.find("=") + 1));
        } else if (line.find("LSYS") != std::string::npos) {
            config.LSYS = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("K") != std::string::npos) {
            config.K = std::stod(line.substr(line.find("=") + 1));
        } else if (line.find("gamma") != std::string::npos) {
            config.gamma = std::stod(line.substr(line.find("=") + 1));
        } else if (line.find("cut_off_radius") != std::string::npos) {
            config.cut_off_radius = std::stod(line.substr(line.find("=") + 1));
        } else if (line.find("particle_num") != std::string::npos) {
            config.particle_num = std::stoi(line.substr(line.find("=") + 1));
        } else if (line.find("gravity") != std::string::npos) {
            config.gravity = std::stod(line.substr(line.find("=") + 1));
        } else if (line.find("particle_positions") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, vector_regex)) {
                    positions.push_back({
                        std::stod(match[1]),
                        std::stod(match[2]),
                        std::stod(match[3])
                    });
                }
            }
        } else if (line.find("particle_velocities") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, vector_regex)) {
                    velocities.push_back({
                        std::stod(match[1]),
                        std::stod(match[2]),
                        std::stod(match[3])
                    });
                }
            }
        } else if (line.find("particle_masses") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, number_regex)) {
                    masses.push_back(std::stod(match[0]));
                }
            }
        }
        else if (line.find("particle_radiuses") != std::string::npos) {
            while (std::getline(infile, line) && line.find("};") == std::string::npos) {
                std::smatch match;
                if (std::regex_search(line, match, number_regex)) {
                    radiuses.push_back(std::stod(match[0]));
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
        p.radius = radiuses[i];
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
       << "  Mass: " << particle.mass << "\n"
       << " Radius: " << particle.radius;
    return os;
}

std::ostream& operator<<(std::ostream& os, const Config& config) {
    os << "Time Steps: " << config.num_time_step << "\n";
    os << "Time Step Length: " << config.time_step_length << "\n";
    os << "K: " << config.K << "\n";
    os << "gamma: " << config.gamma << "\n";
    os << "LSYS: " << config.LSYS << "\n";
    os << "cut off radius: " << config.cut_off_radius << "\n";
    os << "gravity: " << config.gravity << "\n";
    os << "Particle Count: " << config.particle_num << "\n";


    for (size_t i = 0; i < config.particles.size(); ++i) {
        os << "Particle " << i + 1 << ":\n" << config.particles[i] << "\n";
    }

    return os;
}
