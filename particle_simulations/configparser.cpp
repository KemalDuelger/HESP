#include "configparser.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>

// Helper function to parse vector data
std::vector<std::vector<float>> parseVectorList(const std::string& value) {
    std::vector<std::vector<float>> vectors;
    std::regex vectorRegex(R"(\(([^)]+)\))");
    std::smatch match;

    std::string::const_iterator searchStart(value.cbegin());
    while (std::regex_search(searchStart, value.cend(), match, vectorRegex)) {
        std::vector<float> vec;
        std::istringstream vecStream(match[1].str());
        std::string component;
        while (std::getline(vecStream, component, ',')) {
            vec.push_back(std::stof(component));
        }
        vectors.push_back(vec);
        searchStart = match.suffix().first;
    }

    return vectors;
}

// Function to parse the configuration file
void parseConfigFile(const std::string& filename, 
                     std::map<std::string, std::string>& config, 
                     std::vector<std::vector<float>>& particlePositions, 
                     std::vector<std::vector<float>>& particleVelocities) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        // Remove spaces
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        // Skip empty lines or comments
        if (line.empty() || line[0] == '/' || line[0] == '#') {
            continue;
        }

        // Parse key-value pairs
        std::size_t pos = line.find('=');
        if (pos != std::string::npos) {
            std::string key = line.substr(0, pos);
            std::string value = line.substr(pos + 1);

            // Remove trailing semicolon
            if (!value.empty() && value.back() == ';') {
                value.pop_back();
            }

            config[key] = value;
        }
    }

    inputFile.close();

    // Parse particle positions and velocities
    if (config.find("particle_positions") != config.end()) {
        particlePositions = parseVectorList(config["particle_positions"]);
    }
    if (config.find("particle_velocities") != config.end()) {
        particleVelocities = parseVectorList(config["particle_velocities"]);
    }
}