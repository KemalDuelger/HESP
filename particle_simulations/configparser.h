#pragma once

#include <string>
#include <map>
#include <vector>

// Function to parse the configuration file
void parseConfigFile(const std::string& filename, 
                     std::map<std::string, std::string>& config, 
                     std::vector<std::vector<float>>& particlePositions, 
                     std::vector<std::vector<float>>& particleVelocities);