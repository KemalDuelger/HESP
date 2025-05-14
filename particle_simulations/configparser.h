#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <string>
#include <vector>
#include <ostream>

struct Vector3 {
    float x = 0;
    float y = 0;
    float z = 0;

    Vector3(float x,float y,float z) : x(x), y(y), z(z) {}
};

struct Particle {
    Vector3 position;
    Vector3 velocity;
    float mass;
    Vector3 force;
};

struct Config {
    int num_time_step = 0;
    int time_step_length = 0;
    float sigma = 0.0f;
    float epsilon = 0.0f;
    int particle_num = 0;
    std::vector<Particle> particles;
};

class ConfigParser {
public:
    static Config parse(const std::string& filename);
};

// Stream-Ausgaben
std::ostream& operator<<(std::ostream& os, const Vector3& vec);
std::ostream& operator<<(std::ostream& os, const Particle& particle);
std::ostream& operator<<(std::ostream& os, const Config& config);

#endif // CONFIGPARSER_H
