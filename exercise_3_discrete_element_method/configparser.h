#ifndef CONFIGPARSER_H
#define CONFIGPARSER_H

#include <string>
#include <vector>
#include <ostream>

struct Vector3 {
    double x;
    double y;
    double z;

    Vector3(double x = 0,double y = 0,double z = 0) : x(x), y(y), z(z) {}
};

struct Particle {
    Vector3 position;
    Vector3 velocity;
    double mass;
    Vector3 force;
    Vector3 acceleration;
    double radius;
};

struct Config {
    int num_time_step = 0;
    double time_step_length = 0;
    double K = 0.0;
    double gamma = 0.0;
    int particle_num = 0.0;
    double cut_off_radius = 0.0;
    int LSYS = 0;
    double gravity = 0.0;
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
