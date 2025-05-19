#include <iostream>
#include <map>
#include <vector>
#include "configparser.h"
#include <cmath>
#include <fstream>

void computeForces(std::vector<Particle> particles, Config config, int LSYS) {
    // Berechnung der Kräfte für alle Partikel
    for (size_t i = 0; i < particles.size(); ++i) {
        particles[i].force = {0, 0, 0};

        for (size_t j = 0; j < particles.size(); ++j) {
            if (i == j) continue;

            double dx = particles[j].position.x - particles[i].position.x;
            double dy = particles[j].position.y - particles[i].position.y;
            double dz = particles[j].position.z - particles[i].position.z;

            // Minimum-Image-Kriterium anwenden
            dx -= LSYS * round(dx / LSYS);
            dy -= LSYS * round(dy / LSYS);
            dz -= LSYS * round(dz / LSYS);

            double abstand = sqrt(dx*dx + dy*dy + dz*dz);

            // Hier Kraft berechnen, z.B. Lennard-Jones, etc.
            particles[i].force.x = 24 * config.epsilon * pow(config.sigma / abstand, 6) *(2 * pow(config.sigma / abstand, 6) -1) * dx / (abstand*abstand);
            particles[i].force.y = 24 * config.epsilon * pow(config.sigma / abstand, 6) *(2 * pow(config.sigma / abstand, 6) -1) * dy / (abstand*abstand);
            particles[i].force.z = 24 * config.epsilon * pow(config.sigma / abstand, 6) *(2 * pow(config.sigma / abstand, 6) -1) * dz / (abstand*abstand);

        }

        particles[i].force.y += 9.81 * particles[i].mass;
    }
}

void updatePositionAndHalfStepVelocity(std::vector<Particle>& particles, float time_step_length) {
    for (size_t i = 0; i < particles.size(); ++i) {
        // xi (t+ dt) rechnen
        particles[i].position.x += particles[i].velocity.x * time_step_length + 0.5 * particles[i].acceleration.x * time_step_length * time_step_length;
        particles[i].position.y += particles[i].velocity.y * time_step_length + 0.5 * particles[i].acceleration.y * time_step_length * time_step_length;    
        particles[i].position.z += particles[i].velocity.z * time_step_length + 0.5 * particles[i].acceleration.z * time_step_length * time_step_length;    
        
        // vi (t +1/2 dt) rechnen
        particles[i].velocity.x += 0.5 * particles[i].acceleration.x * time_step_length;
        particles[i].velocity.y += 0.5 * particles[i].acceleration.y * time_step_length;
        particles[i].velocity.z += 0.5 * particles[i].acceleration.z * time_step_length;

    }   

}

void updateAccelerationAndFullStepVelocity(std::vector<Particle>& particles, float time_step_length) {
    for (size_t i = 0; i < particles.size(); ++i) {
        // ai (t + dt) rechnen
        particles[i].acceleration.x = particles[i].force.x / particles[i].mass;
        particles[i].acceleration.y = particles[i].force.y / particles[i].mass;
        particles[i].acceleration.z = particles[i].force.z / particles[i].mass;
        // vi (t + dt) rechnen
        particles[i].velocity.x += 0.5 * particles[i].acceleration.x * time_step_length;
        particles[i].velocity.y += 0.5 * particles[i].acceleration.y * time_step_length;
        particles[i].velocity.z += 0.5 * particles[i].acceleration.z * time_step_length;
    }
}

void writeToVTK(const std::vector<Particle>& particles, int step) {
    std::string filename = "output_" + std::to_string(step) + ".vtk";
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Fehler beim Öffnen der Datei: " << filename << std::endl;
        return;
    }

    ofs << "# vtk DataFile Version 3.0\n";
    ofs << "Particle data\n";
    ofs << "ASCII\n";
    ofs << "DATASET POLYDATA\n";
    ofs << "POINTS " << particles.size() << " float\n";
    for (const auto& p : particles) {
        ofs << p.position.x << " " << p.position.y << " " << p.position.z << "\n";
    }
    ofs << "VERTICES " << particles.size() << " " << particles.size() * 2 << "\n";
    for (size_t i = 0; i < particles.size(); ++i) {
        ofs << "1 " << i << "\n";
    }
    ofs.close();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configfile>" << std::endl;
        return 1;
    }

    Config config = ConfigParser::parse(argv[1]);
    
    std::vector<Particle> particles = config.particles;

    const int numSteps = 10000;       // Anzahl der Zeitschritte
    const int LSYS = 100;            // Systemgröße (z.B. Boxgröße)
    for (int t = 0; t < numSteps; ++t) {
        // 1. Positionen aktualisieren (Velocity-Verlet Schritt 1)
        updatePositionAndHalfStepVelocity(particles, config.time_step_length);

        // 2. Neue Kräfte berechnen
        computeForces(particles, config, LSYS);

        // 3. Geschwindigkeiten aktualisieren (Velocity-Verlet Schritt 2)
        updateAccelerationAndFullStepVelocity(particles, config.time_step_length);

        // 4. Optional: Ausgabe alle n Schritte
        if (t % 100 == 0) {
            writeToVTK(particles, t);
        }
    }

    return 0;
}