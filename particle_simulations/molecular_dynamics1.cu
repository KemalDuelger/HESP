#include <iostream>
#include <map>
#include <vector>
#include "configparser.h"

void computeForces(std::vector<Particle> particles, Config config) {
    for( auto particle : particles) {
        particle.force.y = 9.81*particle.mass;

    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <configfile>" << std::endl;
        return 1;
    }

    Config config = ConfigParser::parse(argv[1]);
    
    std::vector<Particle> particles = config.particles;

    const int numSteps = 10000;       // Anzahl der Zeitschritte

    for (int t = 0; t < numSteps; ++t) {
        // 1. Kraftberechnung (basierend auf Position)
        computeForces(particles, config);

        // 2. Positionen aktualisieren (Velocity-Verlet Schritt 1)
        updatePositions(particles, config.time_step_length);

        // 3. Neue Kräfte berechnen (nach neuen Positionen)
        computeForces(particles, config);

        // 4. Geschwindigkeiten aktualisieren (Velocity-Verlet Schritt 2)
        updateVelocities(particles, config.time_step_length);

        // 5. Optional: Ausgabe alle n Schritte
        if (t % 100 == 0) {
            writeToVTK(particles, t);
        }
    }

    return 0;
}