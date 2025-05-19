#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cuda_runtime.h>
#include "configparser.h"


// CUDA-Kernel für computeForces
__global__ void computeForcesKernel(Particle* particles, int n, double epsilon, double sigma, int LSYS) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    for (int j = 0; j < n; ++j) {
        if (i == j) continue;
        double dx = particles[j].position.x - particles[i].position.x;
        double dy = particles[j].position.y - particles[i].position.y;
        double dz = particles[j].position.z - particles[i].position.z;
        // Minimum-Image-Kriterium
        dx -= LSYS * round(dx / LSYS);
        dy -= LSYS * round(dy / LSYS);
        dz -= LSYS * round(dz / LSYS);
        double abstand = sqrt(dx*dx + dy*dy + dz*dz);
        if (abstand == 0) continue;
        double lj = 24 * epsilon * pow(sigma / abstand, 6) * (2 * pow(sigma / abstand, 6) - 1) / (abstand*abstand);
        particles[i].force.x += lj * dx;
        particles[i].force.y += lj * dy;
        particles[i].force.z += lj * dz;
    }
    particles[i].force.y -= 9.81 * particles[i].mass;
}

// CUDA-Kernel für updatePositionAndHalfStepVelocity
__global__ void updatePositionAndHalfStepVelocityKernel(Particle* particles, int n, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    particles[i].position.x += particles[i].velocity.x * dt + 0.5 * particles[i].acceleration.x * dt * dt;
    particles[i].position.y += particles[i].velocity.y * dt + 0.5 * particles[i].acceleration.y * dt * dt;
    particles[i].position.z += particles[i].velocity.z * dt + 0.5 * particles[i].acceleration.z * dt * dt;

    particles[i].velocity.x += 0.5 * particles[i].acceleration.x * dt;
    particles[i].velocity.y += 0.5 * particles[i].acceleration.y * dt;
    particles[i].velocity.z += 0.5 * particles[i].acceleration.z * dt;
}

// CUDA-Kernel für updateAccelerationAndFullStepVelocity
__global__ void updateAccelerationAndFullStepVelocityKernel(Particle* particles, int n, double dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    particles[i].acceleration.x = particles[i].force.x / particles[i].mass;
    particles[i].acceleration.y = particles[i].force.y / particles[i].mass;
    particles[i].acceleration.z = particles[i].force.z / particles[i].mass;

    particles[i].velocity.x += 0.5 * particles[i].acceleration.x * dt;
    particles[i].velocity.y += 0.5 * particles[i].acceleration.y * dt;
    particles[i].velocity.z += 0.5 * particles[i].acceleration.z * dt;
}

// Host-Funktion für VTK-Ausgabe (wie gehabt)
void writeToVTK(const std::vector<Particle>& particles, int step) {
    std::string filename = "vtk_files/output_" + std::to_string(step) + ".vtk";
    std::ofstream ofs(filename);
    if (!ofs) {
        std::cerr << "Fehler beim Öffnen der Datei: " << filename << std::endl;
        return;
    }

    ofs << "# vtk DataFile Version 4.0\n";
    ofs << "hesp visualization file\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    ofs << "POINTS " << particles.size() << " double\n";
    for (const auto& p : particles) {
        ofs << p.position.x << " " << p.position.y << " " << p.position.z << "\n";
    }
    ofs << "CELLS 0 0\nCELL_TYPES 0\n";
    ofs << "POINT_DATA " << particles.size() << "\n";
    ofs << "SCALARS m double\nLOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        ofs << p.mass << "\n";
    }
    ofs << "VECTORS v double\n";
    for (const auto& p : particles) {
        ofs << p.velocity.x << " " << p.velocity.y << " " << p.velocity.z << "\n";
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
    int n = particles.size();
    const int LSYS = 100;

    // Speicher auf GPU anlegen und kopieren
    Particle* d_particles;
    cudaMalloc(&d_particles, n * sizeof(Particle));
    cudaMemcpy(d_particles, particles.data(), n * sizeof(Particle), cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n + threads - 1) / threads;

    for (int t = 0; t < config.num_time_step; ++t) {
        updatePositionAndHalfStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length);
        cudaDeviceSynchronize();

        computeForcesKernel<<<blocks, threads>>>(d_particles, n, config.epsilon, config.sigma, LSYS);
        cudaDeviceSynchronize();

        updateAccelerationAndFullStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length);
        cudaDeviceSynchronize();

        // Ausgabe alle 100 Schritte
        if (t % 100 == 0) {
            cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
            writeToVTK(particles, t);
        }
    }

    cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
    cudaFree(d_particles);

    return 0;
}