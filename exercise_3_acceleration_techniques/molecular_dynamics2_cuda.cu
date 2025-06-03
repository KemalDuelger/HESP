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

    Particle& pi = particles[i];
    pi.force.x = 0.0;
    pi.force.y = 0.0;
    pi.force.z = 0.0;

    for (int j = 0; j < n; ++j) {
        if (i == j) continue;

        Particle& pj = particles[j];
        double dx = pj.position.x - pi.position.x;
        double dy = pj.position.y - pi.position.y;
        double dz = pj.position.z - pi.position.z;
        // Minimum-Image-Kriterium
        dx -= LSYS * round(dx / LSYS);
        dy -= LSYS * round(dy / LSYS);
        dz -= LSYS * round(dz / LSYS);

        double r2 = dx * dx + dy * dy + dz * dz;
        double r = sqrt(r2);

        double s_over_r = sigma / r;
        double s_over_r6 = pow(s_over_r, 6);
        double lj_scalar = (24 * epsilon * s_over_r6 * (2 * s_over_r6 - 1)) / r2;

        pi.force.x += lj_scalar * dx;
        pi.force.y += lj_scalar * dy;
        pi.force.z += lj_scalar * dz;
    }

    // Optional: Schwerkraft in -y-Richtung hinzufügen
    // pi.force.y -= 9.81 * pi.mass;
}


// Cuda-Methode fuer LSYS Berechnungen
__device__ double wrapPosition(double pos, int limit) {
    if (pos >= limit) return pos - limit;
    if (pos < 0.0)    return pos + limit;
    return pos;
}

// CUDA-Kernel für updatePositionAndHalfStepVelocity
__global__ void updatePositionAndHalfStepVelocityKernel(Particle* particles, int n, double dt, int LSYS) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    Particle& p = particles[i];

    double dt2 = 0.5 * dt * dt;

    // Position update
    p.position.x += p.velocity.x * dt + p.acceleration.x * dt2;
    p.position.y += p.velocity.y * dt + p.acceleration.y * dt2;
    p.position.z += p.velocity.z * dt + p.acceleration.z * dt2;

    // Apply periodic boundary conditions
    p.position.x = wrapPosition(p.position.x, LSYS);
    p.position.y = wrapPosition(p.position.y, LSYS);
    p.position.z = wrapPosition(p.position.z, LSYS);

    // Half-step velocity update
    double half_dt = 0.5 * dt;
    p.velocity.x += p.acceleration.x * half_dt;
    p.velocity.y += p.acceleration.y * half_dt;
    p.velocity.z += p.acceleration.z * half_dt;
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
        std::cerr << "Error opening file: " << filename << std::endl;
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
    std::cout << config << std::endl;

    // Speicher auf GPU anlegen und kopieren
    Particle* d_particles;
    cudaMalloc(&d_particles, n * sizeof(Particle));
    cudaMemcpy(d_particles, particles.data(), n * sizeof(Particle), cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n + threads - 1) / threads;
    

    for (int t = 0; t < config.num_time_step; ++t) {
        updatePositionAndHalfStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length, config.LSYS);
        cudaDeviceSynchronize();

        computeForcesKernel<<<blocks, threads>>>(d_particles, n, config.epsilon, config.sigma, config.LSYS);
        cudaDeviceSynchronize();

        updateAccelerationAndFullStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length);
        cudaDeviceSynchronize();

        // Ausgabe alle 100 Schritte
        if (t % 10 == 0) {
            cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
            writeToVTK(particles, t);
        }
    }

    cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
    cudaFree(d_particles);

    return 0;
}