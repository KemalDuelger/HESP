#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cuda_runtime.h>
#include "configparser.h"

void buildNeighbourList(const std::vector<Particle>& particles, double r_list, int LSYS,
                        std::vector<int>& neighbourList,
                        std::vector<int>& neighbourListStarts,
                        std::vector<int>& neighbourListLengths) {
    neighbourList.clear();
    int n = particles.size();
    for (int i = 0; i < n; ++i) {
        neighbourListStarts[i] = neighbourList.size();
        int count = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double dx = particles[j].position.x - particles[i].position.x;
            double dy = particles[j].position.y - particles[i].position.y;
            double dz = particles[j].position.z - particles[i].position.z;
            dx -= LSYS * round(dx / LSYS);
            dy -= LSYS * round(dy / LSYS);
            dz -= LSYS * round(dz / LSYS);

            double r2 = dx * dx + dy * dy + dz * dz;
            if (r2 < r_list * r_list) {
                neighbourList.push_back(j);
                ++count;
            }
        }
        neighbourListLengths[i] = count;
    }
}


// CUDA-Kernel für computeForces
__global__ void computeForcesNeighbourListKernel(Particle* particles, int n,
    double epsilon, double sigma, int LSYS, double cut_off_radius,
    const int* neighbourList, const int* neighbourListStarts, const int* neighbourListLengths) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    Particle& pi = particles[i];
    pi.force = {0.0, 0.0, 0.0};

    int start = neighbourListStarts[i];
    int length = neighbourListLengths[i];

    for (int k = 0; k < length; ++k) {
        int j = neighbourList[start + k];
        Particle& pj = particles[j];

        double dx = pj.position.x - pi.position.x;
        double dy = pj.position.y - pi.position.y;
        double dz = pj.position.z - pi.position.z;

        dx -= LSYS * round(dx / LSYS);
        dy -= LSYS * round(dy / LSYS);
        dz -= LSYS * round(dz / LSYS);

        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > cut_off_radius * cut_off_radius) continue;

        double r = sqrt(r2);
        double s_over_r = sigma / r;
        double s_over_r6 = pow(s_over_r, 6);
        double lj_scalar = (24 * epsilon * s_over_r6 * (2 * s_over_r6 - 1)) / r2;

        pi.force.x -= lj_scalar * dx;
        pi.force.y -= lj_scalar * dy;
        pi.force.z -= lj_scalar * dz;
    }
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
void writeToVTK(const std::vector<Particle>& particles, int step, int LSYS) {
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

    // Neighbour List Speicher vorbereiten (Host)
    std::vector<int> neighbourList;
    std::vector<int> neighbourListStarts(config.particle_num);
    std::vector<int> neighbourListLengths(config.particle_num);


    // Speicher auf GPU anlegen und kopieren
    Particle* d_particles;
    cudaMalloc(&d_particles, n * sizeof(Particle));
    cudaMemcpy(d_particles, particles.data(), n * sizeof(Particle), cudaMemcpyHostToDevice);

    int* d_neighbourList;
    int* d_neighbourListStarts;
    int* d_neighbourListLengths;

    cudaMalloc(&d_neighbourList, config.particle_num * config.particle_num * sizeof(int)); // worst case
    cudaMalloc(&d_neighbourListStarts, config.particle_num * sizeof(int));
    cudaMalloc(&d_neighbourListLengths, config.particle_num * sizeof(int));

    int threads = 128;
    int blocks = (n + threads - 1) / threads;
    

    for (int t = 0; t < config.num_time_step; ++t) {
        updatePositionAndHalfStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length, config.LSYS);
        cudaDeviceSynchronize();

        if (t % 10 == 0) {
            cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
            
            buildNeighbourList(particles, config.cut_off_radius * 1.1, config.LSYS,
                            neighbourList, neighbourListStarts, neighbourListLengths);

            cudaMemcpy(d_neighbourList, neighbourList.data(), neighbourList.size() * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_neighbourListStarts, neighbourListStarts.data(), n * sizeof(int), cudaMemcpyHostToDevice);
            cudaMemcpy(d_neighbourListLengths, neighbourListLengths.data(), n * sizeof(int), cudaMemcpyHostToDevice);
        }

        computeForcesNeighbourListKernel<<<blocks, threads>>>(d_particles, n, config.epsilon, config.sigma, config.LSYS,
                                                            config.cut_off_radius, d_neighbourList,
                                                            d_neighbourListStarts, d_neighbourListLengths);
        cudaDeviceSynchronize();

        updateAccelerationAndFullStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length);
        cudaDeviceSynchronize();

        // Ausgabe alle 100 Schritte
        if (t % 10 == 0) {
            cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
            writeToVTK(particles, t, config.LSYS);
        }
    }

    cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
    cudaFree(d_particles);
    cudaFree(d_neighbourList);
    cudaFree(d_neighbourListStarts);
    cudaFree(d_neighbourListLengths);

    return 0;
}