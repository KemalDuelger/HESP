#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <cuda_runtime.h>
#include "configparser.h"
#include <chrono>


// Cuda-Kernel für updateNeighborList
__global__ void buildCellListKernel(Particle* particles, int n, int* cell_head, int* cell_next, 
                                    double cutoff_radius, int LSYS, int num_cells_per_dim) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    Particle p = particles[i];
    int ix = (int)(p.position.x / cutoff_radius);
    int iy = (int)(p.position.y / cutoff_radius);
    int iz = (int)(p.position.z / cutoff_radius);

    // Periodische Bedingungen absichern (+ num_cells_per_dim um negative Indizes zu vermeiden)
    ix = (ix + num_cells_per_dim) % num_cells_per_dim;
    iy = (iy + num_cells_per_dim) % num_cells_per_dim;
    iz = (iz + num_cells_per_dim) % num_cells_per_dim;

    int cell_idx = ix + iy * num_cells_per_dim + iz * num_cells_per_dim * num_cells_per_dim;

    // Atomare Verkettung in Linked List
    cell_next[i] = atomicExch(&cell_head[cell_idx], i);
}

// CUDA-Kernel für computeForces
__global__ void computeForcesCellKernel(Particle* particles, int n, double epsilon, double sigma, 
                                        int LSYS, double cutoff_radius,
                                        int* cell_head, int* cell_next,
                                        int num_cells_per_dim) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;

    Particle& pi = particles[i];
    pi.force.x = 0.0;
    pi.force.y = 0.0;   
    pi.force.z = 0.0;

    int ix = (int)(pi.position.x / cutoff_radius);
    int iy = (int)(pi.position.y / cutoff_radius);
    int iz = (int)(pi.position.z / cutoff_radius);

    for (int dx = -1; dx <= 1; ++dx) {
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dz = -1; dz <= 1; ++dz) {
                int nx = (ix + dx + num_cells_per_dim) % num_cells_per_dim;
                int ny = (iy + dy + num_cells_per_dim) % num_cells_per_dim;
                int nz = (iz + dz + num_cells_per_dim) % num_cells_per_dim;

                int neighbor_cell = nx + ny * num_cells_per_dim + nz * num_cells_per_dim * num_cells_per_dim;

                for (int j = cell_head[neighbor_cell]; j != -1; j = cell_next[j]) {
                    if (i == j) continue;

                    Particle pj = particles[j];
                    double dx = pj.position.x - pi.position.x;
                    double dy = pj.position.y - pi.position.y;
                    double dz = pj.position.z - pi.position.z;

                    dx -= LSYS * round(dx / LSYS);
                    dy -= LSYS * round(dy / LSYS);
                    dz -= LSYS * round(dz / LSYS);

                    double r2 = dx * dx + dy * dy + dz * dz;
                    if (r2 > cutoff_radius * cutoff_radius) continue;

                    double r = sqrt(r2);
                    double s_over_r = sigma / r;
                    double s_over_r6 = pow(s_over_r, 6);
                    double lj_scalar = (24 * epsilon * s_over_r6 * (2 * s_over_r6 - 1)) / r2;

                    pi.force.x -= lj_scalar * dx;
                    pi.force.y -= lj_scalar * dy;
                    pi.force.z -= lj_scalar * dz;
                }
            }
        }
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
    
    // Zeitmessung starten
    auto start = std::chrono::high_resolution_clock::now();
    
    // Linked List für Verlet Neighbor List
    int num_cells_per_dim = std::ceil(config.LSYS / config.cut_off_radius);
    int total_cells = num_cells_per_dim * num_cells_per_dim * num_cells_per_dim;

    std::vector<int> h_cell_head(total_cells, -1); // cell → erster Partikel
    std::vector<int> h_cell_next(n, -1);           // partikel → nächster in Zelle

    int* d_cell_head;
    int* d_cell_next;
    cudaMalloc(&d_cell_head, total_cells * sizeof(int));
    cudaMalloc(&d_cell_next, n * sizeof(int));

    // Speicher auf GPU anlegen und kopieren
    Particle* d_particles;
    cudaMalloc(&d_particles, n * sizeof(Particle));
    cudaMemcpy(d_particles, particles.data(), n * sizeof(Particle), cudaMemcpyHostToDevice);

    int threads = 128;
    int blocks = (n + threads - 1) / threads;
    
    cudaError_t err;

    for (int t = 0; t < config.num_time_step; ++t) {

        // Zellzuordnung aktualisieren (jedes Mal nötig)
        cudaMemset(d_cell_head, -1, total_cells * sizeof(int));
        buildCellListKernel<<<blocks, threads>>>(d_particles, n, d_cell_head, d_cell_next,
                                             config.cut_off_radius, config.LSYS, num_cells_per_dim);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel Error: " << cudaGetErrorString(err) << std::endl;
        }

        updatePositionAndHalfStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length, config.LSYS);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel Error: " << cudaGetErrorString(err) << std::endl;
        }
        
        computeForcesCellKernel<<<blocks, threads>>>(d_particles, n, config.epsilon, config.sigma, config.LSYS,
                                                 config.cut_off_radius, d_cell_head, d_cell_next, num_cells_per_dim);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel Error: " << cudaGetErrorString(err) << std::endl;
        }

        updateAccelerationAndFullStepVelocityKernel<<<blocks, threads>>>(d_particles, n, config.time_step_length);
        cudaDeviceSynchronize();
        err = cudaGetLastError();
        if (err != cudaSuccess) {
            std::cerr << "Kernel Error: " << cudaGetErrorString(err) << std::endl;
        }
        // Ausgabe alle 100 Schritte
        if (t % 10 == 0) {
            cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
            writeToVTK(particles, t, config.LSYS);
        }
    }

    cudaMemcpy(particles.data(), d_particles, n * sizeof(Particle), cudaMemcpyDeviceToHost);
    cudaFree(d_particles);

        // Zeitmessung beenden
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Ausführungszeit: " << diff.count() << " Sekunden\n";

    return 0;
}