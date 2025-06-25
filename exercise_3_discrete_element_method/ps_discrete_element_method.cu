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

    // Periodische Bedingungen absichern
    ix = (ix + num_cells_per_dim) % num_cells_per_dim;
    iy = (iy + num_cells_per_dim) % num_cells_per_dim;
    iz = (iz + num_cells_per_dim) % num_cells_per_dim;

    int cell_idx = ix + iy * num_cells_per_dim + iz * num_cells_per_dim * num_cells_per_dim;

    // Atomare Verkettung in Linked List
    cell_next[i] = atomicExch(&cell_head[cell_idx], i);
}

__global__ void computeForcesCellKernel(
    Particle* particles, int n, double gamma, double K, int LSYS, double cutoff_radius,
    int* cell_head, int* cell_next, int num_cells_per_dim, double gravity) {

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i >= n) return;

    Particle& pi = particles[i];

    pi.force.x = 0.0;
    pi.force.y = 0.0;   
    pi.force.z = 0.0;

    int ix = (int)(pi.position.x / cutoff_radius);
    int iy = (int)(pi.position.y / cutoff_radius);
    int iz = (int)(pi.position.z / cutoff_radius);



    // Nachbarsuche (keine Periodizität mehr!)
   for (int dx = -1; dx <= 1; ++dx) {

        int nx = ix + dx;

        if (nx < 0 || nx >= num_cells_per_dim) continue;

        for (int dy = -1; dy <= 1; ++dy) {

            int ny = iy + dy;

            if (ny < 0 || ny >= num_cells_per_dim) continue;

            for (int dz = -1; dz <= 1; ++dz) {

                int nz = iz + dz;

                if (nz < 0 || nz >= num_cells_per_dim) continue; The resultant expression for the force is quite simple:

                int neighbor_cell = nx + ny * num_cells_per_dim + nz * num_cells_per_dim * num_cells_per_dim;

                for (int j = cell_head[neighbor_cell]; j != -1; j = cell_next[j]) {

                    if (i == j) continue;

                    Particle pj = particles[j];

                    double dx = pi.position.x - pj.position.x;
                    double dy = pi.position.y - pj.position.y;
                    double dz = pi.position.z - pj.position.z;

                    double sphere_diameter = 2*pi.radius;
                    double xij = sqrt(dx*dx + dy*dy + dz*dz);
                    double overlap_cond = sphere_diameter - xij;
                    double epsilon_x = 0;
                    double epsilon_y = 0;
                    double epsilon_z = 0;
                    double nx = dx / xij;
                    double ny = dy / xij;
                    double nz = dz / xij;


                    if (overlap_cond >= 0.0) {
                        epsilon_x = nx * overlap_cond;
                        epsilon_y = ny * overlap_cond;
                        epsilon_z = nz * overlap_cond;
                    }

                    double epsilon_der_x = 0;
                    double epsilon_der_y = 0;
                    double epsilon_der_z = 0;

                    if (sphere_diameter - fabs(dx) >= 0)
                    {
                        epsilon_der_x = -nx * (nx * (pi.velocity.x - pj.velocity.x));
                    }

                    if (sphere_diameter - fabs(dy) >= 0)
                    {
                        epsilon_der_y = -ny * (ny * (pi.velocity.y - pj.velocity.y));
                    }

                    if (sphere_diameter - fabs(dz) >= 0)
                    {
                        epsilon_der_z = -nz * (nz * (pi.velocity.z - pj.velocity.z));
                    }

                    pi.force.x += K* epsilon_x + gamma * epsilon_der_x;
                    pi.force.y += K* epsilon_y + gamma * epsilon_der_y; 
                    pi.force.z += K* epsilon_z + gamma * epsilon_der_z;

                    /*

                    if (overlap >= 0.0) { // Kontakt!

                        // Einheitsvektor

                        double nx = dx / dist;

                        double ny = dy / dist;

                        double nz = dz / dist;



                        // Relativgeschwindigkeit in Richtung des Kontakts

                        double dvx = pi.velocity.x - pj.velocity.x;

                        double dvy = pi.velocity.y - pj.velocity.y;

                        double dvz = pi.velocity.z - pj.velocity.z;

                        double v_rel = dvx * nx + dvy * ny + dvz * nz;



                        // Federkraft + Dämpfung

                        double f_scalar = K * overlap + gamma * v_rel;



                        pi.force.x += f_scalar * nx;

                        pi.force.y += f_scalar * ny;

                        pi.force.z += f_scalar * nz;

                    }

                    */

                }

            }

        }

    }



    double pos[3] = {pi.position.x, pi.position.y, pi.position.z};
    double vel[3] = {pi.velocity.x, pi.velocity.y, pi.velocity.z};
    double* force[3] = {&pi.force.x, &pi.force.y, &pi.force.z};


    for (int d = 0; d < 3; ++d) {
        double overlap = pi.radius - pos[d];
        if (overlap > 0.0) {
            double v_rel = -vel[d];
            double f = K * overlap + gamma * v_rel;
            *force[d] += f;
        }

        overlap = pi.radius - (LSYS - pos[d]);

        if (overlap > 0.0) {
            double v_rel = -vel[d];
            double f = K * overlap + gamma * v_rel;
            *force[d] -= f;
        }

    }

    // Gravitation (z.B. in -y Richtung)

    pi.force.y += gravity * pi.mass;

}

    // --- Wandkontakte (Box: 0 <= x,y,z <= LSYS) ---
    double r = pi.radius;


    // Linke Wand (x = 0)
    double overlap = r - pi.position.x;

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.x;
        double f = K * overlap + gamma * v_rel;
        pi.force.x += f;

    }

    // Rechte Wand (x = LSYS)
    overlap = r - (LSYS - pi.position.x);

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.x;
        double f = K * overlap + gamma * v_rel;
        pi.force.x -= f;

    }

    // Untere Wand (y = 0)
    overlap = r - pi.position.y;

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.y;
        double f = K * overlap + gamma * v_rel;
        pi.force.y += f;

    }

    // Obere Wand (y = LSYS)
    overlap = r - (LSYS - pi.position.y);

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.y;
        double f = K * overlap + gamma * v_rel;
        pi.force.y -= f;

    }

    // Vorderwand (z = 0)
    overlap = r - pi.position.z;

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.z;
        double f = K * overlap + gamma * v_rel;
        pi.force.z += f;
    }

    // Rückwand (z = LSYS)
    overlap = r - (LSYS - pi.position.z);

    if (overlap > 0.0) {
        double v_rel = -pi.velocity.z;
        double f = K * overlap + gamma * v_rel;
        pi.force.z -= f;

    }

    // Gravitation (z.B. in -y Richtung)
    pi.force.y += gravity * pi.mass;

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

    // Periodische Randbedingungen
    // Sicherstellen, dass Partikel innerhalb des Systems bleiben
    p.position.x = fmax(p.radius, fmin(LSYS - p.radius, p.position.x));
    p.position.y = fmax(p.radius, fmin(LSYS - p.radius, p.position.y));
    p.position.z = fmax(p.radius, fmin(LSYS - p.radius, p.position.z));


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
    ofs << "SCALARS radius double\nLOOKUP_TABLE default\n";
    for (const auto& p : particles) {
        ofs << p.radius << "\n";
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
        
        computeForcesCellKernel<<<blocks, threads>>>(d_particles, n, config.gamma, config.K, config.LSYS,
                                                    config.cut_off_radius, d_cell_head, d_cell_next, num_cells_per_dim,
                                                    config.gravity);
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