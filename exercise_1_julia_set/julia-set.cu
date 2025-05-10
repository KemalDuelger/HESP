#include <iostream>
#include <complex>
#include <vector>
#include <cuda_runtime.h>
#include "lodepng.h" // Für das Speichern als PNG

// CUDA-Kernel zur Berechnung der Julia-Menge
__global__ void julia_kernel(int width, int height, double x_min, double x_max, double y_min, double y_max, 
                             double escape_radius, int max_iter, double c_real, double c_imag, int* image) {
    int px = blockIdx.x * blockDim.x + threadIdx.x; // Pixel x-Koordinate
    int py = blockIdx.y * blockDim.y + threadIdx.y; // Pixel y-Koordinate

    if (px < width && py < height) {
        // Skalierung der Pixelkoordinaten auf den komplexen Zahlenraum y [-2, 2] x [-2, 2] nach Vorgabe
        double x = x_min + (x_max - x_min) * px / (width - 1);
        double y = y_min + (y_max - y_min) * py / (height - 1);
        double z_real = x;
        double z_imag = y;

        // Iteration der Julia-Funktion
        // Kernel Methode kennt keine std::complex bzw std Lib, daher manuelle Berechnung
        int iter = 0;
        while (z_real * z_real + z_imag * z_imag < escape_radius * escape_radius && iter < max_iter) {
            double temp = z_real * z_real - z_imag * z_imag + c_real;
            z_imag = 2.0 * z_real * z_imag + c_imag;
            z_real = temp;
            ++iter;
        }

        // Speichern der Iterationsanzahl im Bild (Graustufen-Mapping)
        image[py * width + px] = static_cast<int>(255.0 * iter / max_iter); // Skaliert auf 0-255
    }
}

int main() {
    // Konstanten für die Julia-Menge
    const int width = 2048;
    const int height = 2048;
    const double x_min = -2.0, x_max = 2.0;
    const double y_min = -2.0, y_max = 2.0;
    const double escape_radius = 20.0;
    const int max_iter = 128;
    const std::complex<double> c(-0.744, 0.148);

    // Speicher für das Bild auf der Host-Seite
    std::vector<int> image_host(width * height, 0);

    // Speicher für das Bild auf der GPU
    int* image_device;
    cudaMalloc(&image_device, width * height * sizeof(int));

    // Block- und Grid-Größe für CUDA, wobei nur x und y verwendet werden und z ignoriert wird
    dim3 block_size(16, 16);
    dim3 grid_size((width + block_size.x - 1) / block_size.x, (height + block_size.y - 1) / block_size.y);


    // CUDA-Kernel starten
    julia_kernel<<<grid_size, block_size>>>(width, height, x_min, x_max, y_min, y_max, escape_radius, max_iter, 
                                            c.real(), c.imag(), image_device);

    // Ergebnisse von der GPU auf die Host-Seite kopieren
    cudaMemcpy(image_host.data(), image_device, width * height * sizeof(int), cudaMemcpyDeviceToHost);

    // GPU-Speicher freigeben
    cudaFree(image_device);

    // Konvertieren der Iterationsanzahl in ein Graustufen-Bild
    std::vector<unsigned char> image_data(width * height); // Nur Graustufen (1 Byte pro Pixel)
    for (int i = 0; i < width * height; ++i) {
        image_data[i] = static_cast<unsigned char>(image_host[i]); // Graustufenwert (0-255)
    }

    // Speichern des Bildes als PNG
    const char* filename = "julia_set_grayscale.png";
    unsigned error = lodepng::encode(filename, image_data, width, height, LCT_GREY, 8);
    if (error) {
        std::cerr << "Fehler beim Speichern des Bildes: " << lodepng_error_text(error) << std::endl;
    } else {
        std::cout << "Bild erfolgreich gespeichert als " << filename << std::endl;
    }

    return 0;
}