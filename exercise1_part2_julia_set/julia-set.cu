#include <iostream>
#include <complex>
#include <vector>
#include <cuda_runtime.h>
#include "lodepng.h" // For saving as PNG

// CUDA method for calculating the Julia set
__global__ void julia_kernel(int width, int height, double x_min, double x_max, double y_min, double y_max, 
                             double max_value, int max_iter, double c_real, double c_imag, int* image) {
    int px = blockIdx.x * blockDim.x + threadIdx.x; // Pixel x-coordinate
    int py = blockIdx.y * blockDim.y + threadIdx.y; // Pixel y-coordinate

    if (px < width && py < height) {
        // Scale pixel coordinates to the complex plane y [-2, 2] x [-2, 2] as specified
        double x = x_min + (x_max - x_min) * px / (width - 1);
        double y = y_min + (y_max - y_min) * py / (height - 1);
        double z_real = x;
        double z_imag = y;

        // Iteration of the Julia function
        // GPU method does not support std::complex or the standard library, so manual calculation is used
        int iter = 0;
        while (z_real * z_real + z_imag * z_imag < max_value * max_value && iter < max_iter) {
            double temp = z_real * z_real - z_imag * z_imag + c_real;
            z_imag = 2.0 * z_real * z_imag + c_imag;
            z_real = temp;
            ++iter;
        }

        // Store the number of iterations in the image (grayscale mapping)
        image[py * width + px] = static_cast<int>(255.0 * iter / max_iter); // Scaled to 0-255
    }
}

int main() {
    // Constants for the Julia set
    const int width = 1024;
    const int height = 1024;
    const double x_min = -2.0, x_max = 2.0;
    const double y_min = -2.0, y_max = 2.0;
    const double max_value = 20.0;
    const int max_iter = 128;
    const std::complex<double> c(0.355, 0.355);

    // Memory for the image on the host side
    std::vector<int> image_host(width * height, 0);

    // Memory for the image on the GPU
    int* image_device;
    cudaMalloc(&image_device, width * height * sizeof(int));

    // Block and grid size for CUDA, using only x and y while ignoring z
    dim3 block_size(16, 16);
    dim3 grid_size((width + block_size.x - 1) / block_size.x, (height + block_size.y - 1) / block_size.y);

    // Launch the CUDA kernel
    julia_kernel<<<grid_size, block_size>>>(width, height, x_min, x_max, y_min, y_max, max_value, max_iter, 
                                            c.real(), c.imag(), image_device);

    // Copy results from the GPU to the host side
    cudaMemcpy(image_host.data(), image_device, width * height * sizeof(int), cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(image_device);

    // Convert the number of iterations into a grayscale image
    std::vector<unsigned char> image_data(width * height); // Grayscale only (1 byte per pixel)
    for (int i = 0; i < width * height; ++i) {
        image_data[i] = static_cast<unsigned char>(image_host[i]); // Grayscale value (0-255)
    }

    // Save the image as PNG
    const char* filename = "julia_set_grayscale.png";
    unsigned error = lodepng::encode(filename, image_data, width, height, LCT_GREY, 8);
    if (error) {
        std::cerr << "Error saving the image: " << lodepng_error_text(error) << std::endl;
    } else {
        std::cout << "Image successfully saved as " << filename << std::endl;
    }

    return 0;
}