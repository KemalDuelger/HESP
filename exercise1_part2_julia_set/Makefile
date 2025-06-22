# Makefile for mingw32-make on Windows

# Compiler
NVCC = nvcc.exe

# Source and output
SRC = julia-set.cu lodepng.cpp
OUT = julia-set.exe

# Default target
all: $(OUT)

$(OUT): $(SRC)
	$(NVCC) $(SRC) -o $(OUT)

# Clean target
clean:
	del /Q $(OUT)
