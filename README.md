# 2PSP: C++ Implementation

This repository contains the official C++ source code for the paper:

> [**Two-Pass Shock Propagation for Stable Stacking with Gauss-Seidel**] > Ziyan Xiong, Andrew Leach, Griffith Thomas, and Shinjiro Sueda
> *SCA*, 2025
> [Paper](https://dl.acm.org/doi/10.1145/3747867)

## Introduction

This code implements the [**Two-Pass Shock Propagation Method**] presented in our paper. The project is written in C++ and relies on several external libraries for tasks such as linear algebra, and mesh collsioin detecion.

The main components of the repository are:
* `/src`: Contains all the C++ source code (`.cpp`) and headers (`.h`).
* `/ShapeFiles`: Includes sample input meshes.
* `/Results`: The default directory where simulation results are saved.

---

## Build Instructions

This project uses **CMake** to manage the build process. The following dependencies are required.

### Dependencies
* **CMake** (version 3.10 or higher)
* A modern C++ compiler (GCC, Clang, or MSVC) that supports C++17.
* **[Eigen]**: For linear algebra operations.
* **[coal]**: For mesh collsion detection.

### Compiling
To build the project, clone the repository and use the following commands:

```bash
# Install coal
conda install coal -c conda-forge

# Set the path to Coal 
export COAL_ROOT='path to installed folder of coal'

# Set the path to Eigen
export EIGEN3_INCLUDE_DIR='path to eigen folder'

mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j4
```
An executable named `2PSP` will be created in the `build` directory.

---

## How to Run

After successfully building the project, you can run the executable from the `build` directory.

### Expected Results

Running the simulation will show the running time and average iteration count in the command line and produce output files in the `Results` directory.

1.  A txt file (`Results\Scene_id\Body_States_2PSP.txt`) for each scene contating the states for each object in the scene in each frame.
2.  A txt file (`Results\time_report.txt`) containing timing information for each scene in ms.