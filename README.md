# Optics and Electro-Optics Final Project
This project explores advanced concepts in optics, such as light propagation, polarization, lens modeling, and Fourier analysis, using MATLAB. It includes theoretical explanations, simulations, and visualizations to analyze various optical phenomena.

## Table of Contents
- [Overview](#overview)
- [Technologies Used](#technologies-used)
- [File Structure](#file-structure)
- [How to Run the Project](#how-to-run-the-project)
- [Results](#results)
- [License](#license)
- [Contributors](#contributors)

---

## Overview
### Question 1: Polarization and Brewster Angles
- **Objective**: Analyze light transmission through polarizers and glass plates to determine conditions for polarization and calculate Brewster angles.
- **Key Steps**:
  1. Simulate light interaction with a polarizer.
  2. Calculate transmission coefficients based on angle and refractive index.

### Question 2: Lens Modeling
- **Objective**: Model and simulate wavelength-dependent lens behavior using ray tracing and ABCD matrices.
- **Key Steps**:
  1. Simulate lens focal lengths based on refractive indices.
  2. Analyze ray progression through the optical system.

### Question 3: Fourier and Fresnel Transformations
- **Objective**: Perform diffraction and aperture analysis using Fourier transforms and Fresnel propagation.
- **Key Steps**:
  1. Simulate circular apertures and diffraction patterns.
  2. Analyze intensity distributions using custom MATLAB functions.

---

## Technologies Used
- **Language**: MATLAB
- **Tools**: Fourier Transform, Fresnel Propagation, ABCD Matrices
- **Visualization**: Surface Plots, Heatmaps, Custom Visualizations
- **Mathematical Models**: Polarization, Ray Tracing, Diffraction Analysis

---

## File Structure
```plaintext
.
├── Electro-optics part A.pdf             # Answers and code for Part A
├── Electro-optics part B.pdf             # Answers and code for Part B
├── Electro-optics_part1-Instructions.pdf # Instructions for Part 1
├── Electro-optics_part2-Instructions.pdf # Instructions for Part 2
├── Electro-optics_Q1.m                   # MATLAB script for Question 1
├── Electro-optics_Q2.m                   # MATLAB script for Question 2
├── Electro-optics_Q3.m                   # MATLAB script for Question 3
├── circ.m                                # Function to generate circular apertures
├── propFresnel.m                         # Function for Fresnel propagation
├── F.m                                   # Custom Fourier transform function
├── iF.m                                  # Custom inverse Fourier transform function
├── jinc.m                                # Jinc function for Fraunhofer diffraction
├── LICENSE                               # License file for the repository
└── README.md                             # Project documentation
```

## How to Run the Project

### Prerequisites
- MATLAB installed on your system.

### Instructions
1. Clone the repository:
   ```bash
   git clone https://github.com/12danielLL/electro-optics-project.git
   cd electro-optics-project
- run('Electro-optics_Q1.m')
- run('Electro-optics_Q2.m')
- run('Electro-optics_Q3.m')

## Visualize the Results
The scripts generate various plots and visualizations, such as diffraction patterns, ray tracing, and intensity distributions.

---

## Results

### Question 1: Polarization Analysis
- Simulated transmission coefficients for various polarizer angles.
- Calculated Brewster angle for optimal polarization.

### Question 2: Lens Behavior
- Modeled lens focal lengths for multiple wavelengths.
- Simulated ray tracing through lenses, showing clear progression.

### Question 3: Diffraction and Transformations
- Generated circular aperture diffraction patterns.
- Visualized intensity distributions using Fourier and Fresnel transforms.

---

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contributors
- **Daniel Brooker**  
  GitHub: [12danielLL](https://github.com/12danielLL)


## Acknowledgments

- Special thanks to Prof. Zeev Zalevsky for guidance and to course assistants Nadav Shvero and Sami Apsal for their support.
