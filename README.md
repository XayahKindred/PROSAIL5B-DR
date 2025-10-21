# PROSAIL5B-DR
An extended radiative transfer model incorporating dust retention effects on vegetation spectra (MATLAB implementation)

## Overview
The **PROSAIL5B-DR** model is an extended version of the classic PROSAIL5B radiative transfer model, designed to include the optical effects of **dust retention on vegetation leaves and canopies**. It provides a physically based framework to simulate spectral reflectance under dust-influenced conditions, which are common in urban and industrial environments.
This implementation maintains full compatibility with the original PROSAIL5B (PROSPECT + SAIL) structure while introducing additional parameters to represent dust layer absorption and scattering. The modifications are minimally intrusive, ensuring smooth integration into existing MATLAB workflows without altering the fundamental radiative transfer equations.

## File Structure
PROSAIL5B-DR/
├── main_PROSAIL_DR.m # Main driver script
├── PRO4SAIL.m # Core canopy radiative transfer module (SAIL)
├── prospect_DR.m # Leaf optical module (PROSPECT with dust)
├── dataSpec_DR.m # Spectral properties of dust materials
├── k_dust.mat # Spectral absorption and scattering coefficients of dust
├── Refl_CAN.txt # Example output: canopy reflectance spectra
│
├── campbell.m # Leaf angle distribution (Campbell function)
├── volscatt.m # Volume scattering function
├── tav.m, dcum.m, dladgen.m # Numerical utility functions
├── Jfunc1.m, Jfunc2.m, Jfunc3.m # Radiative integration functions
└── (other supporting MATLAB scripts)

## Citation
If you use this model in your research, please cite:
Yang, W. (2025). PROSAIL5B-DR: An extended radiative transfer model incorporating dust retention effects on vegetation spectra [Software]. Zenodo. https://doi.org/10.5281/zenodo.17402826.

## Contact
Author: Wei Yang
Institution: Shanghai Normal University
Email: 1000549339@smail.shnu.edu.cn OR yangw1663@gmail.com
For technical questions, collaboration inquiries, or suggestions regarding model improvement, please contact the author directly.

