# README

## Overview  
This repository contains the implementation, simulation studies, and empirical applications accompanying the paper **Confounded Ordinary Differential Equation**.  
The project is organized into two main components:  

1. **Simulation**: Codes and results for numerical simulation experiments.  
2. **RealData**: Codes and results for empirical applications with real-world datasets.  

---

## Table of Contents  
- [README](#readme)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [1. Simulation](#1-simulation)
    - [(a) Simulation_Linear](#a-simulation_linear)
    - [(b) Simulation_Nonlinear](#b-simulation_nonlinear)
  - [2. RealData](#2-realdata)
    - [(a) ncov2019](#a-ncov2019)
    - [(b) GeneData](#b-genedata)

---

## 1. Simulation  
This directory contains the codes and outputs for numerical simulation studies, divided into two subdirectories:

### (a) Simulation_Linear  
Contains implementation and results for linear ODE numerical simulations.

- **`Linear_simulation.R`**: Implementation of linear ODE simulation; run this file to reproduce linear ODE simulation results.  
- **`Output/`**: Stores all simulation outputs from linear ODE experiments.

**Usage**: Run `Linear_simulation.R` to replicate the linear ODE simulation study.

### (b) Simulation_Nonlinear  
Contains implementation and results for nonlinear ODE numerical simulations.

- **`Nonlinear_simulation.R`**: Implementation of nonlinear ODE simulation; run this file to reproduce nonlinear ODE simulation results.  
- **`Output/`**: Stores all simulation outputs from nonlinear ODE experiments.

**Usage**: Run `Nonlinear_simulation.R` to replicate the nonlinear ODE simulation study.

---

## 2. RealData  
This directory contains applications to real-world datasets. It is divided into two submodules:  

### (a) ncov2019  
Contains COVID-19 pandemic data analysis.

- **`Preprocess_ncov.R`**: Preprocessing code for COVID-19 infection data. The original data is sourced from the R package `ncov2019`.
- **`ncov_est1.R`**: Network estimation using different methods on the complete dataset.  
- **`Preprocess_ncov2.R`**: Preprocessing code after removing data from France and Germany. The original data is sourced from the R package `ncov2019`.
- **`ncov_est2.R`**: Network estimation using different methods on the modified dataset (without France and Germany).  
- **`Output/`**: Stores estimation results for pandemic data.
  - `ncov_result1.RData`: Network estimation results from complete data.
  - `ncov_result2.RData`: Network estimation results from modified data (without France and Germany).

**Usage**: Run the scripts in the following order to replicate the COVID-19 data analysis:
1. `Preprocess_ncov.R`
2. `ncov_est1.R` 
3. `Preprocess_ncov2.R`
4. `ncov_est2.R`

### (b) GeneData  
Contains yeast gene regulatory network data analysis.

- **`Preprocess1.R`**: Preprocessing code for yeast gene expression data. The original data is sourced from the R package `GEOquery`.
- **`Rugmatrix_est1.R`**: Gene regulatory network estimation using different methods on the complete dataset.  
- **`Preprocess2.R`**: Preprocessing code after removing GAL3 gene expression data. The original data is sourced from the R package `GEOquery`.
- **`Rugmatrix_est2.R`**: Gene regulatory network estimation using different methods on the modified dataset (without GAL3).  
- **`Output/`**: Stores estimation results for gene expression data.
  - `Generesult1.RData`: Gene regulatory network estimation results from complete data.
  - `Generesult2.RData`: Gene regulatory network estimation results from modified data (without GAL3).

**Usage**: Run the scripts in the following order to replicate the gene regulatory network analysis:
1. `Preprocess1.R`
2. `Rugmatrix_est1.R`
3. `Preprocess2.R` 
4. `Rugmatrix_est2.R`

