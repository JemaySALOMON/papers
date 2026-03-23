---
title: "README for the article: 'Joint modeling of direct and social genetic effects in mono-genotypic and pluri-specific groups: case study in intercrops'"
author: "J. Salomon"
date: "19/03/2026"
---

# Overview
This repository contains supplementary information and code for the article:

**"Joint modeling of direct and social genetic effects in mono-genotypic and pluri-specific groups: case study in intercrops"**.

It includes simulation, inference, and plotting scripts to reproduce the main figures in the article.

# Supplementary Information
- **suppl-S1.pdf**: Default parameters used for the simulations.
- **suppl-S2.pdf**: Supplementary Table 1 and Supplementary Figures 1–12.

# R folder (`R/`)
The `R/` folder contains `.Rmd` files to run simulations, perform inference, gather results, and generate plots:

- `simulPhenos.Rmd` – Set up and run simulation experiments.
- `Inference.Rmd` – Perform statistical inference on simulated or real data.
- `GatherResults.Rmd` – Collect and summarize inference results.
- `Plots.Rmd` – Generate figures and visualizations from the results.
- `InferwithK_vs_InferwithID.Rmd` – Generate plots for Figure 5 and Suppl. Figure 8 

# src folder (`src/`)
The `src/` folder contains functions and compiled code for inference:

- `utils.R` and `utils.hpp` – Utility functions.
- `MMTMB.hpp` and `MMTMB.cpp` – Core C++ implementation for the model.
- `MMTMB.o` and `MMTMB.so` – Compiled object and shared library files.


# data folder(`data/`)
The `data/` folder contains simulated data for the main figures of the articles.


# How to work with the scripts

## 1. Install Dependencies
The scripts were built using **R version 4.4.1**. Install required packages as follows:

```r
Rscript install.R
```

## Reproducing the Figures of the Study

The main figures of the article can be reproduced by running the following scripts:

- `Plots.Rmd`
- `InferwithK_vs_inferwitID.Rmd`

Pre-rendered `.html` files are also provided for quick viewing without rerunning the analyses.
---

## Reproducing the Scenarios (Data Simulation)

### 1. Run `simulPhenos.Rmd`

This script creates a `paper/` directory containing a subdirectory `simul/` with all simulation objects.

### 2. Run `Inference.Rmd`

This script generates an `infer/` directory inside `paper/`, containing all inference objects.

### 3. Run `GatherResults.Rmd`

This script creates an `eval/` directory inside `paper/`, where the summarized inference results are stored.

### 4. Run `Plots.Rmd`

This script generates the plots for the considered scenarios.

**Note:** Parameters and file paths may need to be adjusted to match the specific scenarios, as the same script is used to produce the main figures of the article.

## Using Docker

For experienced users familiar with Docker and CLI, I've prebuilt an  **r-lab** image for use in Docker containers. Everything is already installed. See the **Docker.pdf** file for instructions.
