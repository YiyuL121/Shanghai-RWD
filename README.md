# Shanghai-RWD Analysis

## Overview

This repository contains supporting code and model outputs for the results reported in *"Quantifying Heterogeneous Antiviral Efficacy of Nirmatrelvir–Ritonavir in SARS-CoV-2 Patients Using Real-World Data and Viral Dynamics Modeling."*

## Code

All analyses are implemented in R and organized into the following sections:

- **Part 1:** Propensity score matching
- **Part 2:** TREIV-PKPD model fitting (Monolix)
- **Part 3:** Association between treatment efficacy and patient characteristics
- **Part 4:** Sensitivity analysis

## Data Availability

Due to patient privacy constraints, individual-level clinical data (viral load trajectories, demographics, and treatment records) are not publicly available.

To support reproducibility, we provide the following Monolix model outputs, which are sufficient to replicate all downstream analyses from Part 2 onward:

- **`populationParameters.txt`** — Fixed-effect and variance parameter estimates from the TREIV-PKPD model
- **`estimatedIndividualParameters.txt`** — Empirical Bayes estimates (EBEs) of individual parameters for all 458 matched patients
- **`simulatedIndividualParameters.txt`** — Simulated individual parameters from the conditional distribution
- **`ind_eff.csv`** — Average individual treatment efficacy estimates for 229 treated individuals

These outputs enable full replication of the efficacy association analysis (Part 3) and sensitivity analysis (Part 4) without access to the raw data.

Researchers who wish to access the individual-level data for the purpose of replicating the model fitting step (Part 2) may contact the corresponding author to discuss a data use agreement.

## Model Specification

The viral dynamics model is a target-cell limited model with an eclipse phase, immune response, and viral production (TREIV-PKPD). Treatment effect is modeled as a reduction in viral production rate. The model was fit using Monolix with the SAEM algorithm. Full model equations and parameterization are described in the manuscript.

## Contact

For data access requests or questions about the analysis, please contact the corresponding author.
