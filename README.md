# Shanghai-RWD Analysis

## Overview

This repository contains supporting code and model outputs for the results reported in *"Quantifying Heterogeneous Antiviral Efficacy of Nirmatrelvir–Ritonavir in SARS-CoV-2 Patients Using Real-World Data and Viral Dynamics Modeling."*

## System Requirements

### Software Dependencies

- **R** (>= 4.2.0)
- **Monolix** (2023R1 or later) — required only for model fitting (Part 2); not needed for downstream analyses
- R packages: `dplyr`, `ggplot2`, `GGally`, `tidyr`, `finalfit`, `emmeans`, `deSolve`, `MatchIt`

### Operating Systems

The code has been tested on:

- Windows 10/11
- macOS 13 (Ventura)

### Hardware

No non-standard hardware is required. A standard desktop or laptop computer is sufficient.

## Installation Guide

### Instructions

1. Install [R](https://cran.r-project.org/) (>= 4.2.0) and optionally [RStudio](https://posit.co/download/rstudio-desktop/).
2. Install required R packages by running:

```r
install.packages(c("dplyr", "ggplot2", "GGally", "tidyr", "finalfit", "emmeans", "deSolve", "MatchIt"))
```

3. Clone or download this repository:

```bash
git clone https://github.com/YiyuL121/Shanghai-RWD.git
```

### Typical Install Time

Under 5 minutes on a standard desktop computer (assuming R is already installed).

## Code

All analyses are implemented in R and organized into the following sections:

- **Part 1:** Propensity score matching
- **Part 2:** TREIV-PKPD model fitting (Monolix)
- **Part 3:** Association between treatment efficacy and patient characteristics
- **Part 4:** Sensitivity analysis (constant-efficacy model)

## Data

### Provided in This Repository

- **`populationParameters.txt`** — Fixed-effect and variance parameter estimates from the TREIV-PKPD model
- **`estimatedIndividualParameters.txt`** — Empirical Bayes estimates (EBEs) of individual parameters for all 458 matched patients
- **`simulatedIndividualParameters.txt`** — Simulated individual parameters from the conditional distribution
- **`ind_eff.csv`** — Average individual treatment efficacy estimates for 229 treated individuals

### Not Provided (Privacy Restrictions)

Due to patient privacy constraints, individual-level clinical data (viral load trajectories, demographics, and treatment records) are not publicly available. Researchers who wish to access the underlying data may contact the corresponding author to discuss a data use agreement.

### Reproducibility by Section

| Section | Reproducible? | Required files |
|---------|:---:|----------------|
| Part 1 — PSM | No | Raw patient data (not provided) |
| Part 2 — Monolix fitting | No | Raw viral load and dosing data (not provided) |
| Part 3 — Efficacy associations | **Yes** | `estimatedIndividualParameters.txt`, `ind_eff.csv` (provided) |
| Part 4 — Sensitivity analysis | No | Individual parameters from a separate constant-efficacy Monolix run (not provided) |

## Demo

### Instructions

To reproduce the efficacy association analysis (Part 3):

1. Open `(3) Efficacy associations.Rmd` in RStudio.
2. Ensure `estimatedIndividualParameters.txt` and `ind_eff.csv` are in the working directory.
3. Run all chunks (or knit the document).

### Expected Output

- Multivariate regression results of treatment efficacy on patient characteristics (age, sex, comorbidity, vaccination status)
- Figure 4: Boxplots of efficacy by age group and vaccination dose with adjusted marginal means
- Supplementary Figures 4–6: Individual parameter distributions, age–dose–efficacy scatter plots, and subgroup-averaged viral load trajectories

### Expected Run Time

Under 1 minute on a standard desktop computer.

## Instructions for Use

### Running on Your Own Data

To apply the analysis pipeline to a new dataset:

1. **Part 1 (PSM):** Prepare a patient-level dataset with treatment assignment, demographics (age, sex), vaccination status, and comorbidity. Run `(1) Propensity Score Matching` to generate a matched cohort.
2. **Part 2 (Monolix fitting):** Format the matched cohort as a Monolix-compatible dataset with viral load observations, dosing records, and regressors (`antiviral`, `t_treat`). Fit the TREIV-PKPD model in Monolix to obtain individual parameter estimates.
3. **Part 3 (Efficacy associations):** Use the Monolix output (`estimatedIndividualParameters.txt`) and computed individual efficacy (`ind_eff.csv`) as inputs to `(3) Efficacy associations.Rmd`.
4. **Part 4 (Sensitivity analysis):** Refit the model with a constant-efficacy specification in Monolix, then run `(4) Sensitivity Analysis` on the resulting individual parameters.

## Model Specification

The viral dynamics model is a target-cell limited model with an eclipse phase, immune response, and treatment-modulated viral production (TREIV-PKPD). Treatment effect is modeled as a reduction in viral production rate driven by a pharmacokinetic–pharmacodynamic (PKPD) sub-model of nirmatrelvir–ritonavir. The model was fit using Monolix with the SAEM algorithm. Full model equations and parameterization are described in the manuscript.

## License

This project is licensed under the MIT License.

## Contact

For data access requests or questions about the analysis, please contact the corresponding author.
