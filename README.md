# Shanghai-RWD

Code and model outputs for **"Host Factors Modulate Nirmatrelvir–Ritonavir Efficacy in COVID-19 Patients: A Viral Dynamics Modeling Study."**

A viral dynamics model with a two-compartment PK/PD model of nirmatrelvir was fitted in Monolix (SAEM) to 458 propensity-score-matched patients (229 treated, 229 untreated) from Shanghai. The selected model (M4) includes age and vaccination as covariates on the infected cell clearance rate (δ) and vaccination as a covariate on the in vivo EC₅₀. The R scripts here take the Monolix outputs and produce the simulations, tables and figures reported in the paper.

## Repository structure

```
code/
  00_matching.R                     cohort construction and propensity-score matching (needs clinical data)
  01_viral_dynamics_simulation.R    Figure 4, individual efficacy metrics, Supplementary Figures 3-4
  02_host_factor_associations.R     Figure 3, Supplementary Table 4, Supplementary Figures 2 and 5
  03_population_parameter_table.R   Tables 2-3
  04_sensitivity_models.R           Supplementary Table 2
data/
  M4_primary/         selected model
  M1_baseline/        sensitivity models (no covariates)
  M2_fixed_upslope/
  M3_direct_effect/
```

Each `data/` folder holds Monolix outputs: `populationParameters.txt`, `estimatedIndividualParameters.txt` and `logLikelihood.txt`. `M4_primary/` also has `simulatedIndividualParameters.txt` (50 draws per patient from the conditional distribution) and `ind_metrics.csv` (individual treatment metrics computed by script 01). Patient IDs are study IDs (1–458).

## Data availability

**Provided:** population parameters, individual parameter estimates (with the categorical covariates used in the analysis: age group, sex, comorbidity, vaccination status, severity, treatment arm and detection-to-treatment interval), simulated individual parameters and likelihood criteria for every model.

**Not provided:** individual-level clinical data (viral load measurements, continuous age, baseline viral load). Because of privacy and ethical restrictions, these are available to qualified researchers on reasonable request to the corresponding author, subject to ethical, institutional and data-sharing requirements.

## Requirements

- R ≥ 4.3 (tested with R 4.6.0 on Windows 11)
- R packages: dplyr, deSolve, ggplot2, emmeans, ggtext, patchwork (and MatchIt, optmatch, rpart for `00_matching.R`)
- Monolix 2024R1 was used for model fitting; it is not needed to run the code here.

```r
install.packages(c("dplyr", "deSolve", "ggplot2", "emmeans", "ggtext", "patchwork"))
```

Installation takes a few minutes on a standard computer. No special hardware is required.

## How to run

`00_matching.R` documents how the matched cohort was built. It needs the clinical data, so it cannot be run from this repository (see below); the Monolix fits were made on its output. Starting from the shared Monolix outputs, run the remaining scripts in order from the repository root:

```r
source("code/01_viral_dynamics_simulation.R")   # about 40 minutes
source("code/02_host_factor_associations.R")    # seconds
source("code/03_population_parameter_table.R")  # seconds
source("code/04_sensitivity_models.R")          # under 1 minute
```

Scripts 02 and 03 use the individual metrics written by script 01. To skip the 40-minute simulation, they fall back to the copy shipped in `data/M4_primary/ind_metrics.csv`. Results are written to `output/tables/` and `output/figures/`.

## Outputs and reproducibility

| Manuscript item | Script | Output | Reproducible from this repository |
|---|---|---|---|
| Matched cohort (Supplementary Figure 1) | 00 | `data_private/patient_data.csv` | No: needs the clinical data |
| Tables 2–3 | 03 | `Tables2_3_population_parameters.csv` | Yes |
| Figure 3 | 02 | `Fig3_delta_efficacy_by_hostfactor.pdf` | Yes |
| Figure 4 | 01 | `Fig4_population_trajectories.pdf` | Yes |
| Supplementary Table 2 (M1–M3) | 04 | `SupplementaryTable2_sensitivity_models.csv` | Yes |
| Supplementary Table 4 | 02 | `SupplementaryTable4_associations.csv` | Yes |
| Supplementary Figures 3–4 | 01 | `SFig3_individual_fits_untreated.pdf`, `SFig4_individual_fits_treated.pdf` | Model fits and 95% bands; observed viral loads need the clinical data |
| Supplementary Figure 5 | 02 | `SFig5_parameter_distributions.pdf` | Yes |
| Supplementary Figure 2 | 02 | `SFig2_age_vaccination_efficacy.pdf` | Needs the clinical data (continuous age) |
| Timing and baseline viral load vs efficacy (Results text) | 02 | `timing_baseline_vl_vs_efficacy.csv` | Needs the clinical data (baseline viral load) |
| Population-model covariate tests | 02 | `population_model_covariate_effects.csv` | Yes |

Table 1, Figures 1–2, the Supplementary Figure 1 flow chart, Supplementary Table 3 and the Monolix fits themselves were prepared outside these scripts and are not produced here.

**With the clinical data** (`combine_df.RData`, one row per PCR test), place it at `data_private/combine_df.RData` (or set `RWD_RAW_DATA`) and run `00_matching.R` first. It writes the matched cohort to `data_private/patient_data.csv` with the same study IDs as `data/`. Scripts 01 and 02 then also produce the items marked as needing clinical data. Another location can be set with `RWD_PATIENT_DATA`; the file needs one row per viral load measurement, with at least `id`, `Day_since_detection`, `VL` (log10 copies/mL), `age` (years) and `init_VL` (first measured viral load, log10 copies/mL).

## Methods notes

- **Cohort and matching.** Viral load (log10 copies/mL) is derived from N-gene Ct as −0.32·Ct + 14.11, with Ct 40 (1.31 log10) as the detection limit. Patients with severe kidney or liver disease, incomplete screening or outcome dates, treatment more than 7 days after detection, or oxygen requirement by the reference date (treatment start, or day 4 after detection for controls) are excluded. Treated and untreated patients are matched 1:1 by optimal matching on a classification-tree (rpart) propensity score, using age group, sex, vaccination (<2 doses, 2 doses, booster), comorbidity and first measured viral load.

- **Simulation.** Viral dynamics (target cells, refractory cells, eclipse and productively infected cells, virus) are coupled to a two-compartment PK model of nirmatrelvir (300 mg twice daily for 5 days). Drug efficacy follows a Hill function of plasma concentration. Fixed constants are listed in Tables 2–3. Each patient's average treatment efficacy (εave) is the mean drug efficacy over the 5-day course, simulated from that patient's mode estimates.
- **Individual fits.** The best fit uses each patient's mode estimates. The 95% band comes from Monolix's simulated individual parameters, with the infection-to-detection delay τ held at its mode.
- **Host-factor associations.** For the three covariate effects in the population model (age and vaccination on δ, vaccination on EC₅₀), Wald tests on the Monolix estimates are reported. All other cells in Supplementary Table 4 come from multivariable linear regression of individual mode estimates on age group, sex, comorbidity, vaccination and subsequent oxygen requirement. All participants are used for log10β, log10π, log10ρ, log10φ, δ and τ; treated participants only for EC₅₀ and εave. δ and EC₅₀ are regressed on the log scale and reported as ratios. Figure 3 shows the same regression's adjusted means (emmeans) for the displayed factor.
- **Sensitivity models.** For M1–M3, which have no covariates, εave is computed for each treated patient (estimated directly in M3) and regressed on the same host factors. ΔAIC is relative to M4.

## License

MIT License. See `LICENSE`.

## Contact

For questions about the code, please open an issue. For data access, please contact the corresponding author.
