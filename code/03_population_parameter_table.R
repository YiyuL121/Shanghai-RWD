# Population parameter estimates of the selected model (M4), with the fixed
# PK/PD and viral dynamics constants: manuscript Tables 2-3.
# Reads ind_metrics.csv from script 01 (or the copy in data/). Run from the repository root.

data_dir <- file.path("data", "M4_primary")
pop_file <- file.path(data_dir, "populationParameters.txt")
out_dir  <- file.path("output", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pop <- read.csv(pop_file, stringsAsFactors = FALSE)
pop$parameter <- trimws(pop$parameter)
get_est <- function(name) as.numeric(pop$value[pop$parameter == name])
get_se  <- function(name) as.numeric(pop$se_sa[pop$parameter == name])

# From script 01 if it has been run; otherwise the copy shipped with the data.
ind_metrics_file <- file.path(out_dir, "ind_metrics.csv")
if (!file.exists(ind_metrics_file)) ind_metrics_file <- file.path(data_dir, "ind_metrics.csv")
ind_metrics <- read.csv(ind_metrics_file, stringsAsFactors = FALSE)
eps_ave_pop <- mean(ind_metrics$eps_ave, na.rm = TRUE) * 100

row_est <- function(parameter, symbol, est_name, omega_name, distribution) {
  data.frame(Parameter = parameter, Symbol = symbol,
             `Population Mean` = signif(get_est(est_name), 3),
             `Standard Error`  = signif(get_se(est_name), 3),
             `SD of Random Effect` = signif(get_est(omega_name), 3),
             Distribution = distribution, Source = "Estimated",
             check.names = FALSE, stringsAsFactors = FALSE)
}
row_fixed <- function(parameter, symbol, value, source = "Fixed") {
  data.frame(Parameter = parameter, Symbol = symbol,
             `Population Mean` = value, `Standard Error` = "--",
             `SD of Random Effect` = "--", Distribution = "--", Source = source,
             check.names = FALSE, stringsAsFactors = FALSE)
}
row_beta <- function(parameter, symbol, est_name) {
  data.frame(Parameter = parameter, Symbol = symbol,
             `Population Mean` = signif(get_est(est_name), 3),
             `Standard Error`  = signif(get_se(est_name), 3),
             `SD of Random Effect` = "--", Distribution = "Fixed effect (covariate)",
             Source = "Estimated", check.names = FALSE, stringsAsFactors = FALSE)
}

viral_dynamics <- rbind(
  row_est("Viral infectivity (log10(RNA copies/mL)-1 day-1)", "log10beta", "log10beta_pop", "omega_log10beta", "Normal"),
  row_est("Viral production rate (log10 day-1)", "log10pi", "log10pi_pop", "omega_log10pi", "Normal"),
  row_est("Refractory to susceptible reversion rate (log10 day-1)", "log10rho", "log10rho_pop", "omega_log10rho", "Normal"),
  row_est("Susceptible to refractory conversion rate (log10 cell-1 day-1)", "log10phi", "log10phi_pop", "omega_log10phi", "Normal"),
  row_est("Infected cell clearance rate (day-1), reference group*", "delta", "delta_pop", "omega_delta", "Lognormal"),
  row_beta("Covariate effect on delta: age >=65 vs 18-64 (log scale)", "beta_age,delta", "beta_delta_age_cat___65"),
  row_beta("Covariate effect on delta: vaccination >=2 vs <2 doses (log scale)", "beta_vac,delta", "beta_delta_Vaccination_cat_More_than_2_doses"),
  row_est("Delay between infection and detection (days)", "tau", "tau_pop", "omega_tau", "Logit[0,20]"),
  row_fixed("Initial viral inoculum (RNA copies/mL)", "V(0)", "1"),
  row_fixed("Viral clearance rate (day-1)", "gamma", "15", "Fixed (Goyal 2020)"),
  row_fixed("Mean eclipse phase duration (days)", "1/k", "0.25", "Fixed (Ke 2021)"),
  row_fixed("Initial number of susceptible cells", "S(0)", "1e7"),
  row_fixed("Initial number of refractory cells", "R(0)", "0"),
  row_fixed("Initial number of productively infected cells", "Ip(0)", "0"),
  row_fixed("Initial number of infected cells in eclipse phase", "Ie(0)", "0"),
  data.frame(Parameter = "Residual error, additive on log10 VL scale", Symbol = "a",
             `Population Mean` = signif(get_est("a"), 3), `Standard Error` = signif(get_se("a"), 3),
             `SD of Random Effect` = "--", Distribution = "Additive", Source = "Estimated",
             check.names = FALSE, stringsAsFactors = FALSE)
)

pkpd <- rbind(
  row_fixed("Absorption rate constant (day-1)", "Ka", "9.98", "Fixed (Esmaeili 2024)"),
  row_fixed("Rate constant of transfer from plasma to lung (day-1)", "KPL", "1.58", "Fixed (Esmaeili 2024)"),
  row_fixed("Rate constant of transfer from lung to plasma (day-1)", "KLP", "1.22", "Fixed (Esmaeili 2024)"),
  row_fixed("Elimination rate constant (day-1)", "KCL", "4.96", "Fixed (Esmaeili 2024)"),
  row_fixed("Plasma volume of distribution (mL)", "Vol", "41,743", "Fixed (Esmaeili 2024)"),
  row_fixed("Maximum efficacy (%)", "Emax", "99.9", "Fixed (Esmaeili 2024)"),
  row_est("In-vivo effective concentration (uM), reference group*", "EC50", "IC50_pop", "omega_IC50", "Lognormal"),
  row_beta("Covariate effect on EC50: vaccination >=2 vs <2 doses (log scale)", "beta_vac,EC50", "beta_IC50_Vaccination_cat_More_than_2_doses"),
  row_fixed("Hill coefficient", "n", "3.16", "Fixed (Esmaeili 2024)"),
  row_fixed("Average treatment efficacy (%)", "eps_ave", signif(eps_ave_pop, 3), "Computed")
)

pop_table <- rbind(
  cbind(Model = "Viral dynamics model", viral_dynamics),
  cbind(Model = "PK/PD model", pkpd)
)

# * delta and EC50 are covariate models (age and/or vaccination shift the
# typical value away from this reference-group mean); covariate effects are
# listed directly below each parameter.

print(pop_table, row.names = FALSE)
write.csv(pop_table, file.path(out_dir, "Tables2_3_population_parameters.csv"), row.names = FALSE)
