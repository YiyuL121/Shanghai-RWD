# Sensitivity analysis (Supplementary Table 2): alternative structural models
# without covariates (M1-M3). Each treated patient's efficacy is computed from
# that model's individual estimates and regressed on host factors; AIC is
# compared with the selected model (M4). Run from the repository root.

library(dplyr)
library(deSolve)

data_root <- "data"
out_dir   <- file.path("output", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Fixed PK/PD constants (Esmaeili 2024; same as script 01) --
Tmin <- 0; Tmax <- 40; step_size <- 0.1
times <- round(seq(Tmin, Tmax, step_size), 1)
S0 <- 1e7; V0 <- 1; k_eclipse <- 4; gamma <- 15; LOD <- 1.31
E_max <- 0.999; K_a <- 9.98; K_PL <- 1.58; K_LP <- 1.22; K_CL <- 4.96
Vol <- 41743; MolMass <- 499.5; hill_n <- 3.16
dose_mg <- 300; dose_interval <- 0.5; treat_days <- 5
state_vars <- c("S", "R", "IE", "IP", "V", "A_G", "A_P", "A_L")

# 1. Simulator (same model as script 01)

make_derivs <- function(beta, pi, phi, rho, delta, on, EC50) {
  function(t, y, parms) {
    with(as.list(y), {
      C_uM <- A_P * (1e6 / (Vol * MolMass))
      eps  <- if (on) E_max * C_uM^hill_n / (EC50^hill_n + C_uM^hill_n) else 0

      dA_G <- -K_a * A_G
      dA_P <-  K_a * A_G + K_LP * A_L - (K_CL + K_PL) * A_P
      dA_L <-  K_PL * A_P - K_LP * A_L

      dS  <- -beta * S * V - phi * IP * S + rho * R
      dR  <-  phi * IP * S - rho * R
      dIE <-  beta * S * V - k_eclipse * IE
      dIP <-  k_eclipse * IE - delta * IP
      dV  <- (1 - eps) * pi * IP - gamma * V

      list(c(dS, dR, dIE, dIP, dV, dA_G, dA_P, dA_L), eps = eps, C_uM = C_uM)
    })
  }
}

simulate_viral_load <- function(beta, pi, phi, rho, delta, tau, EC50,
                                t_treat = NA, treated = 0) {
  treated   <- as.numeric(treated)
  rec       <- if (treated == 1) round(t_treat + tau, 1) else Inf
  treat_end <- rec + treat_days

  run_phase <- function(t_from, t_to, y_start, on, dose_events = NULL) {
    t_grid <- times[times >= t_from & times <= t_to]
    if (length(t_grid) < 2) {
      C_uM <- unname(y_start["A_P"]) * (1e6 / (Vol * MolMass))
      eps  <- if (on) E_max * C_uM^hill_n / (EC50^hill_n + C_uM^hill_n) else 0
      return(data.frame(time = t_from, t(y_start), eps = eps, C_uM = C_uM))
    }
    derivs <- make_derivs(beta, pi, phi, rho, delta, on, EC50)
    as.data.frame(ode(y = y_start, times = t_grid, func = derivs, parms = NULL,
                      events = dose_events, method = "lsoda"))
  }

  y0 <- c(S = S0, R = 0, IE = 0, IP = 0, V = V0, A_G = 0, A_P = 0, A_L = 0)

  if (treated == 0 || !is.finite(rec) || rec >= Tmax) {
    out <- run_phase(Tmin, Tmax, y0, on = FALSE)
  } else {
    pre   <- run_phase(Tmin, rec, y0, on = FALSE)
    y_rec <- unlist(pre[nrow(pre), state_vars])
    dose_times  <- seq(rec, rec + treat_days - dose_interval, by = dose_interval)
    dose_times  <- dose_times[dose_times <= min(treat_end, Tmax)]
    dose_events <- list(data = data.frame(var = "A_G", time = dose_times, value = dose_mg, method = "add"))
    on_treat <- run_phase(rec, treat_end, y_rec, on = TRUE, dose_events = dose_events)
    if (treat_end >= Tmax) {
      out <- rbind(pre, on_treat[-1, ])
    } else {
      y_end <- unlist(on_treat[nrow(on_treat), state_vars])
      post  <- run_phase(treat_end, Tmax, y_end, on = FALSE)
      out   <- rbind(pre, on_treat[-1, ], post[-1, ])
    }
  }
  data.frame(time = out$time, V_log10 = log10(pmax(out$V, 1e-30)), eps = out$eps, C_uM = out$C_uM)
}

compute_eps_ave <- function(sim, t_start) {
  w <- sim[sim$time >= t_start & sim$time <= t_start + treat_days, ]
  if (nrow(w) < 2) return(NA_real_)
  sum(diff(w$time) * (head(w$eps, -1) + tail(w$eps, -1)) / 2) / treat_days
}

# 2. Load one candidate model's Monolix output + host-factor covariates

load_candidate_model <- function(model_dir) {
  pop <- read.csv(file.path(model_dir, "populationParameters.txt"), stringsAsFactors = FALSE)
  pop$parameter <- trimws(pop$parameter)
  get_est <- function(name) as.numeric(pop$value[pop$parameter == name])

  Est <- read.csv(file.path(model_dir, "estimatedIndividualParameters.txt"), stringsAsFactors = FALSE)

  if (!"log10beta_mode" %in% names(Est)) Est$log10beta_mode <- get_est("log10beta_pop")
  if (!"log10pi_mode"   %in% names(Est)) Est$log10pi_mode   <- get_est("log10pi_pop")

  Est$age_cat  <- factor(Est$age_cat, levels = c("18-64", ">=65"))
  Est$sex      <- factor(Est$sex)
  Est$comorbidity_cat <- factor(as.character(Est$comorbidity_cat), levels = c("0", ">=1"))
  Est$vac_group <- factor(ifelse(Est$Vaccination_cat == "More than 2 doses", ">=2 doses", "<2 doses"),
                          levels = c("<2 doses", ">=2 doses"))
  if ("Severity" %in% names(Est)) {
    Est$severity_yn <- factor(ifelse(Est$Severity %in% c("Severe", "Critical", "Death"), "Yes", "No"),
                              levels = c("No", "Yes"))
  }
  Est
}

# 3. Compute each treated patient's efficacy from their own parameters

compute_individual_efficacy <- function(Est) {
  treated <- Est[Est$antiviral == 1, ]
  eff <- data.frame(id = treated$id, eps_ave = NA_real_)

  for (i in seq_len(nrow(treated))) {
    r <- treated[i, ]
    on <- simulate_viral_load(10^r$log10beta_mode, 10^r$log10pi_mode, 10^r$log10phi_mode,
                              10^r$log10rho_mode, r$delta_mode, r$tau_mode, r$IC50_mode,
                              t_treat = r$T_treat, treated = 1)
    eff$eps_ave[i] <- compute_eps_ave(on, round(r$T_treat + r$tau_mode, 1))
  }
  eff
}

# 4. Host-factor association: two-step EBE regression on computed efficacy

ebe_regression <- function(data, outcome, rhs) {
  f <- reformulate(rhs, response = outcome)
  m <- lm(f, data = data)
  s <- summary(m)
  ci <- confint(m)
  data.frame(outcome = outcome, term = rownames(s$coefficients),
             estimate = s$coefficients[, 1], se = s$coefficients[, 2],
             lo = ci[, 1], hi = ci[, 2], p_value = s$coefficients[, 4],
             row.names = NULL) %>%
    filter(term != "(Intercept)")
}

# 5. Run one candidate model end to end

get_aic <- function(model_dir) {
  ll <- read.csv(file.path(model_dir, "logLikelihood.txt"), stringsAsFactors = FALSE)
  ll$importanceSampling[ll$criteria == "AIC"]
}

run_candidate_model <- function(label, model_dir, direct_efficacy_col = NULL,
                                rhs = c("age_cat", "sex", "comorbidity_cat", "vac_group", "severity_yn")) {
  cat("\n===", label, "===\n")
  Est <- load_candidate_model(model_dir)
  treated <- Est[Est$antiviral == 1, ]

  if (!is.null(direct_efficacy_col)) {
    treated$eps_ave <- treated[[direct_efficacy_col]]   # M3 estimates efficacy directly
  } else {
    treated <- treated %>% left_join(compute_individual_efficacy(Est), by = "id")
  }

  res <- ebe_regression(treated, "eps_ave", rhs) %>% mutate(model = label)
  print(res %>% mutate(across(c(estimate, se, lo, hi), ~ round(.x, 3)), p_value = signif(p_value, 3)))
  list(reg = res, aic = get_aic(model_dir))
}

# 6. Candidate models M1-M3 and the selected model M4

results <- list(
  M1 = run_candidate_model("M1 baseline", file.path(data_root, "M1_baseline")),
  M2 = run_candidate_model("M2 fixed upslope (beta, pi fixed)", file.path(data_root, "M2_fixed_upslope")),
  M3 = run_candidate_model("M3 direct treatment effect (no PK/PD)", file.path(data_root, "M3_direct_effect"),
                           direct_efficacy_col = "eps_mode")
)
selected_aic <- get_aic(file.path(data_root, "M4_primary"))

# 7. Supplementary Table 2: AIC, delta-AIC (vs M4) and age/vaccination evidence

sig <- function(reg, term) ifelse(reg$p_value[reg$term == term] < .05, "Yes", "No")
suptab2 <- bind_rows(lapply(names(results), function(nm) {
  r <- results[[nm]]
  data.frame(Model = nm, AIC = round(r$aic, 0), `Delta-AIC` = round(r$aic - selected_aic, 0),
             Age = sig(r$reg, "age_cat>=65"), Vaccination = sig(r$reg, "vac_group>=2 doses"),
             check.names = FALSE)
}))
print(suptab2, row.names = FALSE)
write.csv(suptab2, file.path(out_dir, "SupplementaryTable2_sensitivity_models.csv"), row.names = FALSE)
