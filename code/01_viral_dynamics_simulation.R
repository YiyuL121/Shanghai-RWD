# Viral dynamics and PK/PD simulation from the selected model (M4).
# Writes Figure 4, individual treatment metrics (ind_metrics.csv, used by
# scripts 02 and 03) and the individual fits (Supplementary Figures 3-4).
# Run from the repository root.

library(dplyr)
library(deSolve)
library(ggplot2)
library(patchwork)

# ---- Paths ----------------------------------------------------------------
data_dir     <- file.path("data", "M4_primary")
pop_file     <- file.path(data_dir, "populationParameters.txt")
ebe_file     <- file.path(data_dir, "estimatedIndividualParameters.txt")
sim_file     <- file.path(data_dir, "simulatedIndividualParameters.txt")
# Observed viral loads are not public. If available (see README), they are
# added as points to the individual fits.
patient_file <- Sys.getenv("RWD_PATIENT_DATA", unset = file.path("data_private", "patient_data.csv"))
out_dir      <- file.path("output", "tables")
fig_dir      <- file.path("output", "figures")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(20260811)

# ---- Time grid -----------------------------------------------------------
Tmin      <- 0
Tmax      <- 40
step_size <- 0.1
times     <- round(seq(Tmin, Tmax, step_size), 1)

# ---- Fixed viral-dynamics constants ------------------------------------
S0        <- 1e7    # initial susceptible cells
V0        <- 1      # initial viral load, RNA copies/mL (below assay LOD)
k_eclipse <- 4      # 1 / eclipse-phase length (0.25 d)
gamma     <- 15     # free virion clearance rate, /d
LOD       <- 1.31   # assay detection limit, log10 copies/mL

# ---- Fixed PK/PD constants (Esmaeili 2024) -----------------------------
E_max         <- 0.999
K_a           <- 9.98
K_PL          <- 1.58
K_LP          <- 1.22
K_CL          <- 4.96
Vol           <- 41743
MolMass       <- 499.5
hill_n        <- 3.16
dose_mg       <- 300
dose_interval <- 0.5   # days between doses (q12h)
treat_days    <- 5     # nirmatrelvir-ritonavir course length

# ---- Monolix population parameters -------------------------------------
pop <- read.csv(pop_file, stringsAsFactors = FALSE)
pop$parameter <- trimws(pop$parameter)

get_par <- function(name) {
  v <- pop$value[pop$parameter == name]
  if (length(v) != 1) stop("Parameter '", name, "' not found in ", pop_file)
  as.numeric(v)
}
get_omega <- function(name) {
  v <- pop$value[pop$parameter == name]
  if (length(v) == 0) return(0)   # no random effect estimated for this parameter
  as.numeric(v)
}

# Population fixed effects, natural scale. The Monolix run exports "IC50";
# the manuscript calls the same quantity EC50.
beta_pop  <- 10^get_par("log10beta_pop")
pi_pop    <- 10^get_par("log10pi_pop")
phi_pop   <- 10^get_par("log10phi_pop")
rho_pop   <- 10^get_par("log10rho_pop")
delta_pop <- get_par("delta_pop")
tau_pop   <- get_par("tau_pop")
EC50_pop  <- get_par("IC50_pop")

# Covariate model: delta ~ age + vaccination, EC50 ~ vaccination.
beta_delta_age65 <- get_par("beta_delta_age_cat___65")
beta_delta_vac2  <- get_par("beta_delta_Vaccination_cat_More_than_2_doses")
beta_EC50_vac2   <- get_par("beta_IC50_Vaccination_cat_More_than_2_doses")

# Typical delta / EC50 for a covariate group (age_65, vac_2plus are 0/1 flags).
# The reference group (both 0) reproduces delta_pop / EC50_pop exactly.
# unname() strips names that would otherwise propagate into the ODE output.
group_delta <- function(age_65 = 0, vac_2plus = 0) {
  unname(delta_pop * exp(beta_delta_age65 * age_65 + beta_delta_vac2 * vac_2plus))
}
group_EC50 <- function(vac_2plus = 0) {
  unname(EC50_pop * exp(beta_EC50_vac2 * vac_2plus))
}

# ---- Sample population parameters (for the 95%/50% CI bands) ----------
sample_population_pars <- function(n, age_65 = 0, vac_2plus = 0) {
  logit     <- function(p) log(p / (20 - p))
  inv_logit <- function(x) 20 * exp(x) / (1 + exp(x))   # tau is bounded 0-20 d

  data.frame(
    beta  = 10^rnorm(n, log10(beta_pop), get_omega("omega_log10beta")),
    pi    = 10^rnorm(n, log10(pi_pop),   get_omega("omega_log10pi")),
    phi   = 10^rnorm(n, log10(phi_pop),  get_omega("omega_log10phi")),
    rho   = 10^rnorm(n, log10(rho_pop),  get_omega("omega_log10rho")),
    delta = exp(rnorm(n, log(group_delta(age_65, vac_2plus)), get_omega("omega_delta"))),
    tau   = inv_logit(rnorm(n, logit(tau_pop), get_omega("omega_tau"))),
    EC50  = exp(rnorm(n, log(group_EC50(vac_2plus)), get_omega("omega_IC50")))
  )
}

# ---- TREIV + PK/PD simulator -----------------------------------------
state_vars <- c("S", "R", "IE", "IP", "V", "A_G", "A_P", "A_L")

make_derivs <- function(beta, pi, phi, rho, delta, on, EC50) {
  function(t, y, parms) {
    with(as.list(y), {
      C_uM <- A_P * (1e6 / (Vol * MolMass))                  # plasma drug concentration
      eps  <- if (on) E_max * C_uM^hill_n / (EC50^hill_n + C_uM^hill_n) else 0

      dA_G <- -K_a * A_G
      dA_P <-  K_a * A_G + K_LP * A_L - (K_CL + K_PL) * A_P
      dA_L <-  K_PL * A_P - K_LP * A_L

      # S susceptible ("T" in the paper), R refractory, IE eclipse phase,
      # IP productively infected, V free virus
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
  rec       <- if (treated == 1) round(t_treat + tau, 1) else Inf   # treatment start since infection
  treat_end <- rec + treat_days

  run_phase <- function(t_from, t_to, y_start, on, dose_events = NULL) {
    t_grid <- times[times >= t_from & times <= t_to]

    if (length(t_grid) < 2) {
      # Zero-width phase (a sampled tau rounds rec to exactly t_from): ode()
      # needs >= 2 time points, so report the starting state instead.
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

    dose_times  <- seq(rec, rec + treat_days - dose_interval, by = dose_interval)  # 10 doses
    dose_times  <- dose_times[dose_times <= min(treat_end, Tmax)]
    dose_events <- list(data = data.frame(var = "A_G", time = dose_times,
                                          value = dose_mg, method = "add"))

    on_treat <- run_phase(rec, treat_end, y_rec, on = TRUE, dose_events = dose_events)

    if (treat_end >= Tmax) {
      out <- rbind(pre, on_treat[-1, ])
    } else {
      y_end <- unlist(on_treat[nrow(on_treat), state_vars])
      post  <- run_phase(treat_end, Tmax, y_end, on = FALSE)
      out   <- rbind(pre, on_treat[-1, ], post[-1, ])
    }
  }

  data.frame(time = out$time, V_log10 = log10(pmax(out$V, 1e-30)),
             eps = out$eps, C_uM = out$C_uM, A_P = out$A_P)
}

# ---- Summary metrics ------------------------------------------------
# Average drug efficacy over the 5-day treatment window.
compute_eps_ave <- function(sim, t_start) {
  w <- sim[sim$time >= t_start & sim$time <= t_start + treat_days, ]
  if (nrow(w) < 2) return(NA_real_)
  sum(diff(w$time) * (head(w$eps, -1) + tail(w$eps, -1)) / 2) / treat_days
}

# Shedding duration and AUC above the assay detection limit.
shedding_duration <- function(sim) {
  above <- sim$V_log10 > LOD
  if (!any(above)) return(0)
  i0 <- which(above)[1]
  below_after <- which(!above & seq_along(above) > i0)
  if (!length(below_after)) return(Tmax)   # still shedding at end of window
  sim$time[below_after[1]] - sim$time[i0]
}
auc_above_lod <- function(sim) {
  y <- pmax(sim$V_log10 - LOD, 0)
  sum(diff(sim$time) * (head(y, -1) + tail(y, -1)) / 2)
}

# ---- Shared plot theme -------------------------------------------
theme_vd <- function(base = 14) {
  theme_minimal(base_size = base) +
    theme(axis.text        = element_text(colour = "black"),
          axis.ticks       = element_line(colour = "black"),
          axis.line        = element_line(colour = "black"),
          panel.grid       = element_blank(),
          panel.background = element_blank(),
          strip.text       = element_text(size = base),
          strip.background = element_rect(fill = "grey85", colour = NA),
          legend.position  = "bottom",
          plot.margin      = margin(.5, .5, .5, .5, "cm"))
}

groups <- list(   # age_65, vac_2plus flags per covariate group
  "18-64, <2 doses"  = c(age_65 = 0, vac_2plus = 0),
  "18-64, >=2 doses" = c(age_65 = 0, vac_2plus = 1),
  ">=65, <2 doses"   = c(age_65 = 1, vac_2plus = 0),
  ">=65, >=2 doses"  = c(age_65 = 1, vac_2plus = 1)
)
t_treat_ref <- 4   # median detection -> treatment delay, days (Supp Table 1)

eps_by_group <- data.frame(group = names(groups), delta = NA_real_,
                           EC50 = NA_real_, eps_ave = NA_real_)
for (i in seq_along(groups)) {
  g       <- groups[[i]]
  delta_g <- group_delta(g["age_65"], g["vac_2plus"])
  EC50_g  <- group_EC50(g["vac_2plus"])
  sim <- simulate_viral_load(beta_pop, pi_pop, phi_pop, rho_pop, delta_g, tau_pop, EC50_g,
                             t_treat = t_treat_ref, treated = 1)

  eps_by_group$delta[i]   <- delta_g
  eps_by_group$EC50[i]    <- EC50_g
  eps_by_group$eps_ave[i] <- compute_eps_ave(sim, round(t_treat_ref + tau_pop, 1))
}
print(eps_by_group)

n_draws   <- 1000   # population draws for the 95%/50% CI envelope
scenarios <- data.frame(t_treat = c(0, 5, 0, 5), treated = c(1, 1, 0, 0))

Fit4 <- list()
for (g in seq_len(nrow(scenarios))) {
  t_treat <- scenarios$t_treat[g]
  treated <- scenarios$treated[g]

  best_fit <- simulate_viral_load(beta_pop, pi_pop, phi_pop, rho_pop, delta_pop, tau_pop, EC50_pop,
                                  t_treat = t_treat, treated = treated)

  draws   <- sample_population_pars(n_draws)
  V_draws <- sapply(seq_len(n_draws), function(i) {
    simulate_viral_load(draws$beta[i], draws$pi[i], draws$phi[i], draws$rho[i],
                        draws$delta[i], draws$tau[i], draws$EC50[i],
                        t_treat = t_treat, treated = treated)$V_log10
  })

  Fit4[[g]] <- data.frame(
    time     = best_fit$time,
    best_fit = pmax(pmin(best_fit$V_log10, 11), -2),
    A_P      = best_fit$A_P,
    Min95    = pmax(apply(V_draws, 1, quantile, 0.025, na.rm = TRUE), -2),
    Max95    = pmin(apply(V_draws, 1, quantile, 0.975, na.rm = TRUE), 11),
    Min50    = pmax(apply(V_draws, 1, quantile, 0.25,  na.rm = TRUE), -2),
    Max50    = pmin(apply(V_draws, 1, quantile, 0.75,  na.rm = TRUE), 11),
    treated  = treated,
    t_treat_since_inf = round(t_treat + tau_pop, 1),
    panel    = paste("Treatment day", t_treat)
  )
}
fig4 <- do.call(rbind, Fit4)
fig4$treated <- factor(fig4$treated, levels = c(0, 1))

p4a <- ggplot(fig4, aes(time, best_fit, group = treated)) +
  geom_ribbon(data = filter(fig4, treated == "1"),
              aes(ymin = Min95, ymax = Max95, fill = treated), alpha = .1) +
  geom_ribbon(data = filter(fig4, treated == "1"),
              aes(ymin = Min50, ymax = Max50, fill = treated), alpha = .2) +
  geom_line(aes(colour = treated), linewidth = .75, alpha = .8) +
  geom_vline(aes(xintercept = t_treat_since_inf), colour = "green",
             linetype = "dotted", linewidth = 1) +
  geom_vline(aes(xintercept = t_treat_since_inf + treat_days), colour = "green",
             linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = LOD, colour = "grey", linetype = "dotted", linewidth = 1) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#d62728"),
                      labels = c("Untreated", "Treated")) +
  scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#d62728"), guide = "none") +
  labs(x = "Days after infection",
       y = expression(log[10]~"viral load (RNA copies/mL)"), colour = NULL) +
  coord_cartesian(xlim = c(0, 30), ylim = c(-2, 11), expand = FALSE) +
  theme_vd()

p4b <- ggplot(fig4, aes(time, A_P, group = treated, colour = treated)) +
  geom_line(linewidth = .75, alpha = .8) +
  geom_vline(aes(xintercept = t_treat_since_inf), colour = "green",
             linetype = "dotted", linewidth = 1) +
  geom_vline(aes(xintercept = t_treat_since_inf + treat_days), colour = "green",
             linetype = "dotted", linewidth = 1) +
  facet_wrap(~ panel) +
  scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#d62728"),
                      labels = c("Untreated", "Treated")) +
  labs(x = "Days after infection", y = "Drug amount in plasma (mg)", colour = NULL) +
  coord_cartesian(xlim = c(0, 20)) +
  theme_vd()

Fit4c <- list()
k     <- 1
for (label in names(groups)) {
  g       <- groups[[label]]
  delta_g <- group_delta(g["age_65"], g["vac_2plus"])
  EC50_g  <- group_EC50(g["vac_2plus"])

  for (treated in c(1, 0)) {
    best_fit <- simulate_viral_load(beta_pop, pi_pop, phi_pop, rho_pop, delta_g, tau_pop, EC50_g,
                                    t_treat = t_treat_ref, treated = treated)

    draws   <- sample_population_pars(n_draws, age_65 = g["age_65"], vac_2plus = g["vac_2plus"])
    V_draws <- sapply(seq_len(n_draws), function(i) {
      simulate_viral_load(draws$beta[i], draws$pi[i], draws$phi[i], draws$rho[i],
                          draws$delta[i], draws$tau[i], draws$EC50[i],
                          t_treat = t_treat_ref, treated = treated)$V_log10
    })

    Fit4c[[k]] <- data.frame(
      time     = best_fit$time,
      best_fit = pmax(pmin(best_fit$V_log10, 11), -2),
      Min95    = pmax(apply(V_draws, 1, quantile, 0.025, na.rm = TRUE), -2),
      Max95    = pmin(apply(V_draws, 1, quantile, 0.975, na.rm = TRUE), 11),
      Min50    = pmax(apply(V_draws, 1, quantile, 0.25,  na.rm = TRUE), -2),
      Max50    = pmin(apply(V_draws, 1, quantile, 0.75,  na.rm = TRUE), 11),
      treated  = treated,
      age_lab  = ifelse(g["age_65"] == 1, ">=65", "18-64"),
      vac_lab  = ifelse(g["vac_2plus"] == 1, ">=2 doses", "<2 doses"),
      t_treat_since_inf = round(t_treat_ref + tau_pop, 1)
    )
    k <- k + 1
  }
}
fig4c <- do.call(rbind, Fit4c)
fig4c$treated <- factor(fig4c$treated, levels = c(0, 1))
fig4c$age_lab <- factor(fig4c$age_lab, levels = c("18-64", ">=65"))
fig4c$vac_lab <- factor(fig4c$vac_lab, levels = c("<2 doses", ">=2 doses"))

p4c <- ggplot(fig4c, aes(time, best_fit, group = treated)) +
  geom_ribbon(data = filter(fig4c, treated == "1"),
              aes(ymin = Min95, ymax = Max95, fill = treated), alpha = .1) +
  geom_ribbon(data = filter(fig4c, treated == "1"),
              aes(ymin = Min50, ymax = Max50, fill = treated), alpha = .2) +
  geom_line(aes(colour = treated), linewidth = .75, alpha = .8) +
  geom_vline(aes(xintercept = t_treat_since_inf), colour = "green",
             linetype = "dotted", linewidth = 1) +
  geom_vline(aes(xintercept = t_treat_since_inf + treat_days), colour = "green",
             linetype = "dotted", linewidth = 1) +
  geom_hline(yintercept = LOD, colour = "grey", linetype = "dotted", linewidth = 1) +
  facet_grid(age_lab ~ vac_lab) +
  scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#d62728"),
                      labels = c("Untreated", "Treated")) +
  scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#d62728"), guide = "none") +
  labs(x = "Days after infection",
       y = expression(log[10]~"viral load (RNA copies/mL)"), colour = NULL) +
  coord_cartesian(xlim = c(0, 30), ylim = c(-2, 11), expand = FALSE) +
  theme_vd()

###############################################################################
# Combine into one Figure 4 (A/B/C)
###############################################################################

fig4_combined <- (p4a / p4b / p4c) +
  plot_layout(heights = c(1, 1, 1.6)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(fig_dir, "Fig4_population_trajectories.pdf"), fig4_combined,
       width = 9, height = 16)

Est         <- read.csv(ebe_file, stringsAsFactors = FALSE)
treated_ids <- Est$id[Est$antiviral == 1]

ind_metrics <- data.frame(id = treated_ids, eps_ave = NA_real_,
                          shed_treated = NA_real_, shed_untreated = NA_real_,
                          auc_treated = NA_real_, auc_untreated = NA_real_,
                          peak_treated = NA_real_, peak_untreated = NA_real_)

for (i in seq_along(treated_ids)) {
  r <- Est[Est$id == treated_ids[i], ]

  beta  <- 10^r$log10beta_mode
  pi    <- 10^r$log10pi_mode
  phi   <- 10^r$log10phi_mode
  rho   <- 10^r$log10rho_mode
  delta <- r$delta_mode
  tau   <- r$tau_mode
  EC50  <- r$IC50_mode

  on  <- simulate_viral_load(beta, pi, phi, rho, delta, tau, EC50, t_treat = r$T_treat, treated = 1)
  off <- simulate_viral_load(beta, pi, phi, rho, delta, tau, EC50, t_treat = r$T_treat, treated = 0)

  ind_metrics$eps_ave[i]        <- compute_eps_ave(on, round(r$T_treat + tau, 1))
  ind_metrics$shed_treated[i]   <- shedding_duration(on)
  ind_metrics$shed_untreated[i] <- shedding_duration(off)
  ind_metrics$auc_treated[i]    <- auc_above_lod(on)
  ind_metrics$auc_untreated[i]  <- auc_above_lod(off)
  ind_metrics$peak_treated[i]   <- max(on$V_log10)
  ind_metrics$peak_untreated[i] <- max(off$V_log10)
}

ind_metrics$shed_reduction <- ind_metrics$shed_untreated - ind_metrics$shed_treated
ind_metrics$auc_reduction  <- ind_metrics$auc_untreated  - ind_metrics$auc_treated
ind_metrics$peak_reduction <- ind_metrics$peak_untreated - ind_metrics$peak_treated

write.csv(ind_metrics, file.path(out_dir, "ind_metrics.csv"), row.names = FALSE)

# ---- Supplementary Figures 3-4: individual fits ---------------------------
# Best fit from each person's mode estimates; 95% band from Monolix's
# simulated individual parameters, with tau held at its mode.
sim_pars  <- read.csv(sim_file, stringsAsFactors = FALSE)
min_draws <- min(table(sim_pars$id))
if (min_draws < 20) stop("Only ", min_draws, " simulated draws per person in ", sim_file,
                         " - too few for a 95% band; re-export with more replicates.")

has_obs <- file.exists(patient_file)
if (has_obs) {
  observed_vl <- read.csv(patient_file, stringsAsFactors = FALSE)
} else {
  message("Observed viral loads not found (", patient_file, "); individual fits drawn without data points.")
}

Fit_ind <- vector("list", nrow(Est))
for (i in seq_len(nrow(Est))) {
  r <- Est[i, ]

  best <- simulate_viral_load(10^r$log10beta_mode, 10^r$log10pi_mode, 10^r$log10phi_mode,
                              10^r$log10rho_mode, r$delta_mode, r$tau_mode, r$IC50_mode,
                              t_treat = r$T_treat, treated = r$antiviral)

  draws   <- sim_pars[sim_pars$id == r$id, ]
  V_draws <- sapply(seq_len(nrow(draws)), function(j) {
    simulate_viral_load(10^draws$log10beta[j], 10^draws$log10pi[j], 10^draws$log10phi[j],
                        10^draws$log10rho[j], draws$delta[j], r$tau_mode, draws$IC50[j],
                        t_treat = r$T_treat, treated = r$antiviral)$V_log10
  })
  Min95 <- pmax(apply(V_draws, 1, quantile, 0.025, na.rm = TRUE), -2)
  Max95 <- apply(V_draws, 1, quantile, 0.975, na.rm = TRUE)

  # fit is NA (not floored) below -2, so geom_line() stops drawing there
  # instead of showing a flat clamped segment; Min95/Max95 stay floored,
  # since the ribbon's lower edge is meant to sit at the y-axis cutoff.
  fit_i <- data.frame(id = r$id, time = best$time,
                      fit = ifelse(best$V_log10 < -2, NA_real_, best$V_log10),
                      Min95 = Min95, Max95 = Max95, VL = NA_real_)
  if (has_obs) {
    obs      <- observed_vl[observed_vl$id == r$id, ]
    obs_time <- round(obs$Day_since_detection + round(r$tau_mode, 1), 1)
    fit_i$VL <- obs$VL[match(fit_i$time, obs_time)]
  }
  fit_i$treated <- r$antiviral
  fit_i$t_treat_since_inf <- if (r$antiviral == 1) round(r$T_treat + r$tau_mode, 1) else NA_real_

  Fit_ind[[i]] <- fit_i
}
ind_fits <- do.call(rbind, Fit_ind)

plot_panels <- function(d, file, per_page = 120) {
  # Drop the first 0.5 d: with S(0) fixed and IE(0)=IP(0)=0, V briefly
  # decays before new infections ramp up production, producing a spurious
  # dip-then-rise at the very start of every simulated trajectory.
  d <- filter(d, time >= 0.5)
  ids <- unique(d$id)
  pdf(file, 16, 16)
  for (s in seq(1, length(ids), per_page)) {
    chunk <- filter(d, id %in% ids[s:min(s + per_page - 1, length(ids))])
    print(
      ggplot(chunk, aes(time)) +
        geom_ribbon(aes(ymin = Min95, ymax = Max95, fill = "95% CI"), alpha = .7) +
        geom_line(aes(y = fit, colour = "Best fit"), linewidth = .3) +
        geom_hline(aes(yintercept = LOD, colour = "Detection limit"),
                   linetype = "dotted", linewidth = .6) +
        { if (has_obs) geom_point(aes(y = VL, colour = "Observed"), shape = 16, size = 1) } +
        { if (all(!is.na(chunk$t_treat_since_inf)))
            geom_vline(aes(xintercept = t_treat_since_inf, colour = "Treatment start"),
                       linetype = "dotted", linewidth = .5) } +
        facet_wrap(~ id, ncol = 12) +
        scale_fill_manual(values = c("95% CI" = "pink")) +
        scale_colour_manual(values = c("Best fit" = "red", "Observed" = "blue",
                                       "Treatment start" = "green",
                                       "Detection limit" = "grey")) +
        labs(x = "Days after infection", y = expression(log[10]~"viral load"),
             colour = NULL, fill = NULL) +
        coord_cartesian(xlim = c(0, 30), ylim = c(-2, 10), expand = FALSE) +
        theme_vd(9)
    )
  }
  dev.off()
}

plot_panels(filter(ind_fits, treated == 0), file.path(fig_dir, "SFig3_individual_fits_untreated.pdf"))
plot_panels(filter(ind_fits, treated == 1), file.path(fig_dir, "SFig4_individual_fits_treated.pdf"))
