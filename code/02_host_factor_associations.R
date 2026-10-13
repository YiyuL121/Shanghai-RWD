# Host-factor associations with viral dynamics and treatment efficacy.
# Writes the population-model covariate tests, Figure 3, Supplementary Table 4
# and Supplementary Figures 2 and 5. Reads ind_metrics.csv from script 01 (or the copy in data/).
# Run from the repository root.

library(dplyr)
library(ggplot2)
library(emmeans)
library(ggtext)     # bold/larger panel letters in Supplementary Figure 5
library(patchwork)  # combine multi-panel figures into one file

# ---- Paths ----------------------------------------------------------------
data_dir     <- file.path("data", "M4_primary")
pop_file     <- file.path(data_dir, "populationParameters.txt")
ebe_file     <- file.path(data_dir, "estimatedIndividualParameters.txt")
# Continuous age and baseline viral load are not public. If available (see
# README), Supplementary Figure 2 and the timing regression are also produced.
patient_file <- Sys.getenv("RWD_PATIENT_DATA", unset = file.path("data_private", "patient_data.csv"))
out_dir      <- file.path("output", "tables")
fig_dir      <- file.path("output", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

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

# ---- Load data ------------------------------------------------------------

pop <- read.csv(pop_file, stringsAsFactors = FALSE)
pop$parameter <- trimws(pop$parameter)
get_est <- function(name) as.numeric(pop$value[pop$parameter == name])
get_se  <- function(name) as.numeric(pop$se_sa[pop$parameter == name])

Est <- read.csv(ebe_file, stringsAsFactors = FALSE)

has_patient <- file.exists(patient_file)
if (has_patient) {
  patient <- read.csv(patient_file, stringsAsFactors = FALSE)
  Est <- merge(Est, patient[!duplicated(patient$id), c("id", "age", "init_VL")], by = "id")
} else {
  message("Patient-level data not found (", patient_file, "); skipping Supplementary ",
          "Figure 2 and the timing/baseline viral load regression.")
}

# Covariate flags exactly as fed into the SAEM model.
Est$age_ge65 <- as.integer(Est$age_cat == ">=65")
Est$vac_ge2  <- as.integer(Est$Vaccination_cat == "More than 2 doses")
Est$age_cat  <- factor(Est$age_cat, levels = c("18-64", ">=65"))

# 2-level vaccination split matching the SAEM covariate.
Est$vac_group <- factor(ifelse(Est$vac_ge2 == 1, ">=2 doses", "<2 doses"),
                        levels = c("<2 doses", ">=2 doses"))
Est$comorbidity_cat <- factor(as.character(Est$comorbidity_cat), levels = c("0", ">=1"))
Est$sex <- factor(Est$sex)

# Eventual severity status: Severity is Mild/Severe/Critical/Death; there is no
# category between Mild and Severe, so "No" = Mild and "Yes" = Severe or worse.
Est$severity_yn <- factor(ifelse(Est$Severity %in% c("Severe", "Critical", "Death"), "Yes", "No"),
                          levels = c("No", "Yes"))

# From script 01 if it has been run; otherwise the copy shipped with the data.
ind_metrics_file <- file.path(out_dir, "ind_metrics.csv")
if (!file.exists(ind_metrics_file)) ind_metrics_file <- file.path(data_dir, "ind_metrics.csv")
ind_metrics <- read.csv(ind_metrics_file, stringsAsFactors = FALSE)
treated <- Est %>% filter(antiviral == 1) %>% left_join(ind_metrics, by = "id")

# ---- Covariate effects estimated in the population model (Wald tests) ---

pop_cov <- data.frame(
  parameter    = c("delta", "delta", "EC50"),
  covariate    = c("Age >=65 (vs 18-64)",
                   "Vaccination >=2 doses (vs <2)",
                   "Vaccination >=2 doses (vs <2)"),
  monolix_name = c("beta_delta_age_cat___65",
                   "beta_delta_Vaccination_cat_More_than_2_doses",
                   "beta_IC50_Vaccination_cat_More_than_2_doses"),
  stringsAsFactors = FALSE
)
pop_cov$estimate <- sapply(pop_cov$monolix_name, get_est)
pop_cov$se       <- sapply(pop_cov$monolix_name, get_se)
pop_cov$z        <- pop_cov$estimate / pop_cov$se
pop_cov$p_value  <- 2 * pnorm(-abs(pop_cov$z))
pop_cov$ratio    <- exp(pop_cov$estimate)
pop_cov$ratio_lo <- exp(pop_cov$estimate - 1.96 * pop_cov$se)
pop_cov$ratio_hi <- exp(pop_cov$estimate + 1.96 * pop_cov$se)
pop_cov$interpretation <- c(
  "older age -> slower infected-cell clearance",
  "more vaccination -> faster infected-cell clearance",
  "more vaccination -> lower EC50 (higher apparent potency)"
)

print(pop_cov %>% mutate(across(c(estimate, se, z, ratio, ratio_lo, ratio_hi), ~ round(.x, 3)),
                        p_value = signif(p_value, 3)))
write.csv(pop_cov, file.path(out_dir, "population_model_covariate_effects.csv"), row.names = FALSE)

# ---- Post-hoc multivariable regression of individual estimates (EBEs) ----

ebe_regression <- function(data, outcome,
                           rhs = c("age_cat", "sex", "comorbidity_cat", "vac_group", "severity_yn")) {
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

Est$log_delta_mode    <- log(Est$delta_mode)
treated$log_IC50_mode <- log(treated$IC50_mode)

all_arm_outcomes     <- c("log10beta_mode", "log10pi_mode", "log10phi_mode", "log10rho_mode",
                          "tau_mode", "log_delta_mode")
treated_arm_outcomes <- c("log_IC50_mode", "eps_ave")

ebe_tab <- bind_rows(
  bind_rows(lapply(all_arm_outcomes, function(out) ebe_regression(Est, out) %>% mutate(sample = "all"))),
  bind_rows(lapply(treated_arm_outcomes, function(out) ebe_regression(treated, out) %>% mutate(sample = "treated")))
)

print(ebe_tab %>% mutate(across(c(estimate, se, lo, hi), ~ round(.x, 3)), p_value = signif(p_value, 3)))

# ---- Figure 3: clearance and efficacy by host factor ----------------------

group_panel <- function(data, xvar, yvar, xlab, ylab, model_rhs) {
  m   <- lm(reformulate(model_rhs, response = yvar), data = data)
  emm <- as.data.frame(emmeans(m, specs = xvar))
  names(emm)[1] <- "x"
  d <- data %>% mutate(x = .data[[xvar]], y = .data[[yvar]])

  ggplot(d, aes(x, y)) +
    geom_boxplot(outlier.shape = NA, fill = "grey90", colour = "grey40") +
    geom_point(position = position_jitter(width = .1, height = 0, seed = 20260811), alpha = .25, colour = "darkred") +
    geom_line(data = emm, aes(x, emmean, group = 1), colour = "blue",
              linetype = "dashed", linewidth = .7) +
    geom_errorbar(data = emm, aes(x, y = emmean, ymin = lower.CL, ymax = upper.CL),
                  width = .3, colour = "blue", linewidth = 1) +
    geom_point(data = emm, aes(x, emmean), colour = "blue", size = 3) +
    labs(x = xlab, y = ylab) +
    theme_vd()
}

rhs_treated <- c("age_cat", "sex", "comorbidity_cat", "vac_group", "severity_yn")

p3a <- group_panel(Est, "vac_group", "delta_mode", "Vaccination group",
                   expression("Infected cell clearance rate"~delta~(d^-1)), rhs_treated)
p3b <- group_panel(Est, "age_cat", "delta_mode", "Age group",
                   expression("Infected cell clearance rate"~delta~(d^-1)), rhs_treated)
p3c <- group_panel(treated, "vac_group", "eps_ave", "Vaccination group",
                   "Average treatment efficacy", rhs_treated)
p3d <- group_panel(treated, "age_cat", "eps_ave", "Age group",
                   "Average treatment efficacy", rhs_treated)

fig3_combined <- (p3a + p3b) / (p3c + p3d) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(fig_dir, "Fig3_delta_efficacy_by_hostfactor.pdf"), fig3_combined,
       width = 10, height = 9,
       device = if (capabilities("cairo")) cairo_pdf else pdf)   # cairo_pdf embeds the font for the Greek delta

# ---- Supplementary Figure 5: individual parameter distributions -----------
panel_letter <- function(letter, text) {
  sprintf("<span style='font-size:15pt; font-weight:bold;'>%s.</span> %s", letter, text)
}
sfig5_labels <- c(
  log10beta_mode = panel_letter("A", "Infectivity (log10 beta)"),
  log10pi_mode   = panel_letter("B", "Production rate (log10 pi)"),
  log10rho_mode  = panel_letter("C", "Reversion rate (log10 rho)"),
  log10phi_mode  = panel_letter("D", "Conversion rate (log10 phi)"),
  delta_mode     = panel_letter("E", "Clearance rate (delta)"),
  IC50_mode      = panel_letter("F", "EC50 (treated arm only)")
)
# EC50 is identified only for the treated arm (a control's efficacy term is
# zero throughout), so restrict its histogram to treated individuals - pooling
# all individuals produced two prior-value spikes instead of heterogeneity.
sfig5_dat <- do.call(rbind, lapply(names(sfig5_labels), function(col) {
  vals <- if (col == "IC50_mode") Est[[col]][Est$antiviral == 1] else Est[[col]]
  data.frame(parameter = sfig5_labels[[col]], value = vals)
}))
sfig5_dat$parameter <- factor(sfig5_dat$parameter, levels = unname(sfig5_labels))

p_sfig5 <- ggplot(sfig5_dat, aes(value)) +
  geom_histogram(bins = 30, fill = "skyblue", colour = "black") +
  facet_wrap(~ parameter, scales = "free_x", ncol = 3) +   # x free, y fixed
  labs(x = NULL, y = "Count") +
  theme_vd(12) +
  theme(panel.spacing = unit(1, "cm"),
        strip.text = ggtext::element_markdown(size = 12))
ggsave(file.path(fig_dir, "SFig5_parameter_distributions.pdf"), p_sfig5, width = 11, height = 7)

# ---- Supplementary Figure 2: age x vaccination x efficacy (needs patient data)

if (has_patient) {
  p_sfig2 <- ggplot(treated, aes(age, eps_ave, colour = vac_group)) +
    geom_point() +
    facet_wrap(~ vac_group) +
    labs(x = "Age (years)", y = "Average efficacy", colour = "Vaccination status") +
    theme_vd()
  ggsave(file.path(fig_dir, "SFig2_age_vaccination_efficacy.pdf"), p_sfig2, width = 10, height = 5)
}

# ---- Supplementary Table 4: host factor x parameter ----------------------

fmt_cell <- function(est, lo, hi, p) {
  stars <- ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
  sprintf("%.2f%s (%.2f, %.2f)", est, stars, lo, hi)
}
fmt_ratio_cell <- function(est, se, p) fmt_cell(exp(est), exp(est - 1.96 * se), exp(est + 1.96 * se), p)

param_cols <- c("log10beta_mode", "log10pi_mode", "log10rho_mode", "log10phi_mode",
                "log_delta_mode", "tau_mode", "log_IC50_mode", "eps_ave")
ratio_cols <- c("log_delta_mode", "log_IC50_mode")

contrasts <- c(
  "Age >=65 (vs 18-64)"                    = "age_cat>=65",
  "Sex male (vs female)"                   = "sexMale",
  "Comorbidities >=1 (vs 0)"               = "comorbidity_cat>=1",
  "Vaccination >=2 doses (vs <2 doses)"    = "vac_group>=2 doses",
  "Eventual severity status: Yes (vs No)"  = "severity_ynYes"
)

saem_override <- list(
  "log_delta_mode|Age >=65 (vs 18-64)"                 = pop_cov[pop_cov$monolix_name == "beta_delta_age_cat___65", ],
  "log_delta_mode|Vaccination >=2 doses (vs <2 doses)" = pop_cov[pop_cov$monolix_name == "beta_delta_Vaccination_cat_More_than_2_doses", ],
  "log_IC50_mode|Vaccination >=2 doses (vs <2 doses)"  = pop_cov[pop_cov$monolix_name == "beta_IC50_Vaccination_cat_More_than_2_doses", ]
)

supp_table4 <- data.frame(Covariate = names(contrasts), stringsAsFactors = FALSE)
for (col in param_cols) {
  supp_table4[[col]] <- sapply(names(contrasts), function(cn) {
    key <- paste(col, cn, sep = "|")
    # Three cells (delta~age, delta~vaccination, EC50~vaccination) use the
    # population-model Wald estimate; every other cell is the post-hoc EBE
    # regression. eps_ave is a near-deterministic transform of each patient's
    # EC50 (population PK is fixed), so it is not independent of the EC50 column.
    if (key %in% names(saem_override)) {
      row <- saem_override[[key]]
      return(fmt_ratio_cell(row$estimate, row$se, row$p_value))
    }
    r <- ebe_tab[ebe_tab$outcome == col & ebe_tab$term == contrasts[[cn]], ]
    if (nrow(r) == 0) return(NA_character_)
    if (col %in% ratio_cols) fmt_ratio_cell(r$estimate, r$se, r$p_value)
    else fmt_cell(r$estimate, r$lo, r$hi, r$p_value)
  })
}

names(supp_table4) <- c("Covariate", "log10beta", "log10pi", "log10rho", "log10phi",
                        "delta (ratio)", "tau (days)", "EC50 (ratio)", "eps_ave")

print(supp_table4, row.names = FALSE)
write.csv(supp_table4, file.path(out_dir, "SupplementaryTable4_associations.csv"), row.names = FALSE)

# ---- Timing and baseline viral load vs average treatment efficacy ---------
# Needs baseline viral load from the patient-level data.

if (has_patient) {
  timing_model <- lm(eps_ave ~ age_cat + sex + comorbidity_cat + vac_group + severity_yn +
                       T_treat + init_VL, data = treated)
  timing_s  <- summary(timing_model)$coefficients
  timing_ci <- confint(timing_model)

  timing_out <- data.frame(
    term     = c("T_treat", "init_VL"),
    estimate = timing_s[c("T_treat", "init_VL"), "Estimate"],
    lo       = timing_ci[c("T_treat", "init_VL"), 1],
    hi       = timing_ci[c("T_treat", "init_VL"), 2],
    p_value  = timing_s[c("T_treat", "init_VL"), "Pr(>|t|)"],
    row.names = NULL
  )
  print(timing_out)
  write.csv(timing_out, file.path(out_dir, "timing_baseline_vl_vs_efficacy.csv"), row.names = FALSE)
}
