# Cohort construction and propensity-score matching (Supplementary Figure 1).
# Needs the individual-level clinical data, which are not public (see README);
# shown for transparency. Writes the matched cohort, with study IDs, to
# data_private/patient_data.csv, the file used by scripts 01 and 02.
# Run from the repository root.

library(dplyr)
library(MatchIt)   # method = "optimal" also needs optmatch; distance = "rpart" needs rpart

# ---- Paths ----------------------------------------------------------------
raw_file     <- Sys.getenv("RWD_RAW_DATA", unset = file.path("data_private", "combine_df.RData"))
patient_file <- Sys.getenv("RWD_PATIENT_DATA", unset = file.path("data_private", "patient_data.csv"))
dir.create(dirname(patient_file), recursive = TRUE, showWarnings = FALSE)

load(raw_file)   # combine_df: one row per PCR test, with patient-level fields repeated

set.seed(20260811)

# ---- Derive analysis variables --------------------------------------------
cohort <- combine_df %>%
  mutate(
    age_cat = factor(ifelse(age < 65, "18-64", ">=65"), levels = c("18-64", ">=65")),
    vaccination = factor(case_when(vax_doses %in% c(0, 1) ~ "Notfull",
                                   vax_doses == 2 ~ "Full",
                                   vax_doses == 3 ~ "Booster"),
                         levels = c("Notfull", "Full", "Booster")),
    Severity = factor(Severity, levels = c("Mild", "Severe", "Critical", "Death")),
    comorbidity_sum = coalesce(mh_hyp + mh_dm + mh_cvd + mh_haemato + mh_kidney + mh_neurolo +
                               mh_pulmonary + mh_malignancy + mh_liver + mh_rheumatolo + mh_HIV, 0),
    comorbidity_cat = factor(ifelse(comorbidity_sum == 0, "0", ">=1"), levels = c("0", ">=1")),
    # days from screening detection to treatment; 50 marks untreated patients
    T_treat   = as.numeric(difftime(med_antiviral_date, epi_screen_date, units = "days")),
    T_treat   = if_else(is.na(T_treat), 50, T_treat),
    antiviral = if_else(med_antiviral_whe == 0 | is.na(med_antiviral_whe), 0, med_antiviral_whe),
    Group     = factor(ifelse(antiviral == 0, "Control", "Treatment")),
    # severity at the reference date (treatment start, or day 4 for controls)
    ref_date = if_else(antiviral == 0, as.Date(epi_screen_date) + 4, med_antiviral_date),
    ref_severity = case_when(
      (pgs_oxy == 0) | (pgs_oxy_date - ref_date > 0) ~ "Mild",
      (pgs_oxy == 1 & (pgs_oxy_date - ref_date <= 0)) &
        (pgs_cri == 0 | (pgs_cri == 1 & (pgs_cri_date - ref_date > 0))) ~ "Severe",
      (pgs_cri == 1 & (pgs_cri_date - ref_date <= 0)) &
        (pgs_death == 0 | (pgs_death == 1 & (pgs_death_date - ref_date < 0))) ~ "Critical",
      TRUE ~ "Death"),
    # viral load (log10 copies/mL) from N-gene Ct; Ct 40 (VL 1.31) = below detection
    PCR_date  = as.POSIXct(PCR_date, origin = "1970-01-01", tz = "UTC"),
    PCR_N_.Ct = pmin(PCR_N_.Ct, 40),
    VL        = -0.32 * as.numeric(PCR_N_.Ct) + 14.11,
    censor    = if_else(PCR_N_.Ct == 40, 1, 0),
    Day_since_detection = as.numeric(difftime(PCR_date, epi_screen_date, units = "days"))) %>%
  group_by(hosp_id) %>%
  mutate(init_VL = first(VL[order(Day_since_detection)])) %>%   # first measured viral load
  ungroup() %>%
  # spread same-day replicate tests across the day
  group_by(hosp_id, PCR_date) %>%
  mutate(Day_since_detection = case_when(
    n() == 2 & PCR_seq == max(PCR_seq)    ~ Day_since_detection + 0.5,
    n() == 3 & PCR_seq == median(PCR_seq) ~ Day_since_detection + 0.3,
    n() == 3 & PCR_seq == max(PCR_seq)    ~ Day_since_detection + 0.7,
    TRUE ~ Day_since_detection)) %>%
  ungroup() %>%
  # eligibility: valid viral load; no severe kidney or liver disease; complete
  # screening and outcome dates; treated within 7 days; mild at the reference date
  filter(!is.na(VL),
         !(mh_kidney %in% c(2, 3)) & !(mh_liver == 1),
         !(is.na(epi_screen_date) |
             (pgs_oxy == 1 & is.na(pgs_oxy_date)) |
             (pgs_cri == 1 & is.na(pgs_cri_date))),
         T_treat == 50 | T_treat <= 7,
         ref_severity == "Mild") %>%
  select(hosp_id, age, age_cat, sex, vax_doses, vaccination, comorbidity_sum, comorbidity_cat,
         Severity, antiviral, Group, T_treat, init_VL, VL, censor, Day_since_detection)

# ---- 1:1 optimal matching on a classification-tree propensity score -------
patients <- cohort %>% group_by(hosp_id) %>% slice(1) %>% ungroup()

matched <- matchit(antiviral ~ age_cat + sex + vaccination + comorbidity_cat + init_VL,
                   data = patients, distance = "rpart", method = "optimal", ratio = 1)
print(summary(matched))

matched_ids <- match.data(matched)$hosp_id
cat("Matched:", length(matched_ids), "patients\n")

# ---- Matched cohort with study IDs ----------------------------------------
# Study IDs 1..n follow hospital-ID order, matching the IDs in data/.
id_map <- data.frame(hosp_id = sort(unique(matched_ids)))
id_map$id <- seq_len(nrow(id_map))

patient_data <- cohort %>%
  filter(hosp_id %in% matched_ids) %>%
  mutate(Vaccination_cat = ifelse(vax_doses >= 2, "More than 2 doses", "Less than 2 doses"),
         id = id_map$id[match(hosp_id, id_map$hosp_id)]) %>%
  select(id, everything(), -hosp_id) %>%
  arrange(id, Day_since_detection)

write.csv(patient_data, patient_file, row.names = FALSE)
