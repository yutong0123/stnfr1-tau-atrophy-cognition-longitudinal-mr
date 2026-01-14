HRS Main Analysis
# analysis.R (core models only)

library(dplyr); library(tidyr); library(broom)
library(lme4); library(broom.mixed)

dir.create("output", showWarnings = FALSE)

# ---- load prepared data ----
analysis_wide <- readRDS("data/analysis_wide.rds")

# ---- build logistic dataset ----
logit_df <- analysis_wide %>%
  filter(cognition_status == "Normal") %>%
  mutate(incident_impairment = as.integer(
    cognition2018 %in% c("MCI","Dementia") | cognition2020 %in% c("MCI","Dementia")
  ))

# ---- helper: logistic ----
fit_logit <- function(formula, label){
  glm(formula, family = binomial, data = logit_df) |>
    tidy(exponentiate = TRUE, conf.int = TRUE) |>
    mutate(model = label)
}

logit_out <- bind_rows(
  fit_logit(incident_impairment ~ log_TNFR1 + age + gender, "age_sex"),
  fit_logit(incident_impairment ~ log_TNFR1 + age + gender + race_ethnicity +
              diabetes_combined_alt + hbp_combined_alt + hyperlipidemia + APOE_e4_positive,
            "basic_adj"),
  fit_logit(incident_impairment ~ log_TNFR1 + age + gender + depression_case + race_ethnicity +
              education + bmi_category + diabetes_combined_alt + hbp_combined_alt +
              hyperlipidemia + heart_disease + stroke + marital_status + lung_disease +
              cancer + APOE_e4_positive + currently_smoke + current_drink,
            "full_adj")
)

write.csv(logit_out, "output/logistic_or_ci.csv", row.names = FALSE)

# ---- helper: LMM ----
fit_lmm <- function(prefix, name){
  long_df <- analysis_wide %>%
    pivot_longer(starts_with(prefix), names_to="wave", values_to="score") %>%
    mutate(wave = as.numeric(gsub(prefix, "", wave)),
           years = case_when(wave==13~0, wave==14~2, wave==15~4, TRUE~NA_real_))

  mod <- lmer(score ~ log_TNFR1 + age + gender + depression_case + race_ethnicity + education +
                bmi_category + diabetes_combined_alt + hbp_combined_alt + hyperlipidemia +
                heart_disease + stroke + marital_status + lung_disease + cancer + APOE_e4_positive +
                currently_smoke + current_drink + years + log_TNFR1:years + (1 + years | hhidpn),
              data = long_df)

  write.csv(tidy(mod, effects="fixed", conf.int=TRUE),
            paste0("output/lmm_", name, "_fixed.csv"), row.names=FALSE)
}

fit_lmm("memory_score_w", "memory")
fit_lmm("working_memory_w", "working_memory")
fit_lmm("attention_w", "attention")
fit_lmm("global_cognitive_score_w", "global_cognition")

cat("Done\n")

ADNI Main Analysis
The same workflow was repeated for other outcomes/modalities: merge each longitudinal outcome with baseline TNFR1 and covariates by ID, compute follow-up time in years, keep post-baseline observations with ≥2 visits per subject, fit an LMM scale(outcome) ~ time_yrs*tnfr1_z + covariates + (1 + time_yrs | ID) (add ICV for MRI volumes if needed), and export the tnfr1_z and time_yrs:tnfr1_z effects.
# ADNI main analysis (example): TNFR1 (tnfr1_z) -> CSF tau
library(dplyr)
library(lme4)
library(broom.mixed)

dir.create("output", showWarnings = FALSE)

# Assumption: all tables already use ID (no rid)
# csf_long: ID, examdate, examdate_bl, tau
# tnfr1_bl: ID, tnfr1_z
# covar_bl: ID, age_z, sex, apoe4, educ_z
# strata:   ID, group5

dat <- csf_long %>%
  inner_join(tnfr1_bl, by = "ID") %>%
  inner_join(covar_bl, by = "ID") %>%
  inner_join(strata %>% select(ID, group5), by = "ID") %>%
  mutate(time_yrs = as.numeric(difftime(examdate, examdate_bl, units = "days")) / 365.25) %>%
  filter(time_yrs >= 0, !is.na(tau)) %>%
  group_by(ID) %>% filter(n() >= 2) %>% ungroup()

m_tau <- lmer(
  scale(tau) ~ time_yrs * tnfr1_z + age_z + sex + apoe4 + educ_z + (1 + time_yrs | ID),
  data = dat, REML = TRUE
)

tidy(m_tau, effects = "fixed", conf.int = TRUE) %>%
  filter(term %in% c("tnfr1_z", "time_yrs:tnfr1_z")) %>%
  write.csv("output/adni_example_tau_terms.csv", row.names = FALSE)

cat("Done: output/adni_example_tau_terms.csv\n")
