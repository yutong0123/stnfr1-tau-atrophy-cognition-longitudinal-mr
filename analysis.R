# analysis.R (core models only)

library(dplyr); library(tidyr); library(broom)
library(lme4); library(broom.mixed)

dir.create("output", showWarnings = FALSE)

# ---- load prepared data ----
analysis_wide <- readRDS("data/analysis_wide.rds")

# 3) Logistic regression (incident impairment by 2018/2020)
logit_df <- analysis_wide %>%
  filter(cognition_status == "Normal") %>%
  mutate(incident_impairment = ifelse(
    cognition2018 %in% c("MCI", "Dementia") | cognition2020 %in% c("MCI", "Dementia"), 1, 0
  ))

m_age_sex <- glm(incident_impairment ~ log_TNFR1 + age + gender,
                 family = binomial, data = logit_df)

m_basic <- glm(incident_impairment ~ log_TNFR1 + age + gender + race_ethnicity +
                 diabetes_combined_alt + hbp_combined_alt + hyperlipidemia + APOE_e4_positive,
               family = binomial, data = logit_df)

m_full <- glm(incident_impairment ~ log_TNFR1 + age + gender + depression_case + race_ethnicity +
                education + bmi_category + diabetes_combined_alt + hbp_combined_alt +
                hyperlipidemia + heart_disease + stroke + marital_status + lung_disease +
                cancer + APOE_e4_positive + currently_smoke + current_drink,
              family = binomial, data = logit_df)

bind_rows(
  tidy(m_age_sex, exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "age_sex"),
  tidy(m_basic,   exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "basic_adj"),
  tidy(m_full,    exponentiate = TRUE, conf.int = TRUE) %>% mutate(model = "full_adj")
) %>% write.csv("output/logistic_or_ci.csv", row.names = FALSE)

# 4) Longitudinal LMM (4 outcomes)
fit_lmm <- function(prefix, name){
  long_df <- analysis_wide %>%
    pivot_longer(starts_with(prefix), names_to="wave", values_to="score") %>%
    mutate(wave = as.numeric(gsub(prefix, "", wave)),
           years = dplyr::case_when(wave==13~0, wave==14~2, wave==15~4, TRUE~NA_real_))

  mod <- lmer(score ~ log_TNFR1 + age + gender + depression_case + race_ethnicity + education +
                bmi_category + diabetes_combined_alt + hbp_combined_alt + hyperlipidemia +
                heart_disease + stroke + marital_status + lung_disease + cancer + APOE_e4_positive +
                currently_smoke + current_drink + years + log_TNFR1:years + (1 + years | hhidpn),
              data = long_df)

  broom.mixed::tidy(mod, effects="fixed", conf.int=TRUE) %>%
    write.csv(paste0("output/lmm_", name, "_fixed.csv"), row.names=FALSE)
}

fit_lmm("memory_score_w", "memory")
fit_lmm("working_memory_w", "working_memory")
fit_lmm("attention_w", "attention")
fit_lmm("global_cognitive_score_w", "global_cognition")

cat("Done\n")
