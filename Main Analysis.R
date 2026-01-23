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

# ---- logistic ----
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

# ---- LMM ----
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

#MR Main Analysis
 # run_mr_online.R (one-page, concise & reproducible)
stopifnot(nzchar(Sys.getenv("OPENGWAS_JWT")))

suppressPackageStartupMessages({
  library(ieugwasr); library(TwoSampleMR)
  library(dplyr); library(purrr); library(readr)
})

# ---- CONFIG ----
chr_gene <- 12L; gene_start <- 6437923L; gene_end <- 6451280L; cis_kb <- 500L
win_start <- gene_start - cis_kb*1000L; win_end <- gene_end + cis_kb*1000L
exposure_id <- "ebi-a-GCST90012015"; p_thr <- 5e-8
clump_kb <- 10000; clump_r2 <- 0.001; pop <- "EUR"
outcomes <- c("ubm-b-216","ubm-b-199","ubm-b-217","ubm-b-200","ubm-b-351","ubm-b-384",
              "ubm-b-401","ubm-b-368","ubm-b-196","ubm-b-213")

# ---- HELPERS ----
with_retry <- function(fun, tries=2) {
  for (i in 1:tries) { x <- try(fun(), TRUE); if (!inherits(x,"try-error")) return(x); Sys.sleep(0.3*i) }
  NULL
}

run_one <- function(exp_all, out_id) {
  exp_cis <- exp_all %>%
    mutate(chr_num = suppressWarnings(as.integer(chr.exposure))) %>%
    filter(chr_num==chr_gene, pos.exposure>=win_start, pos.exposure<=win_end)
  if (nrow(exp_cis)==0) return(NULL)

  exp_cis <- with_retry(function() clump_data(exp_cis, clump_kb=clump_kb, clump_r2=clump_r2, pop=pop), 2)
  if (is.null(exp_cis) || nrow(exp_cis)==0) return(NULL)

  out_dat <- with_retry(function() extract_outcome_data(snps=exp_cis$SNP, outcomes=out_id), 2)
  if (is.null(out_dat) || nrow(out_dat)==0) return(NULL)

  dat0 <- harmonise_data(exp_cis, out_dat, action=2) %>% add_metadata()

  steiger_any <- NA; steiger_pmin <- NA
  try({
    dr <- TwoSampleMR::directionality_test(dat0)
    if ("correct_causal_direction" %in% names(dr)) steiger_any <- any(dr$correct_causal_direction, na.rm=TRUE)
    if (is.na(steiger_any) && "steiger_dir" %in% names(dr)) steiger_any <- any(dr$steiger_dir, na.rm=TRUE)
    if ("steiger_pval" %in% names(dr)) { steiger_pmin <- suppressWarnings(min(dr$steiger_pval, na.rm=TRUE)); if (!is.finite(steiger_pmin)) steiger_pmin <- NA_real_ }
  }, silent=TRUE)

  dat <- steiger_filtering(dat0)
  if (!("steiger_dir" %in% names(dat))) return(NULL)
  dat <- dat %>% filter(steiger_dir | is.na(steiger_dir))

  nsnp <- length(unique(dat$SNP)); if (nsnp==0) return(NULL)
  methods <- if (nsnp==1) "mr_wald_ratio" else if (nsnp==2) c("mr_ivw","mr_weighted_median") else c("mr_ivw","mr_weighted_median","mr_egger_regression")

  mr(dat, method_list=methods) %>%
    mutate(model="main", nsnp=nsnp, ci_low=b-1.96*se, ci_high=b+1.96*se,
           steiger_any=steiger_any, steiger_pmin=steiger_pmin)
}

# ---- RUN ----
exp_all <- with_retry(function() extract_instruments(outcomes=exposure_id, p1=p_thr, clump=FALSE), 2)
stopifnot(!is.null(exp_all), nrow(exp_all)>0, all(c("chr.exposure","pos.exposure") %in% names(exp_all)))

res <- map(outcomes, ~ tryCatch(run_one(exp_all, .x), error=function(e) NULL)) %>%
  compact() %>% bind_rows() %>%
  filter(method %in% c("Wald ratio","Inverse variance weighted","Weighted median")) %>%
  select(model, method, nsnp, exposure, outcome, b, se, ci_low, ci_high, pval, steiger_any, steiger_pmin)

write_tsv(res, "main_mr.tsv")
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
                                                                                                                                                                                                                                                                                                                     
