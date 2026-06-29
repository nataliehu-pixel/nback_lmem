# Created on Sat Apr 11 09:28:39 2026
# erinanncarpenter, nataliehu

rm(list=ls(all=TRUE))
# ---- Libraries ---- #
# install.packages("lme4")
# install.packages("tidyverse")
# install.packages("gamlss")
library(lme4)
library(tidyverse)
library(ggplot2)
library(gamlss)
library(emmeans)
library(lmerTest)

# ---- Setup ---- #
# Set path
code_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
root_dir <- dirname(code_dir)
derivatives_dir <- file.path(root_dir, "derivatives/nback-v3.0")
dir.create(derivatives_dir, showWarnings = FALSE, recursive = TRUE)

# Load processed data
fname <- "group_task-nback_desc-proc_beh.csv"
data_path <- file.path(derivatives_dir, fname)
df <- read_csv(data_path)

# ---- Helper: Prepare Cohort Data ---- #
# Shared columns needed for accuracy models
cols_needed <- c("response_accuracy", "task_cond_gen", "language_group", "age", 
                 "edu", "subject", "trial", "session")

prepare_cohort <- function(df, cohort_filter) {
  df_cohort <- df %>%
    filter(.data$cohort == cohort_filter) %>%
    drop_na(all_of(cols_needed)) %>%
    mutate(
      task_cond_gen = relevel(
        factor(task_cond_gen, levels = c("0back", "1back", "2back")),
        ref = "0back"
      ),
      language_group = relevel(
        factor(language_group, levels = c("bilingual", "monolingual")),
        ref = "monolingual"
      ),
      session = relevel(
        factor(session, levels = c("1", "2")),
        ref = "1"
      ),
      age_c = age - mean(age, na.rm = TRUE),
      edu_c = edu - mean(edu, na.rm = TRUE),
      mpo_c = mpo - mean(mpo, na.rm =  TRUE)
    )
  
  return(df_cohort)
}

# ---- Helper: Clean RT Data ---- #
clean_rt <- function(df_cohort) {
  df_rt <- df_cohort %>%
    filter(!is.na(response_time_ms)) %>%
    filter(response_time_ms >= 200) %>%
    group_by(subject, task_cond_gen) %>%
    mutate(
      rt_mean   = mean(response_time_ms),
      rt_sd     = sd(response_time_ms),
      rt_cutoff = rt_mean + 3 * rt_sd
    ) %>%
    filter(response_time_ms <= rt_cutoff) %>%
    ungroup() %>%
    select(-rt_mean, -rt_sd, -rt_cutoff) %>%
    mutate(log_rt = log(response_time_ms))
  
  return(df_rt)
}

# ---- Helper: Extract Ex-Gaussian Parameters ---- #
# Fits an ex-Gaussian distribution per subject x condition using gamlss
# and returns mu (normal mean), sigma (normal SD), and tau (exponential component)

# load cleaned rt data:


fname <- "group_task-nback_correct-trials-only_beh.csv"
data_path <- file.path(derivatives_dir, fname)
df <- read_csv(data_path)

extract_exgaussian <- function(df_rt) {
  df_rt %>%
    group_by(subject, task_cond_gen, language_group, age, age_c, mpo, mpo_c, edu, edu_c, session) %>%
    summarise(
      exg = list(withCallingHandlers(
        tryCatch({
          fit <- gamlss(
            response_time_ms ~ 1,
            family = exGAUS(),
            trace  = FALSE
          )
          data.frame(
            mu    = exp(fit$mu.coefficients),
            sigma = exp(fit$sigma.coefficients),
            tau   = exp(fit$nu.coefficients)
          )
        }, error = function(e) data.frame(mu = NA, sigma = NA, tau = NA)),
        warning = function(w) {
          if (grepl("not yet converged", conditionMessage(w))) {
            invokeRestart("muffleWarning")
          }
        }
      )),
      .groups = "drop"
    ) %>%
    unnest(exg) %>%
    filter(if_all(c(mu, sigma, tau), ~ !is.na(.) & is.finite(.)))
}




# ---- Shared Plot Settings ---- #
conditions     <- c("0back", "1back", "2back")
language_groups <- c("bilingual", "monolingual")
title_map      <- c("0back" = "0-back", "1back" = "1-back", "2back" = "2-back")
colors         <- c("monolingual" = "steelblue", "bilingual" = "coral")

# ===========================================================================
# ---- PSA Models: Monolingual vs. Bilingual ---- #
# ===========================================================================
df_psa    <- prepare_cohort(df, "PSA")
df_psa_rt <- clean_rt(df_psa)

# ---- PSA Accuracy Model ---- #
model_psa_acc <- glmer(
  response_accuracy ~ task_cond_gen + language_group + age_c + 
    task_cond_gen:language_group + task_cond_gen:age_c  + (1 | subject),
  data   = df_psa,
  family = binomial
)
print(summary(model_psa_acc))
saveRDS(model_psa_acc, file.path(derivatives_dir, "psa_accuracy_model.rds"))

# ---- PSA RT Model (Log-transformed) ---- #
# Maximal random effects structure (Barr et al., 2013): start with
# random intercept + slopes for condition and session, no correlations (||)
# If convergence fails, reduce by dropping session slope, then condition slope

model_psa_rt <- lmer(
  log_rt ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen + session || subject),
  data = df_psa_rt
)

# If singular or non-convergence, fall back step-by-step:
model_psa_rt <- lmer(
  log_rt ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen || subject),
  data = df_psa_rt
)


model_psa_rt <- lmer(
  log_rt ~ task_cond_gen * language_group * age_c + session * stimulus_type + edu_c + mpo_c  + (1  | subject),  
  data = df_psa_rt
)

print(summary(model_psa_rt))
saveRDS(model_psa_rt, file.path(derivatives_dir, "psa_rt_logrt_model.rds"))

# model estimated maarginal means: 
psa_emms = emmeans(model_psa_rt, 
        ~ task_cond_gen * language_group,
        pbkrtest.limit = 5000)

emm_pwc = contrast(psa_emms, method='pairwise') 

emm_pwc_df = as.data.frame(emm_pwc)
View(emm_pwc_df)


dplyr::last_dplyr_warnings()

trial_counts = df_psa_rt %>%
  count(subject, session, task_cond_gen) %>%
  arrange(n)

# ----- EXPORT RT LMER RESULTS -----

library(dplyr)
library(broom.mixed)
library(flextable)
library(officer)


# Format parameter names

format_param <- function(param) {
  param <- gsub("\\(Intercept\\)", "Intercept", param)
  
  param <- gsub("task_cond_gen1back", "1-back vs. 0-back", param)
  param <- gsub("task_cond_gen2back", "2-back vs. 0-back", param)
  
  param <- gsub("language_groupmonolingual", "Monolingual vs. bilingual", param)
  
  param <- gsub("age_c", "Age", param)
  param <- gsub("edu_c", "Education", param)
  
  param <- gsub("stimulus_typestar", "Stimulus type: star", param)
  
  param <- gsub(":", " × ", param)
  
  return(param)
}

fix_minus <- function(x) {
  gsub("-", "\u2212", x)
}


# Extract fixed effects


final_df <- broom.mixed::tidy(
  model_psa_rt,
  effects = "fixed",
  conf.int = TRUE,
  conf.level = 0.95
) %>%
  mutate(
    Parameter = format_param(term),
    p_value = p.value
  ) %>%
  select(
    Parameter,
    Estimate = estimate,
    SE = std.error,
    df,
    t = statistic,
    p_value,
    CI_low = conf.low,
    CI_high = conf.high
  )


# Format table


table_formatted <- final_df %>%
  mutate(
    Estimate = fix_minus(sprintf("%.3f", Estimate)),
    SE = fix_minus(sprintf("%.3f", SE)),
    df = sprintf("%.1f", df),
    t = fix_minus(sprintf("%.2f", t)),
    p_value = case_when(
      is.na(p_value) ~ "—",
      p_value < .001 ~ "< .001",
      TRUE ~ sprintf("%.3f", p_value)
    ),
    `95% CI` = paste0(
      "[",
      fix_minus(sprintf("%.3f", CI_low)),
      ", ",
      fix_minus(sprintf("%.3f", CI_high)),
      "]"
    )
  ) %>%
  select(
    Parameter,
    `β` = Estimate,
    SE,
    df,
    t,
    p = p_value,
    `95% CI`
  )

# Define output path


out_dir <- file.path(derivatives_dir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "psa_rt_lmer_output_table.docx")


# Export flextable


ft <- flextable(table_formatted) %>%
  autofit() %>%
  theme_booktabs() %>%
  fontsize(size = 10, part = "all") %>%
  flextable::font(fontname = "Times New Roman", part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = "Parameter", align = "left", part = "body")

doc <- read_docx()

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      "Table X.",
      prop = fp_text(
        font.family = "Times New Roman",
        bold = TRUE,
        font.size = 11
      )
    )
  ),
  style = "Normal"
)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      "Linear mixed-effects model predicting log-transformed reaction time in the post-stroke aphasia group.",
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 11
      )
    )
  ),
  style = "Normal"
)

doc <- body_add_flextable(doc, ft)

doc <- body_add_fpar(
  doc,
  fpar(
    ftext(
      "Note.",
      prop = fp_text(
        font.family = "Times New Roman",
        italic = TRUE,
        font.size = 11
      )
    ),
    ftext(
      " β = fixed-effect estimate; SE = standard error; df = degrees of freedom; CI = confidence interval. The outcome variable was log-transformed reaction time. Degrees of freedom and p-values were estimated using the Satterthwaite approximation when the model was fit with lmerTest.",
      prop = fp_text(
        font.family = "Times New Roman",
        font.size = 11
      )
    )
  ),
  style = "Normal"
)

doc <- body_end_section_landscape(doc)

print(doc, target = out_file)

# ---- PSA RT Model (Ex-Gaussian) ---- 
# Extract mu, sigma, tau per subject x condition, then model each separately
df_psa_exg <- extract_exgaussian(df_psa)

model_psa_rt_mu <- lm(
  mu ~ task_cond_gen * age_c + language_group,
  data = df_psa_exg
)

summary(model_psa_rt_mu)


model_psa_rt_sigma <- lmer(
  sigma ~ task_cond_gen * language_group * age_c + edu_c +
    (1 | subject),
  data = df_psa_exg
)
model_psa_rt_tau <- lmer(
  tau ~ task_cond_gen * language_group * age_c + edu_c +
    (1 | subject),
  data = df_psa_exg
)

print(tryCatch(summary(model_psa_rt_mu),    error = function(e) paste("mu model summary failed:", conditionMessage(e))))
print(tryCatch(summary(model_psa_rt_sigma), error = function(e) paste("sigma model summary failed:", conditionMessage(e))))
print(tryCatch(summary(model_psa_rt_tau),   error = function(e) paste("tau model summary failed:", conditionMessage(e))))

saveRDS(model_psa_rt_mu,    file.path(derivatives_dir, "psa_rt_exg_mu_model.rds"))
saveRDS(model_psa_rt_sigma, file.path(derivatives_dir, "psa_rt_exg_sigma_model.rds"))
saveRDS(model_psa_rt_tau,   file.path(derivatives_dir, "psa_rt_exg_tau_model.rds"))


t.test(mu ~ language_group, data = df_psa_exg)

cor(df_psa_exg$age, df_psa_exg$mu, method="pearson")
plot( df_psa_exg$age, df_psa_exg$mu)

# ---- PSA Prediction Grids ---- #
age_range_psa <- seq(min(df_psa$age_c), max(df_psa$age_c), length.out = 100)

make_pred_grid <- function(age_range, age_mean) {
  grid <- expand.grid(
    task_cond_gen  = conditions,
    language_group = language_groups,
    age_c          = age_range,
    edu_c          = 0
  )
  grid$task_cond_gen  <- factor(grid$task_cond_gen,  levels = c("0back", "1back", "2back"))
  grid$language_group <- factor(grid$language_group, levels = c("bilingual", "monolingual"))
  grid$age            <- grid$age_c + age_mean
  return(grid)
}

pred_grid_psa_acc <- make_pred_grid(age_range_psa, mean(df_psa$age))
pred_grid_psa_rt  <- make_pred_grid(age_range_psa, mean(df_psa$age))

# ---- PSA Accuracy Visualization ---- #
pred_grid_psa_acc$predicted <- predict(model_psa_acc, newdata = pred_grid_psa_acc, type = "response", re.form = NA)

subject_acc_psa <- df_psa %>%
  group_by(subject, age, language_group, task_cond_gen) %>%
  summarise(mean_acc = mean(response_accuracy, na.rm = TRUE), .groups = "drop")

p_psa_acc <- ggplot(pred_grid_psa_acc, aes(x = age, y = predicted, color = language_group)) +
  geom_point(data = subject_acc_psa, aes(x = age, y = mean_acc), alpha = 0.5, size = 1.5) +
  geom_line() +
  facet_wrap(~ task_cond_gen, labeller = labeller(task_cond_gen = title_map)) +
  scale_color_manual(values = colors, labels = c("bilingual" = "Bilingual", "monolingual" = "Monolingual")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Predicted Accuracy by Age and Language Group",
    subtitle = "PSA Cohort — Split by Condition",
    x = "Age", y = "Predicted Accuracy", color = "Language Group"
  ) +
  theme_minimal()

print(p_psa_acc)
ggsave(file.path(derivatives_dir, "psa_accuracy_lmem.png"), plot = p_psa_acc, width = 15, height = 5, dpi = 300)

# ---- PSA RT Visualization ---- #
pred_grid_psa_rt$predicted_rt <- exp(predict(model_psa_rt, newdata = pred_grid_psa_rt, re.form = NA))

subject_rt_psa <- df_psa_rt %>%
  group_by(subject, age, language_group, task_cond_gen) %>%
  summarise(mean_rt = mean(response_time_ms, na.rm = TRUE), .groups = "drop")

p_psa_rt <- ggplot(pred_grid_psa_rt, aes(x = age, y = predicted_rt, color = language_group)) +
  geom_point(data = subject_rt_psa, aes(x = age, y = mean_rt), alpha = 0.5, size = 1.5) +
  geom_line() +
  facet_wrap(~ task_cond_gen, labeller = labeller(task_cond_gen = title_map)) +
  scale_color_manual(values = colors, labels = c("bilingual" = "Bilingual", "monolingual" = "Monolingual")) +
  labs(
    title    = "Predicted Response Time by Age and Language Group",
    subtitle = "PSA Cohort — Split by Condition",
    x = "Age", y = "Predicted Response Time (ms)", color = "Language Group"
  ) +
  theme_minimal()

print(p_psa_rt)
ggsave(file.path(derivatives_dir, "psa_rt_lmem.png"), plot = p_psa_rt, width = 15, height = 5, dpi = 300)

# ===========================================================================
# ---- CN Models: Monolingual vs. Bilingual ---- #
# ===========================================================================
df_cn    <- prepare_cohort(df, "CN")
df_cn_rt <- clean_rt(df_cn)

# ---- CN Accuracy Model ---- #
model_cn_acc <- glmer(
  response_accuracy ~ task_cond_gen * language_group * age_c + edu_c + (1 | subject),
  data   = df_cn,
  family = binomial
)
print(summary(model_cn_acc))
saveRDS(model_cn_acc, file.path(derivatives_dir, "cn_accuracy_model.rds"))

# ---- CN RT Model (Log-transformed) ---- #
model_cn_rt <- lmer(
  log_rt ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen + session || subject),
  data = df_cn_rt
)

print(summary(model_cn_rt))
saveRDS(model_cn_rt, file.path(derivatives_dir, "cn_rt_logrt_model.rds"))

# ---- CN RT Model (Ex-Gaussian) ---- #
df_cn_exg <- extract_exgaussian(df_cn_rt)

model_cn_rt_mu <- lmer(
  mu ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen || subject),
  data = df_cn_exg
)
model_cn_rt_sigma <- lmer(
  sigma ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen || subject),
  data = df_cn_exg
)
model_cn_rt_tau <- lmer(
  tau ~ task_cond_gen * language_group * age_c + edu_c +
    (1 + task_cond_gen || subject),
  data = df_cn_exg
)

print(summary(model_cn_rt_mu))
print(summary(model_cn_rt_sigma))
print(summary(model_cn_rt_tau))

saveRDS(model_cn_rt_mu,    file.path(derivatives_dir, "cn_rt_exg_mu_model.rds"))
saveRDS(model_cn_rt_sigma, file.path(derivatives_dir, "cn_rt_exg_sigma_model.rds"))
saveRDS(model_cn_rt_tau,   file.path(derivatives_dir, "cn_rt_exg_tau_model.rds"))

# ---- CN Prediction Grids ---- #
age_range_cn <- seq(min(df_cn$age_c), max(df_cn$age_c), length.out = 100)

pred_grid_cn_acc <- make_pred_grid(age_range_cn, mean(df_cn$age))
pred_grid_cn_rt  <- make_pred_grid(age_range_cn, mean(df_cn$age))

# ---- CN Accuracy Visualization ---- #
pred_grid_cn_acc$predicted <- predict(model_cn_acc, newdata = pred_grid_cn_acc, type = "response", re.form = NA)

subject_acc_cn <- df_cn %>%
  group_by(subject, age, language_group, task_cond_gen) %>%
  summarise(mean_acc = mean(response_accuracy, na.rm = TRUE), .groups = "drop")

p_cn_acc <- ggplot(pred_grid_cn_acc, aes(x = age, y = predicted, color = language_group)) +
  geom_point(data = subject_acc_cn, aes(x = age, y = mean_acc), alpha = 0.5, size = 1.5) +
  geom_line() +
  facet_wrap(~ task_cond_gen, labeller = labeller(task_cond_gen = title_map)) +
  scale_color_manual(values = colors, labels = c("bilingual" = "Bilingual", "monolingual" = "Monolingual")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Predicted Accuracy by Age and Language Group",
    subtitle = "CN Cohort — Split by Condition",
    x = "Age", y = "Predicted Accuracy", color = "Language Group"
  ) +
  theme_minimal()

print(p_cn_acc)
ggsave(file.path(derivatives_dir, "cn_accuracy_lmem.png"), plot = p_cn_acc, width = 15, height = 5, dpi = 300)

# ---- CN RT Visualization ---- #
pred_grid_cn_rt$predicted_rt <- exp(predict(model_cn_rt, newdata = pred_grid_cn_rt, re.form = NA))

subject_rt_cn <- df_cn_rt %>%
  group_by(subject, age, language_group, task_cond_gen) %>%
  summarise(mean_rt = mean(response_time_ms, na.rm = TRUE), .groups = "drop")

p_cn_rt <- ggplot(pred_grid_cn_rt, aes(x = age, y = predicted_rt, color = language_group)) +
  geom_point(data = subject_rt_cn, aes(x = age, y = mean_rt), alpha = 0.5, size = 1.5) +
  geom_line() +
  facet_wrap(~ task_cond_gen, labeller = labeller(task_cond_gen = title_map)) +
  scale_color_manual(values = colors, labels = c("bilingual" = "Bilingual", "monolingual" = "Monolingual")) +
  labs(
    title    = "Predicted Response Time by Age and Language Group",
    subtitle = "CN Cohort — Split by Condition",
    x = "Age", y = "Predicted Response Time (ms)", color = "Language Group"
  ) +
  theme_minimal()

print(p_cn_rt)
ggsave(file.path(derivatives_dir, "cn_rt_lmem.png"), plot = p_cn_rt, width = 15, height = 5, dpi = 300)

