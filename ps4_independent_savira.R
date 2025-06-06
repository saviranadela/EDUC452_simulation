library(dplyr)
library(ggplot2)

df <- read.csv("gilbert_meta_34.csv")

df <- df %>%
  group_by(id, cluster_id, block_id, wave, treat, cov_female, cov_marital, cov_attend_school, cov_childcare, cov_permission_work) %>%
  summarise(phq9 = sum(resp, na.rm = TRUE), .groups = "drop")

n_individuals <- n_distinct(df$id)
n_waves <- n_distinct(df$wave)

summary_stats <- df %>%
  group_by(wave) %>%
  summarise(
    mean_phq = mean(phq9, na.rm = TRUE),
    sd_phq = sd(phq9, na.rm = TRUE),
    n_obs = n()
  )

table(df$treat) / n_individuals


#DGM

n <- 1000
n_waves <- 5
waves <- 0:3
prop_treat <- 0.6

beta_0 <- 14
beta_time <- -0.5
beta_treat <- -1.0
sigma_id <- 3  
sigma_noise <- 2        

set.seed(123)
ids <- 1:n
treat_vec <- rbinom(n, 1, prob = prop_treat)

sim_df <- expand.grid(id = ids, wave = waves) %>%
  mutate(
    treat = rep(treat_vec, each = length(waves)),
    alpha_id = rep(rnorm(n, 0, sigma_id), each = length(waves)),
    phq_total = beta_0 + beta_time * wave + beta_treat * wave * treat + alpha_id + rnorm(n * length(waves), 0, sigma_noise)
  )

sim_mcar <- sim_df %>%
  mutate(
    missing_mcar = rbinom(n(), 1, prob = 0.3),  # 30% missing randomly
    phq_obs = ifelse(missing_mcar == 1, NA, phq_total)
  )

library(data.table)
setDT(sim_df)
sim_mar <- copy(sim_df)
sim_mar[, lag_phq := shift(phq_total, 1, type = "lag"), by = id]

sim_mar[, missing_mar := rbinom(.N, 1, prob = plogis(0.3 * (lag_phq - 10)))]
sim_mar[, phq_obs := ifelse(missing_mar == 1 & wave != 0, NA, phq_total)]

sim_mnar <- sim_df %>%
  mutate(
    missing_mnar = rbinom(n(), 1, prob = plogis(0.3 * (phq_total - 10))),
    phq_obs = ifelse(missing_mnar == 1 & wave != 0, NA, phq_total)  # keep wave 0 observed
  )

#DAM

library(lme4)

model_mcar <- lmer(phq_obs ~ wave * treat + (1 | id), data = sim_mcar, REML = FALSE)
summary(model_mcar)

model_mar <- lmer(phq_obs ~ wave * treat + (1 | id), data = sim_mar, REML = FALSE)
summary(model_mar)

model_mnar <- lmer(phq_obs ~ wave * treat + (1 | id), data = sim_mnar, REML = FALSE)
summary(model_mnar)

model_full <- lmer(phq_total ~ wave * treat + (1 | id), data = sim_df, REML = FALSE)
summary(model_full)
