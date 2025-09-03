####### Loading required packages #######
library(pseudo)
library(geepack)
library(MASS)
library(survival)
library(adjustedCurves)
library(dplyr)
library(ggplot2)
library(ggh4x)

####### Simulation of data #######
generate_data <- function(n, confounding_type) {
  mu <- rep(0, 3)
  sigma <- matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3)
  Z <- mvrnorm(n, mu, sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  Z3 <- Z[, 3]
  
  logit_p1 <- (Z1 + Z2 + Z3) / 3
  p1 <- exp(logit_p1) / (1 + exp(logit_p1))
  D <- rbinom(n, 1, p1)
  
  if (confounding_type == "weak") {
    hazard <- exp(-3 + 1*D + 0.5 * Z1 - 0.5 * Z2)
  } else if (confounding_type == "strong") {
    hazard <- exp(-3 + 1*D + 0.5 * Z1 + 0.5 * Z2)
  }
  T_event <- rexp(n, rate = hazard)
  
  T_censor <- rexp(n, rate = exp(-3.5))
  time <- pmin(T_event, T_censor)
  status <- as.numeric(T_event <= T_censor)
  
  return(data.frame(time, status, D, Z1, Z2, Z3))
}


####### True survival estimates #######
get_true_survival <- function(t, d_level, confounding_type) {
  n_large <- 1000000
  mu <- rep(0, 3)
  sigma <- matrix(c(1, 0.2, 0, 0.2, 1, 0, 0, 0, 1), nrow = 3)
  Z <- mvrnorm(n_large, mu, sigma)
  Z1 <- Z[, 1]
  Z2 <- Z[, 2]
  Z3 <- Z[, 3]
  D <- rep(d_level, n_large)
  
  if (confounding_type == "weak") {
    hazard <- exp(-3 + 1*D + 0.5 * Z1 - 0.5 * Z2)
  } else if (confounding_type == "strong") {
    hazard <- exp(-3 + 1*D + 0.5 * Z1 + 0.5 * Z2)
  }
  
  T_event_potential <- rexp(n_large, rate = hazard)
  km_fit <- survfit(Surv(T_event_potential, rep(1, n_large)) ~ 1)
  summary_km <- summary(km_fit, times = t)
  
  return(summary_km$surv)
}

####### Un-Adjusted KM #######
unadjusted_km <- function(data, t) {
  km_fit <- survfit(Surv(time, status) ~ 1, data = data[data$D == 1, ])
  return(summary(km_fit, times = t)$surv)
}

####### IPTW-Pseudo (Manually coded) #######
iptw_pseudo_estimator_manual <- function(data, t, ps_model_correct) {
  if (ps_model_correct) {
    ps_model <- glm(D ~ Z1 + Z2 + Z3, data = data, family = "binomial")
  } else {
    ps_model <- glm(D ~ Z1 + Z3, data = data, family = "binomial")
  }
  p1_hat <- predict(ps_model, type = "response")
  
  data_d1 <- data[data$D == 1, ]
  p1_hat_d1 <- p1_hat[data$D == 1]
  
  pseudo_obs_list <- pseudo::pseudosurv(time = data_d1$time, 
                                        event = data_d1$status, 
                                        tmax = t)
  
  iptw_estimate <- sapply(1:length(t), function(i) {
    weighted.mean(pseudo_obs_list$pseudo[, i], w = 1 / p1_hat_d1)
  })
  
  return(iptw_estimate)
}

####### IPTW-Pseudo (R Package) #######
iptw_pseudo_estimator_adjustedsurv <- function(data, t, ps_model_correct){
  
  if (ps_model_correct) {
    ps_model <- glm(D ~ Z1 + Z2 + Z3, data = data, family = "binomial")
  } else {
    ps_model <- glm(D ~ Z1 + Z3, data = data, family = "binomial")
  }
  
  data$D <- as.factor(data$D)
  
  package_adj_surv <- adjustedsurv(data = data,
                                   variable = "D",
                                   ev_time = "time",
                                   event = "status",
                                   method = "iptw_pseudo",
                                   treatment_model = ps_model,
                                   times = t,
                                   conf_int = FALSE)
  
  package_results <- package_adj_surv$adj[package_adj_surv$adj$group == 1, "surv"]
  
  return(package_results)
}

####### DR-Pseudo (Manually coded) ####### 
aiptw_pseudo_estimator_manual <- function(data, t, ps_model_correct, or_model_correct) {
  n <- nrow(data)
  
  if (ps_model_correct) {
    ps_model <- glm(D ~ Z1 + Z2 + Z3, data = data, family = "binomial")
  } else {
    ps_model <- glm(D ~ Z1 + Z3, data = data, family = "binomial")
  }
  p1_hat <- predict(ps_model, type = "response")
  
  pseudo_obs <- pseudosurv(time = data$time, event = data$status, tmax = t)$pseudo
  
  data_d1 <- data[data$D == 1, ]
  
  if (or_model_correct) {
    or_model_cox <- coxph(Surv(time, status) ~ Z1 + Z2, data = data_d1)
  } else {
    or_model_cox <- coxph(Surv(time, status) ~ Z1 + Z3, data = data_d1)
  }

  base_surv <- survfit(or_model_cox, newdata = data.frame(Z1=0, Z2=0, Z3=0)) # Centered covariates
  
  s0_t <- summary(base_surv, times = t)$surv
  
  if (or_model_correct) {
    lp <- or_model_cox$coefficients["Z1"] * data$Z1 + or_model_cox$coefficients["Z2"] * data$Z2
  } else {
    lp <- or_model_cox$coefficients["Z1"] * data$Z1 + or_model_cox$coefficients["Z3"] * data$Z3
  }
  
  m1_hat_matrix <- sapply(s0_t, function(s) s^exp(lp))
  
  dr_estimate <- sapply(1:length(t), function(i) {
    term1_num <- data$D * pseudo_obs[, i]
    term2_num <- (data$D - p1_hat) * m1_hat_matrix[, i]
    full_term <- (term1_num - term2_num) / p1_hat
    
    mean(full_term)
  })
  
  dr_se_estimate <- sapply(1:length(t), function(i) {
    term1_num <- data$D * pseudo_obs[, i]
    term2_num <- (data$D - p1_hat) * m1_hat_matrix[, i]
    full_term <- (term1_num - term2_num) / p1_hat
    
    sqrt((1/(length(full_term))^2)*sum((full_term - mean(full_term))^2))
  })
  
  return(list(point = dr_estimate,
              se = dr_se_estimate))
}

####### DR-Pseudo (R Package) ####### 
aiptw_pseudo_estimator_adjustedsurv <- function(data, t, ps_model_correct, or_model_correct){
  
  if (ps_model_correct) {
    ps_model <- glm(D ~ Z1 + Z2 + Z3, data = data, family = "binomial")
  } else {
    ps_model <- glm(D ~ Z1 + Z3, data = data, family = "binomial")
  }
  
  if (or_model_correct) {
    or_covars <- c("Z1", "Z2")
  } else {
    or_covars <- c("Z1", "Z3")
  }
  
  data$D <- as.factor(data$D)
  
  package_adj_surv <- adjustedsurv(data = data,
                                   variable = "D",
                                   ev_time = "time",
                                   event = "status",
                                   method = "aiptw_pseudo",
                                   treatment_model = ps_model,
                                   outcome_vars = or_covars,
                                   times = t,
                                   conf_int = TRUE)
  
  package_results <- package_adj_surv$adj[package_adj_surv$adj$group == 1, "surv"]
  se_package_results <- package_adj_surv$adj[package_adj_surv$adj$group == 1, "se"]
  
  return(list(point = package_results,
              se = se_package_results))
}



####### Bootstrap SE estimator #######
estimate_se_bootstrap <- function(data, t, estimator_func, B, seeds, ...) {
  n <- nrow(data)
  boot_results <- matrix(NA, nrow = B, ncol = length(t))
  
  for (j in 1:B) {
    set.seed(seeds[j])
    boot_indices <- sample(1:n, size = n, replace = TRUE)
    boot_data <- data[boot_indices, ]
    boot_results[j, ] <- estimator_func(data = boot_data, t = t, ...)$point
    
    if (j %% 250 == 0) cat(paste("boot iteration:", j, "/", B, "\n"))
  }
  
  se <- apply(boot_results, 2, sd, na.rm = TRUE)
  
  return(se)
}

####### Calculating metrics (relative and non-relative) #######
calculate_metrics_relative <- function(estimates, true_values) {
  diff_mat <- sweep(estimates, 2, true_values, FUN = "-")
  div_mat <- sweep(diff_mat, 2, true_values, FUN = "/")
  relative_bias <- apply(div_mat, 2, mean)
  relative_sd <- apply(estimates, 2, sd) / true_values
  
  relative_mse <- relative_bias^2 + relative_sd^2
  
  return(list(bias = relative_bias,
              sd = relative_sd,
              mse = relative_mse))
}

calculate_metrics_non_relative <- function(estimates, true_values) {
  diff_mat <- sweep(estimates, 2, true_values, FUN = "-")
  bias <- apply(diff_mat, 2, mean)
  sd <- apply(estimates, 2, sd)
  
  mse <- bias^2 + sd^2
  
  return(list(bias = bias,
              sd = sd,
              mse = mse))
}

####### Testing point estimate comparisons #######
test_simulator_point <- function(num_sim,
                                 num_subject,
                                 time_points,
                                 confounding_type,
                                 ps_model_correct,
                                 or_model_correct){
  
  unadjusted_km_results <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  iptw_pseudo_manual_results <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  iptw_pseudo_package_results <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  aiptw_pseudo_manual_results <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  aiptw_pseudo_package_results <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  
  for (i in 1:num_sim){
    set.seed(i)
    dt <- generate_data(n = num_subject, confounding_type = confounding_type)
    unadjusted_km_results[i, ] <- unadjusted_km(data = dt, t = time_points)
    iptw_pseudo_manual_results[i, ] <- iptw_pseudo_estimator_manual(data = dt, t = time_points, 
                                                                    ps_model_correct = ps_model_correct)
    iptw_pseudo_package_results[i, ] <- iptw_pseudo_estimator_adjustedsurv(data = dt, t = time_points, 
                                                                           ps_model_correct = ps_model_correct)
    aiptw_pseudo_manual_results[i, ] <- aiptw_pseudo_estimator_manual(data = dt, t = time_points, 
                                                                      ps_model_correct = ps_model_correct, 
                                                                      or_model_correct = or_model_correct)$point
    aiptw_pseudo_package_results[i, ] <- aiptw_pseudo_estimator_adjustedsurv(data = dt, t = time_points, 
                                                                             ps_model_correct = ps_model_correct, 
                                                                             or_model_correct = or_model_correct)$point
    
    if (i %% 1 == 0) cat(paste("sim iteration:", i, "/", num_sim, "\n"))
    
    }
  
  return(list(unadjusted_km = unadjusted_km_results,
              iptw_pseudo_manual = iptw_pseudo_manual_results,
              iptw_pseudo_package = iptw_pseudo_package_results,
              aiptw_pseudo_manual = aiptw_pseudo_manual_results,
              aiptw_pseudo_package = aiptw_pseudo_package_results))
}


####### Testing SE comparisons #######
test_simulator_se <- function(num_sim,
                              num_subject,
                              time_points,
                              confounding_type,
                              ps_model_correct,
                              or_model_correct){
  
  aiptw_pseudo_manual_model <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  aiptw_pseudo_package_model <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  aiptw_pseudo_manual_bootstrap <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  aiptw_pseudo_package_bootstrap <- data.frame(t1 = numeric(0), t2 = numeric(0), t3 = numeric(0))
  
  for (i in 1:num_sim){
    set.seed(i)
    dt <- generate_data(n = num_subject, confounding_type = confounding_type)
    
    aiptw_pseudo_manual_model[i, ] <- aiptw_pseudo_estimator_manual(data = dt, t = time_points, 
                                                                      ps_model_correct = ps_model_correct, 
                                                                      or_model_correct = or_model_correct)$se
    aiptw_pseudo_package_model[i, ] <- aiptw_pseudo_estimator_adjustedsurv(data = dt, t = time_points, 
                                                                             ps_model_correct = ps_model_correct, 
                                                                             or_model_correct = or_model_correct)$se
    
    num_boot <- 500
  
    set.seed(1235757)
    seeds <- sample.int(.Machine$integer.max, size = num_boot, replace = FALSE)
      
    aiptw_pseudo_manual_bootstrap[i, ] <- estimate_se_bootstrap(data = dt,
                                                                t = time_points,
                                                                estimator_func = aiptw_pseudo_estimator_manual,
                                                                B = num_boot,
                                                                seeds = seeds,
                                                                ps_model_correct = ps_model_correct, 
                                                                or_model_correct = or_model_correct)
    
    aiptw_pseudo_package_bootstrap[i, ] <- estimate_se_bootstrap(data = dt,
                                                                 t = time_points,
                                                                 estimator_func = aiptw_pseudo_estimator_adjustedsurv,
                                                                 B = num_boot,
                                                                 seeds = seeds,
                                                                 ps_model_correct = ps_model_correct, 
                                                                 or_model_correct = or_model_correct)
    
    if (i %% 1 == 0) cat(paste("sim iteration:", i, "/", num_sim, "\n"))
    
  }
  
  return(list(aiptw_pseudo_manual_model = aiptw_pseudo_manual_model,
              aiptw_pseudo_package_model = aiptw_pseudo_package_model,
              aiptw_pseudo_manual_bootstrap = aiptw_pseudo_manual_bootstrap,
              aiptw_pseudo_package_bootstrap = aiptw_pseudo_package_bootstrap))
}


