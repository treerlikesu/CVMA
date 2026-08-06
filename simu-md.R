library(mlogit)
library(dplyr)
library(tidyr)
library(Rsolnp)
library(quadprog)
library(truncnorm)
library(randtoolbox)
library(triangle)
library(nloptr)
library(MASS)
library(evd)
library(abind)
library(mixtools)

par(mfrow = c(2, 3), mar = c(3, 2, 1, 1))

n_set <- c(800)
nn <- length(n_set)
prefix <- "Feature"
choices <- 3
random_beta_mean <- rep(1, 8)
features <- length(random_beta_mean)
feature_to_use <- paste0(prefix, 1:features)
num_ran <- 2
feature_random <- feature_to_use[1:num_ran]

type_list <- c("n", "t", "u", "zbt", "zbu", "ln")
rpars_can_list <- lapply(type_list, function(type) {
  setNames(rep(type, length(feature_random)), feature_random)
})
rpars <- do.call(c, rpars_can_list)

nme <- length(type_list)
ncm_list <- c("CVMA-MD", type_list)
ncm <- length(ncm_list)

K <- 5
rep <- 100
max_attempts_per_scenario <- rep * 10

# Scenario 2 settings in the manuscript plus a positive log-normal case.
# For S7, mean and sd are specified on the original coefficient scale.
scenario_list <- list(
  list(name = "S1_normal_low", type = "normal", mean = c(1, 1), sd = c(0.5, 0.5)),
  list(name = "S2_normal_medium", type = "normal", mean = c(1, 1), sd = c(1, 1)),
  list(name = "S3_normal_high", type = "normal", mean = c(1, 1), sd = c(2.5, 2.5)),
  list(name = "S4_one_sided", type = "one_sided"),
  list(name = "S5_mixed_shapes", type = "mixed_shapes"),
  list(name = "S6_correlated", type = "correlated", mean = c(1, 1), sd = 2.5, rho = 0.5),
  list(name = "S7_lognormal", type = "lognormal", mean = c(1, 1), sd = c(0.5, 0.5))
)
ns <- length(scenario_list)

softmax <- function(x) {
  x <- x - max(x)
  exp(x) / sum(exp(x))
}

correct_na <- function(matrix) {
  if (!any(is.na(matrix))) {
    return(matrix)
  } else {
    for (i in 1:nrow(matrix)) {
      row_na_count <- sum(is.na(matrix[i, ]))
      if (row_na_count > 1) {
        return(matrix)
      } else if (row_na_count == 1) {
        na_col_index <- which(is.na(matrix[i, ]))
        matrix[i, na_col_index] <- 1 - sum(matrix[i, !is.na(matrix[i, ])], na.rm = TRUE)
      }
    }
    return(matrix)
  }
}

# Random coefficient generator for Scenario 2.
sample_random_beta <- function(n, scenario) {
  out <- matrix(0, nrow = n, ncol = num_ran)
  
  if (scenario$type == "normal") {
    out[, 1] <- rnorm(n, mean = scenario$mean[1], sd = scenario$sd[1])
    out[, 2] <- rnorm(n, mean = scenario$mean[2], sd = scenario$sd[2])
  }
  
  if (scenario$type == "mixed_shapes") {
    raw_tri <- rtriangle(n, a = 0, b = 8, c = 1)
    out[, 1] <- 1 + 2.5 * raw_tri
    out[, 2] <- rnorm(n, mean = 1, sd = 2.5)
  }
  
  if (scenario$type == "one_sided") {
    out[,1] <- runif(n, 0, 1)
    out[,2] <- runif(n, 0, 1)
  }
  
  if (scenario$type == "correlated") {
    Sigma <- matrix(c(
      scenario$sd^2, scenario$rho * scenario$sd^2,
      scenario$rho * scenario$sd^2, scenario$sd^2
    ), nrow = 2, byrow = TRUE)
    out <- MASS::mvrnorm(n = n, mu = scenario$mean, Sigma = Sigma)
  }
  
  if (scenario$type == "lognormal") {
    # Convert the desired original-scale moments to rlnorm parameters.
    sdlog <- sqrt(log1p((scenario$sd / scenario$mean)^2))
    meanlog <- log(scenario$mean) - 0.5 * sdlog^2
    for (k in seq_len(num_ran)) {
      out[, k] <- rlnorm(n, meanlog = meanlog[k], sdlog = sdlog[k])
    }
  }
  
  return(out)
}

true_prob <- function(X_3d, scenario, n_draws = 1000) {
  n <- dim(X_3d)[1]
  P <- matrix(0, nrow = n, ncol = choices)
  fixed_beta <- rep(1, features - num_ran)
  
  for (i in 1:n) {
    beta_random <- sample_random_beta(n_draws, scenario)
    beta_draws <- cbind(beta_random, matrix(fixed_beta, nrow = n_draws, ncol = features - num_ran, byrow = TRUE))
    utilities <- X_3d[i, , ] %*% t(beta_draws)
    draw_probs <- apply(utilities, 2, softmax)
    P[i, ] <- rowMeans(draw_probs)
  }
  P
}

gen_X_correlated <- function(n, choices, features, rho = 0.3) {
  Sigma_x <- outer(1:features, 1:features, function(i, j) rho ^ abs(i - j))
  X_mat <- MASS::mvrnorm(n * choices, mu = rep(0, features), Sigma = Sigma_x)
  X_3d <- array(0, dim = c(n, choices, features))
  idx <- 1
  for (i in 1:n) {
    for (j in 1:choices) {
      X_3d[i, j, ] <- X_mat[idx, ]
      idx <- idx + 1
    }
  }
  X_3d
}

set.seed(777)
formula_string <- paste(feature_to_use, collapse = " + ")
formula_object <- as.formula(paste("choice", "~", formula_string, "|0"))

mspe_by_rep <- array(
  NA_real_,
  dim = c(nn, ns, rep, ncm),
  dimnames = list(n = n_set, scenario = vapply(scenario_list, `[[`, "", "name"),
                  replication = seq_len(rep), method = ncm_list)
)
mean_mspe <- array(
  NA_real_,
  dim = c(nn, ns, ncm),
  dimnames = list(n = n_set, scenario = vapply(scenario_list, `[[`, "", "name"),
                  method = ncm_list)
)
weights_by_rep <- array(
  NA_real_,
  dim = c(nn, ns, rep, nme),
  dimnames = list(n = n_set, scenario = vapply(scenario_list, `[[`, "", "name"),
                  replication = seq_len(rep), distribution = type_list)
)
mean_weights <- array(
  NA_real_,
  dim = c(nn, ns, nme),
  dimnames = list(n = n_set, scenario = vapply(scenario_list, `[[`, "", "name"),
                  distribution = type_list)
)
a <- Sys.time()

for (ni in 1:nn) {
  n <- n_set[ni]
  P <- array(0, dim = c(n, choices))
  U <- array(0, dim = c(n, choices))
  X_3d <- array(0, dim = c(n, choices, features))
  successful_simulations <- 1
  g <- 0
  
  for (sig in seq_len(ns)) {
    scenario_mspe <- matrix(NA_real_, nrow = rep, ncol = ncm,
                            dimnames = list(NULL, ncm_list))
    scenario_weights <- matrix(NA_real_, nrow = rep, ncol = nme,
                               dimnames = list(NULL, type_list))
    successful_simulations <- 1
    g <- 0
    
    while (successful_simulations < rep + 1 && g < max_attempts_per_scenario) {
      print(paste(ni, sig, successful_simulations, scenario_list[[sig]]$name))
      result <- tryCatch({
        set.seed(g + 777)
        g <- g + 1
        train_indices <- sort(sample(seq_len(n), size = 0.5 * n, replace = FALSE))
        test_indices <- setdiff(1:n, train_indices)
        n_train <- length(train_indices)
        n_test <- length(test_indices)
        
        random_beta <- matrix(rep(random_beta_mean, each = n), nrow = n, byrow = FALSE)
        X_3d <- gen_X_correlated(n, choices, features, rho = 0.3)
        
        rand_beta_2 <- sample_random_beta(n, scenario_list[[sig]])
        random_beta[, 1:num_ran] <- rand_beta_2
        
        epsilon <- matrix(rgumbel(n * choices), nrow = n, ncol = choices)
        for (i in 1:n) {
          U[i, ] <- X_3d[i, , ] %*% random_beta[i, ] + epsilon[i, ]
        }
        
        y <- matrix(0, n, choices)
        y_sta <- apply(U, 1, which.max)
        y <- sapply(1:choices, function(x) as.integer(y_sta == x))
        
        X_long <- data.frame()
        for (i in 1:n) {
          person_df <- data.frame(
            Person = rep(i, choices),
            Alternative = 1:choices,
            Feature1 = X_3d[i, , 1],
            Feature2 = X_3d[i, , 2],
            Feature3 = X_3d[i, , 3],
            Feature4 = X_3d[i, , 4],
            Feature5 = X_3d[i, , 5],
            Feature6 = X_3d[i, , 6],
            Feature7 = X_3d[i, , 7],
            Feature8 = X_3d[i, , 8],
            choice = as.integer(y_sta[i] == 1:choices)
          )
          X_long <- rbind(X_long, person_df)
        }
        
        data_mlogit <- mlogit.data(
          X_long,
          shape = "long",
          chid.var = "Person",
          choice = "choice",
          alt.var = "Alternative"
        )
        
        P <- true_prob(X_3d, scenario_list[[sig]], n_draws = 1000)
        if (!all(is.finite(P))) {
          print(paste("Skipping iteration due to non-finite P at iteration", g))
          next
        }
        
        P_test <- P[test_indices, ]
        train_data <- data_mlogit[data_mlogit$Person %in% train_indices, ]
        test_data <- data_mlogit[data_mlogit$Person %in% test_indices, ]
        
        train_mat <- y[train_indices, ]
        
        candidates_pred <- array(0, dim = c(n_test, choices, nme))
        for (f in 1:nme) {
          model <- mlogit(
            formula_object,
            train_data,
            rpar = rpars_can_list[[f]],
            halton = NA,
            R = 500,
            seed = 123
          )
          candidates_pred[, , f] <- correct_na(predict(model, newdata = test_data))
        }
        
        if (!all(is.finite(candidates_pred))) {
          print(paste("Skipping iteration due to non-finite candidates_pred at iteration", g))
          next
        }
        
        one_fold_size <- n_train %/% K
        fold_sizes <- rep(one_fold_size, K)
        fold_sizes[seq_len(n_train %% K)] <- fold_sizes[seq_len(n_train %% K)] + 1
        folds <- split(1:n_train, rep(seq_len(K), fold_sizes))
        tilde_p <- array(0, dim = c(n_train, choices, nme))
        
        for (k in 1:K) {
          fold_persons <- train_indices[folds[[k]]]
          fold_train_data <- train_data[!(train_data$Person %in% fold_persons), ]
          fold_test_data <- train_data[train_data$Person %in% fold_persons, ]
          for (f in 1:nme) {
            train_model <- mlogit(
              formula_object,
              fold_train_data,
              rpar = rpars_can_list[[f]],
              halton = NA,
              R = 500,
              seed = 123
            )
            tilde_p[folds[[k]], , f] <- correct_na(predict(train_model, newdata = fold_test_data))
          }
        }
        
        if (!all(is.finite(tilde_p))) {
          print(paste("Skipping iteration due to non-finite tilde_p at iteration", g))
          next
        }
        
        tilde_p_2d <- matrix(aperm(tilde_p, c(2, 1, 3)), nrow = n_train * choices, ncol = nme, byrow = FALSE)
        Dmat <- t(tilde_p_2d) %*% tilde_p_2d
        dvec <- t(tilde_p_2d) %*% as.vector(t(train_mat))
        Amat <- t(rbind(matrix(1, nrow = 1, ncol = nme), diag(nme), -diag(nme)))
        bvec <- rbind(1, matrix(0, nrow = nme, ncol = 1), matrix(-1, nrow = nme, ncol = 1))
        result <- solve.QP(Dmat = Dmat, dvec = dvec, Amat = Amat, bvec = bvec, meq = 1)
        
        w <- as.matrix(result$solution)
        w <- w * (w > 0)
        wm <- w / sum(w)
        
        cv_pred <- apply(candidates_pred, c(1, 2), function(x, wm) sum(x * wm), wm)
        
        scenario_weights[successful_simulations, ] <- wm
        
        ncm_array <- array(0, dim = c(n_test, choices, ncm))
        ncm_array[, , 1] <- cv_pred
        
        for (i in 2:ncm) {
          ncm_array[, , i] <- candidates_pred[, , (i - 1)]
        }
        
        for (i in 1:ncm) {
          scenario_mspe[successful_simulations, i] <- mean((P_test - ncm_array[, , i])^2)
        }
        
        if (successful_simulations > 1) {
          print(paste0("successful simulations = ", successful_simulations))
          current_mean_mspe <- colMeans(scenario_mspe[1:successful_simulations, , drop = FALSE])
          print(current_mean_mspe)
          print(paste0("minimum MSPE: ", ncm_list[which.min(current_mean_mspe)]))
          print(colMeans(scenario_weights[1:successful_simulations, , drop = FALSE]))
        }
        
        successful_simulations <- successful_simulations + 1
      }, error = function(e) {
        cat("Error encountered:", conditionMessage(e), "\n")
      })
    }
    
    if (successful_simulations < rep + 1) {
      warning(
        "Scenario ", scenario_list[[sig]]$name,
        " stopped after ", g, " attempts with ",
        successful_simulations - 1, " successful replications."
      )
    }
    
    mspe_by_rep[ni, sig, , ] <- scenario_mspe
    mean_mspe[ni, sig, ] <- colMeans(scenario_mspe)
    weights_by_rep[ni, sig, , ] <- scenario_weights
    mean_weights[ni, sig, ] <- colMeans(scenario_weights)
    
    print(paste0("scenario=", scenario_list[[sig]]$name, "; n=", n_set[ni],
                 "; minimum MSPE: ", ncm_list[which.min(mean_mspe[ni, sig, ])]))
  }
}

# save results
# setwd('D:\\A工作\\2024_2\\CVMA\\simulation')
# save(mspe_by_rep, mean_mspe, weights_by_rep, mean_weights, file = "Scenario2_results_v2.RData")

##################################################################################
## plot

library(ggplot2)
library(dplyr)
library(reshape2)

# Run this script after the simulation has created these two objects.
required_objects <- c("mspe_by_rep", "mean_weights")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), envir = .GlobalEnv)
]
if (length(missing_objects) > 0) {
  stop(
    "Run the simulation first. Missing object(s): ",
    paste(missing_objects, collapse = ", ")
  )
}

scenario_labels <- c(
  S1_normal_low = "S1: Normal (SD = 0.5)",
  S2_normal_medium = "S2: Normal (SD = 1)",
  S3_normal_high = "S3: Normal (SD = 2.5)",
  S4_one_sided = "S4: One-sided",
  S5_mixed_shapes = "S5: Mixed-shaped",
  S6_correlated = "S6: Correlated normal",
  S7_lognormal = "S7: Log-normal (SD = 0.5)"
)

method_levels <- dimnames(mspe_by_rep)$method
distribution_levels <- dimnames(mean_weights)$distribution

method_colors <- setNames(rep("#E995C9", length(method_levels)), method_levels)
method_colors["CVMA-MD"] <- "#00B8B2"


# -------------------------------------------------------------------------
# Figure 1: MSPE boxplots over simulation replications
# -------------------------------------------------------------------------
mspe_data <- reshape2::melt(
  mspe_by_rep,
  varnames = c("SampleSize", "ScenarioCode", "Replication", "Method"),
  value.name = "MSPE"
) %>%
  filter(!is.na(MSPE)) %>%
  mutate(
    Scenario = factor(
      as.character(ScenarioCode),
      levels = names(scenario_labels),
      labels = unname(scenario_labels)
    ),
    Method = factor(as.character(Method), levels = method_levels),
    Sample = paste0("Sample Size = ", SampleSize)
  )

mspe_plot <- ggplot(
  mspe_data,
  aes(x = Method, y = MSPE, color = Method, fill = Method)
) +
  geom_boxplot(
    width = 0.8,
    alpha = 0.30,
    linewidth = 0.7,
    outlier.size = 0.45
  ) +
  facet_wrap(~Scenario, scales = "free_y", ncol = 4) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  scale_fill_manual(values = method_colors, drop = FALSE) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 5),
    expand = expansion(mult = c(0.09, 0.09))
  ) +
  labs(x = NULL, y = "MSPE") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.spacing = grid::unit(0.8, "cm"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

print(mspe_plot)

ggsave(
  filename = file.path(figure_dir, "dmpsesimu_S1_S7.pdf"),
  plot = mspe_plot,
  dpi = 300,
  width = 28,
  height = 16,
  units = "cm"
)

# -------------------------------------------------------------------------
# Figure 2: mean CVMA-MD weights
# -------------------------------------------------------------------------
weight_data <- reshape2::melt(
  mean_weights,
  varnames = c("SampleSize", "ScenarioCode", "Distribution"),
  value.name = "Weight"
) %>%
  filter(!is.na(Weight)) %>%
  mutate(
    Scenario = factor(
      as.character(ScenarioCode),
      levels = names(scenario_labels),
      labels = unname(scenario_labels)
    ),
    Distribution = factor(
      as.character(Distribution),
      levels = distribution_levels
    ),
    Sample = paste0("Sample Size = ", SampleSize)
  ) %>%
  group_by(SampleSize, Scenario) %>%
  mutate(
    EqualWeight = 1 / n(),
    BarColor = if_else(Weight > EqualWeight, "#00B8B2", "#E995C9")
  ) %>%
  ungroup()

weight_plot <- ggplot(
  weight_data,
  aes(x = Distribution, y = Weight, fill = BarColor)
) +
  geom_col(width = 0.8) +
  facet_wrap(~Scenario, scales = "free_x", ncol = 4) +
  scale_fill_identity() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = "Distribution", y = "Weight") +
  coord_flip() +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    panel.spacing = grid::unit(0.8, "cm"),
    strip.text = element_text(size = 10),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9)
  )

print(weight_plot)

ggsave(
  filename = file.path(figure_dir, "wdissimu_S1_S7.pdf"),
  plot = weight_plot,
  dpi = 600,
  width = 28,
  height = 16,
  units = "cm"
)

