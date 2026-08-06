## R code for CVMA-RC

required_packages <- c(
  "mlogit", "dplyr", "quadprog", "MASS", "evd", "lmtest",
  "xgboost", "nnet", "ggplot2", "reshape2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    "Install the following package(s) before running the simulation: ",
    paste(missing_packages, collapse = ", ")
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

## ======================================================================
## Configuration
## ======================================================================
n_set <- c(200, 400, 600)
choices <- 3
random_beta_mean <- rep(1, 8)
features <- length(random_beta_mean)
feature_names <- paste0("Feature", seq_len(features))

scenario_list <- list(
  list(name = "S1", sd = c(2, 1.5, 1, 0.5, 0.5, 0, 0, 0)),
  list(name = "S2", sd = c(2, 1.5, 0.5, 0.5, 0, 0, 0, 0)),
  list(name = "S3", sd = c(2, 1, 0.5, 0.5, 0, 0, 0, 0)),
  list(name = "S4", sd = c(2, 0.5, 0.5, 0.5, 0.5, 0, 0, 0)),
  list(name = "S5", sd = c(2, 0.5, 0.5, 0.5, 0, 0, 0, 0)),
  list(name = "S6", sd = c(2, 0.5, 0.5, 0, 0, 0, 0, 0)),
  list(name = "S7", sd = c(2, 0, 0, 0, 0, 0, 0, 0))
)

sample_indices <- seq_along(n_set)
scenario_indices <- seq_along(scenario_list)

K <- 5
n_rep <- 100
max_attempts_per_scenario <- 10 * n_rep
mixing_draws <- 250
true_probability_draws <- 1000

# rho_x = 0 reproduces the independent-covariate setting in the original code.
# Set rho_x = 0.3 to use the AR(1)-correlated covariates in the Scenario 2 script.
rho_x <- 0

xgb_config <- list(
  nrounds = 150,
  eta = 0.05,
  max_depth = 6,
  min_child_weight = 1,
  subsample = 0.8,
  colsample_bytree = 0.8
)

nn_config <- list(
  size = 10,
  decay = 1e-4,
  maxit = 300
)

candidate_names <- feature_names
rpars_can_list <- lapply(candidate_names, function(feature) {
  setNames("n", feature)
})
rpars <- setNames(rep("n", features), feature_names)

n_candidates <- length(candidate_names)
equal_weights <- rep(1 / n_candidates, n_candidates)

method_names <- c(
  "CVMA-RC", "CV Model", "EW", "TRUE", "ALL", "MNL", "MT Model",
  "XGBoost", "Neural Network"
)

n_methods <- length(method_names)
n_sample_sizes <- length(n_set)
n_scenarios <- length(scenario_list)

formula_object <- as.formula(
  paste("choice ~", paste(feature_names, collapse = " + "), "| 0")
)

## ======================================================================
## Helper functions
## ======================================================================
softmax <- function(x) {
  x <- x - max(x)
  exp(x) / sum(exp(x))
}

validate_probabilities <- function(prob, n_expected, n_choices, eps = 1e-12) {
  if (is.null(dim(prob))) {
    if (length(prob) != n_expected * n_choices) {
      stop("Prediction vector has an unexpected length.")
    }
    prob <- matrix(prob, nrow = n_expected, ncol = n_choices, byrow = TRUE)
  } else {
    prob <- as.matrix(prob)
  }
  
  if (length(dim(prob)) != 2 || any(dim(prob) != c(n_expected, n_choices))) {
    stop(
      "Prediction matrix has dimensions ", paste(dim(prob), collapse = " x "),
      "; expected ", n_expected, " x ", n_choices, "."
    )
  }
  
  for (i in seq_len(n_expected)) {
    missing <- !is.finite(prob[i, ])
    if (sum(missing) == 1 && all(is.finite(prob[i, !missing]))) {
      prob[i, missing] <- 1 - sum(prob[i, !missing])
    }
  }
  
  if (!all(is.finite(prob))) {
    stop("Predicted probabilities contain non-finite values.")
  }
  
  prob <- pmax(prob, eps)
  row_totals <- rowSums(prob)
  if (any(!is.finite(row_totals)) || any(row_totals <= 0)) {
    stop("Predicted probabilities have invalid row sums.")
  }
  
  prob / row_totals
}

gen_X_correlated <- function(n, n_choices, n_features, rho) {
  Sigma_x <- outer(
    seq_len(n_features), seq_len(n_features),
    function(i, j) rho ^ abs(i - j)
  )
  X_matrix <- MASS::mvrnorm(
    n * n_choices,
    mu = rep(0, n_features),
    Sigma = Sigma_x
  )
  
  X_3d <- array(0, dim = c(n, n_choices, n_features))
  for (k in seq_len(n_features)) {
    X_3d[, , k] <- matrix(X_matrix[, k], nrow = n, byrow = TRUE)
  }
  X_3d
}

make_long_data <- function(X_3d, chosen_alternative) {
  n <- dim(X_3d)[1]
  n_choices <- dim(X_3d)[2]
  n_features <- dim(X_3d)[3]
  
  feature_data <- lapply(seq_len(n_features), function(k) {
    as.vector(t(X_3d[, , k]))
  })
  feature_data <- as.data.frame(feature_data)
  names(feature_data) <- feature_names
  
  long_data <- data.frame(
    Person = rep(seq_len(n), each = n_choices),
    Alternative = rep(seq_len(n_choices), times = n),
    feature_data,
    choice = as.integer(
      rep(chosen_alternative, each = n_choices) ==
        rep(seq_len(n_choices), times = n)
    )
  )
  
  mlogit::mlogit.data(
    long_data,
    shape = "long",
    chid.var = "Person",
    choice = "choice",
    alt.var = "Alternative"
  )
}

true_prob <- function(X_3d, beta_mean, beta_sd, n_draws) {
  n <- dim(X_3d)[1]
  Sigma_beta <- diag(beta_sd^2, nrow = length(beta_sd))
  
  beta_draws <- MASS::mvrnorm(
    n = n_draws,
    mu = beta_mean,
    Sigma = Sigma_beta
  )
  beta_draws <- matrix(beta_draws, nrow = n_draws, ncol = length(beta_mean))
  
  probabilities <- matrix(0, nrow = n, ncol = choices)
  for (i in seq_len(n)) {
    utilities <- X_3d[i, , ] %*% t(beta_draws)
    draw_probabilities <- apply(utilities, 2, softmax)
    probabilities[i, ] <- rowMeans(draw_probabilities)
  }
  
  validate_probabilities(probabilities, n, choices)
}

make_ml_matrix <- function(X_3d) {
  n <- dim(X_3d)[1]
  n_choices <- dim(X_3d)[2]
  n_features <- dim(X_3d)[3]
  
  output <- matrix(NA_real_, nrow = n, ncol = n_choices * n_features)
  output_names <- character(n_choices * n_features)
  
  for (j in seq_len(n_choices)) {
    columns <- ((j - 1) * n_features + 1):(j * n_features)
    output[, columns] <- X_3d[, j, ]
    output_names[columns] <- paste0("Alt", j, "_", feature_names)
  }
  
  colnames(output) <- output_names
  output
}

standardize_train_test <- function(train_x, test_x) {
  center <- colMeans(train_x)
  scale <- apply(train_x, 2, stats::sd)
  scale[!is.finite(scale) | scale == 0] <- 1
  
  list(
    train = sweep(sweep(train_x, 2, center, "-"), 2, scale, "/"),
    test = sweep(sweep(test_x, 2, center, "-"), 2, scale, "/")
  )
}

fit_xgboost <- function(train_x, train_y, test_x, config, seed) {
  dtrain <- xgboost::xgb.DMatrix(
    data = train_x,
    label = as.integer(train_y) - 1L
  )
  dtest <- xgboost::xgb.DMatrix(data = test_x)
  
  params <- list(
    objective = "multi:softprob",
    eval_metric = "mlogloss",
    num_class = choices,
    eta = config$eta,
    max_depth = config$max_depth,
    min_child_weight = config$min_child_weight,
    subsample = config$subsample,
    colsample_bytree = config$colsample_bytree,
    seed = seed,
    nthread = 1
  )
  
  model <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = config$nrounds,
    verbose = 0
  )
  
  prediction <- predict(model, dtest)
  validate_probabilities(prediction, nrow(test_x), choices)
}

fit_neural_network <- function(train_x, train_y, test_x, config, seed) {
  scaled <- standardize_train_test(train_x, test_x)
  targets <- diag(choices)[as.integer(train_y), , drop = FALSE]
  
  set.seed(seed)
  model <- nnet::nnet(
    x = scaled$train,
    y = targets,
    size = config$size,
    decay = config$decay,
    maxit = config$maxit,
    softmax = TRUE,
    trace = FALSE,
    MaxNWts = 10000
  )
  
  prediction <- predict(model, scaled$test, type = "raw")
  validate_probabilities(prediction, nrow(test_x), choices)
}

solve_cv_weights <- function(prediction_array, observed_choice, ridge = 1e-8) {
  n_observations <- dim(prediction_array)[1]
  n_choices <- dim(prediction_array)[2]
  n_models <- dim(prediction_array)[3]
  
  prediction_matrix <- matrix(
    aperm(prediction_array, c(2, 1, 3)),
    nrow = n_observations * n_choices,
    ncol = n_models,
    byrow = FALSE
  )
  
  Dmat <- crossprod(prediction_matrix) + diag(ridge, n_models)
  dvec <- crossprod(prediction_matrix, as.vector(t(observed_choice)))
  Amat <- t(rbind(
    matrix(1, nrow = 1, ncol = n_models),
    diag(n_models),
    -diag(n_models)
  ))
  bvec <- c(1, rep(0, n_models), rep(-1, n_models))
  
  solution <- quadprog::solve.QP(
    Dmat = Dmat,
    dvec = as.vector(dvec),
    Amat = Amat,
    bvec = bvec,
    meq = 1
  )$solution
  
  solution[solution < 0] <- 0
  if (!all(is.finite(solution)) || sum(solution) <= 0) {
    return(rep(1 / n_models, n_models))
  }
  solution / sum(solution)
}

average_predictions <- function(prediction_array, weights) {
  output <- apply(
    prediction_array,
    c(1, 2),
    function(probabilities) sum(probabilities * weights)
  )
  validate_probabilities(output, dim(prediction_array)[1], dim(prediction_array)[2])
}

fit_mlogit_model <- function(train_data, random_features = character(0)) {
  model_args <- list(
    formula = formula_object,
    data = train_data
  )
  
  if (length(random_features) > 0) {
    # Pass the evaluated vector through do.call(). This prevents mlogit from
    # retaining the wrapper's local symbol `random_features` in model$call.
    model_args$rpar <- setNames(
      rep("n", length(random_features)),
      as.character(random_features)
    )
    model_args$halton <- NA
    model_args$R <- mixing_draws
    model_args$seed <- 123
  }
  
  do.call(mlogit::mlogit, model_args)
}

predict_mlogit_model <- function(model, test_data, n_test) {
  validate_probabilities(
    predict(model, newdata = test_data),
    n_test,
    choices
  )
}

mt_random_features <- function(mnl_model, train_data, alpha = 0.05) {
  selected <- character(0)
  test_data <- train_data
  test_data$mnl_probability <- as.vector(t(mnl_model$probabilities))
  
  for (feature in feature_names) {
    auxiliary_name <- paste0("z_", feature)
    augmented <- test_data %>%
      dplyr::group_by(Person) %>%
      dplyr::mutate(
        weighted_mean = sum(mnl_probability * .data[[feature]], na.rm = TRUE),
        !!auxiliary_name := 0.5 * (.data[[feature]] - weighted_mean)^2
      ) %>%
      dplyr::ungroup()
    
    test_data[[auxiliary_name]] <- augmented[[auxiliary_name]]
    augmented_formula <- as.formula(
      paste(
        "choice ~",
        paste(c(feature_names, auxiliary_name), collapse = " + "),
        "| 0"
      )
    )
    augmented_model <- mlogit::mlogit(augmented_formula, data = test_data)
    test_result <- lmtest::waldtest(mnl_model, augmented_model)
    p_value <- as.matrix(test_result)[2, "Pr(>Chisq)"]
    
    if (is.finite(p_value) && p_value < alpha) {
      selected <- c(selected, feature)
    }
  }
  
  selected
}

prediction_metrics <- function(true_probability, predicted_probability, observed_class) {
  predicted_probability <- validate_probabilities(
    predicted_probability,
    nrow(true_probability),
    ncol(true_probability)
  )
  
  predicted_class <- max.col(predicted_probability, ties.method = "first")
  chosen_probability <- predicted_probability[
    cbind(seq_len(nrow(predicted_probability)), as.integer(observed_class))
  ]
  
  c(
    MSPE = mean((true_probability - predicted_probability)^2),
    Accuracy = mean(predicted_class == as.integer(observed_class)),
    LogScore = mean(log(chosen_probability))
  )
}

## ======================================================================
## Result containers
## ======================================================================
scenario_names <- vapply(scenario_list, function(x) x$name, character(1))

result_dim <- c(n_sample_sizes, n_scenarios, n_rep, n_methods)
result_dimnames <- list(
  n = as.character(n_set),
  scenario = scenario_names,
  replication = as.character(seq_len(n_rep)),
  method = method_names
)

mspe_by_rep <- array(NA_real_, dim = result_dim, dimnames = result_dimnames)
accuracy_by_rep <- array(NA_real_, dim = result_dim, dimnames = result_dimnames)
log_score_by_rep <- array(NA_real_, dim = result_dim, dimnames = result_dimnames)

weights_by_rep <- array(
  NA_real_,
  dim = c(n_sample_sizes, n_scenarios, n_rep, n_candidates),
  dimnames = list(
    n = as.character(n_set),
    scenario = scenario_names,
    replication = as.character(seq_len(n_rep)),
    candidate = candidate_names
  )
)

cv_selected_by_rep <- array(
  NA,
  dim = c(n_sample_sizes, n_scenarios, n_rep, n_candidates),
  dimnames = dimnames(weights_by_rep)
)
mt_selected_by_rep <- cv_selected_by_rep

successful_replications <- matrix(
  0L,
  nrow = n_sample_sizes,
  ncol = n_scenarios,
  dimnames = list(n = as.character(n_set), scenario = scenario_names)
)

## ======================================================================
## Main simulation
## ======================================================================
start_time <- Sys.time()

for (ni in sample_indices) {
  n <- n_set[ni]
  
  for (sig in scenario_indices) {
    scenario <- scenario_list[[sig]]
    random_beta_sigma <- diag(scenario$sd^2, nrow = features)
    
    scenario_mspe <- matrix(
      NA_real_, nrow = n_rep, ncol = n_methods,
      dimnames = list(NULL, method_names)
    )
    scenario_accuracy <- scenario_mspe
    scenario_log_score <- scenario_mspe
    scenario_weights <- matrix(
      NA_real_, nrow = n_rep, ncol = n_candidates,
      dimnames = list(NULL, candidate_names)
    )
    scenario_cv_selected <- matrix(
      NA, nrow = n_rep, ncol = n_candidates,
      dimnames = list(NULL, candidate_names)
    )
    scenario_mt_selected <- scenario_cv_selected
    
    successful <- 0L
    attempts <- 0L
    
    while (successful < n_rep && attempts < max_attempts_per_scenario) {
      attempts <- attempts + 1L
      simulation_seed <- 777L + 100000L * ni + 1000L * sig + attempts
      
      tryCatch({
        set.seed(simulation_seed)
        
        train_indices <- sort(sample(
          seq_len(n),
          size = floor(0.5 * n),
          replace = FALSE
        ))
        test_indices <- setdiff(seq_len(n), train_indices)
        n_train <- length(train_indices)
        n_test <- length(test_indices)
        
        X_3d <- gen_X_correlated(n, choices, features, rho_x)
        random_beta <- MASS::mvrnorm(
          n,
          mu = random_beta_mean,
          Sigma = random_beta_sigma
        )
        random_beta <- matrix(random_beta, nrow = n, ncol = features)
        
        epsilon <- matrix(
          evd::rgumbel(n * choices),
          nrow = n,
          ncol = choices
        )
        utilities <- matrix(0, nrow = n, ncol = choices)
        for (i in seq_len(n)) {
          utilities[i, ] <- X_3d[i, , ] %*% random_beta[i, ] + epsilon[i, ]
        }
        
        chosen_alternative <- max.col(utilities, ties.method = "first")
        observed_choice <- sapply(
          seq_len(choices),
          function(j) as.integer(chosen_alternative == j)
        )
        
        true_probability <- true_prob(
          X_3d,
          beta_mean = random_beta_mean,
          beta_sd = scenario$sd,
          n_draws = true_probability_draws
        )
        P_test <- true_probability[test_indices, , drop = FALSE]
        
        data_mlogit <- make_long_data(X_3d, chosen_alternative)
        train_data <- data_mlogit[data_mlogit$Person %in% train_indices, ]
        test_data <- data_mlogit[data_mlogit$Person %in% test_indices, ]
        train_choice <- observed_choice[train_indices, , drop = FALSE]
        test_class <- chosen_alternative[test_indices]
        
        ## --------------------------------------------------------------
        ## Eight one-random-coefficient candidates and CVMA-RC
        ## --------------------------------------------------------------
        candidate_test_probability <- array(
          NA_real_, dim = c(n_test, choices, n_candidates)
        )
        for (f in seq_len(n_candidates)) {
          candidate_model <- fit_mlogit_model(
            train_data,
            random_features = candidate_names[f]
          )
          candidate_test_probability[, , f] <- predict_mlogit_model(
            candidate_model,
            test_data,
            n_test
          )
        }
        
        fold_id <- sample(rep(seq_len(K), length.out = n_train))
        folds <- split(seq_len(n_train), fold_id)
        cv_candidate_probability <- array(
          NA_real_, dim = c(n_train, choices, n_candidates)
        )
        
        for (k in seq_len(K)) {
          fold_positions <- folds[[k]]
          fold_persons <- train_indices[fold_positions]
          fold_train_data <- train_data[
            !(train_data$Person %in% fold_persons),
          ]
          fold_validation_data <- train_data[
            train_data$Person %in% fold_persons,
          ]
          
          for (f in seq_len(n_candidates)) {
            fold_model <- fit_mlogit_model(
              fold_train_data,
              random_features = candidate_names[f]
            )
            cv_candidate_probability[fold_positions, , f] <-
              predict_mlogit_model(
                fold_model,
                fold_validation_data,
                length(fold_positions)
              )
          }
        }
        
        if (!all(is.finite(cv_candidate_probability))) {
          stop("Cross-validation predictions contain non-finite values.")
        }
        
        cv_weights <- solve_cv_weights(cv_candidate_probability, train_choice)
        cvma_probability <- average_predictions(
          candidate_test_probability,
          cv_weights
        )
        ew_probability <- average_predictions(
          candidate_test_probability,
          equal_weights
        )
        
        ## --------------------------------------------------------------
        ## Parametric comparison methods
        ## --------------------------------------------------------------
        mnl_model <- fit_mlogit_model(train_data)
        mnl_probability <- predict_mlogit_model(mnl_model, test_data, n_test)
        
        all_model <- fit_mlogit_model(train_data, feature_names)
        all_probability <- predict_mlogit_model(all_model, test_data, n_test)
        
        true_random_features <- feature_names[scenario$sd > 0]
        true_model <- fit_mlogit_model(train_data, true_random_features)
        true_model_probability <- predict_mlogit_model(
          true_model,
          test_data,
          n_test
        )
        
        cv_random_features <- candidate_names[cv_weights > 1 / n_candidates]
        if (setequal(cv_random_features, feature_names)) {
          cv_model_probability <- all_probability
        } else if (length(cv_random_features) == 0) {
          cv_model_probability <- mnl_probability
        } else {
          cv_model <- fit_mlogit_model(train_data, cv_random_features)
          cv_model_probability <- predict_mlogit_model(
            cv_model,
            test_data,
            n_test
          )
        }
        
        mt_features <- mt_random_features(mnl_model, train_data)
        if (setequal(mt_features, feature_names)) {
          mt_probability <- all_probability
        } else if (length(mt_features) == 0) {
          mt_probability <- mnl_probability
        } else {
          mt_model <- fit_mlogit_model(train_data, mt_features)
          mt_probability <- predict_mlogit_model(mt_model, test_data, n_test)
        }
        
        ## --------------------------------------------------------------
        ## Machine-learning comparison methods
        ## One row per person; alternative-feature blocks are columns.
        ## --------------------------------------------------------------
        ml_matrix <- make_ml_matrix(X_3d)
        ml_train_x <- ml_matrix[train_indices, , drop = FALSE]
        ml_test_x <- ml_matrix[test_indices, , drop = FALSE]
        ml_train_y <- chosen_alternative[train_indices]
        
        xgb_probability <- fit_xgboost(
          ml_train_x,
          ml_train_y,
          ml_test_x,
          config = xgb_config,
          seed = simulation_seed + 1L
        )
        nn_probability <- fit_neural_network(
          ml_train_x,
          ml_train_y,
          ml_test_x,
          config = nn_config,
          seed = simulation_seed + 2L
        )
        
        ## --------------------------------------------------------------
        ## Common out-of-sample evaluation
        ## --------------------------------------------------------------
        method_probability <- array(
          NA_real_,
          dim = c(n_test, choices, n_methods),
          dimnames = list(NULL, NULL, method_names)
        )
        method_probability[, , "CVMA-RC"] <- cvma_probability
        method_probability[, , "CV Model"] <- cv_model_probability
        method_probability[, , "EW"] <- ew_probability
        method_probability[, , "TRUE"] <- true_model_probability
        method_probability[, , "ALL"] <- all_probability
        method_probability[, , "MNL"] <- mnl_probability
        method_probability[, , "MT Model"] <- mt_probability
        method_probability[, , "XGBoost"] <- xgb_probability
        method_probability[, , "Neural Network"] <- nn_probability
        
        replication_number <- successful + 1L
        for (method in method_names) {
          metrics <- prediction_metrics(
            P_test,
            method_probability[, , method],
            test_class
          )
          scenario_mspe[replication_number, method] <- metrics["MSPE"]
          scenario_accuracy[replication_number, method] <- metrics["Accuracy"]
          scenario_log_score[replication_number, method] <- metrics["LogScore"]
        }
        
        scenario_weights[replication_number, ] <- cv_weights
        scenario_cv_selected[replication_number, ] <-
          candidate_names %in% cv_random_features
        scenario_mt_selected[replication_number, ] <-
          candidate_names %in% mt_features
        
        successful <- replication_number
        
        if (successful == 1L || successful %% 10L == 0L) {
          current_mspe <- colMeans(
            scenario_mspe[seq_len(successful), , drop = FALSE],
            na.rm = TRUE
          )
          cat(
            "n=", n,
            "; scenario=", scenario$name,
            "; successful=", successful,
            "; attempts=", attempts,
            "; minimum MSPE=", method_names[which.min(current_mspe)],
            "\n",
            sep = ""
          )
        }
      }, error = function(e) {
        cat(
          "Attempt ", attempts,
          " failed for n=", n,
          ", scenario=", scenario$name,
          ": ", conditionMessage(e), "\n",
          sep = ""
        )
      })
    }
    
    if (successful < n_rep) {
      warning(
        "Scenario ", scenario$name, " with n=", n,
        " stopped after ", attempts, " attempts and ", successful,
        " successful replications."
      )
    }
    
    successful_replications[ni, sig] <- successful
    mspe_by_rep[ni, sig, , ] <- scenario_mspe
    weights_by_rep[ni, sig, , ] <- scenario_weights
    cv_selected_by_rep[ni, sig, , ] <- scenario_cv_selected
    mt_selected_by_rep[ni, sig, , ] <- scenario_mt_selected
  }
}

mean_mspe <- apply(mspe_by_rep, c(1, 2, 4), mean, na.rm = TRUE)
mean_weights <- apply(weights_by_rep, c(1, 2, 4), mean, na.rm = TRUE)
cv_selection_rate <- apply(cv_selected_by_rep, c(1, 2, 4), mean, na.rm = TRUE)
mt_selection_rate <- apply(mt_selected_by_rep, c(1, 2, 4), mean, na.rm = TRUE)

end_time <- Sys.time()
print(end_time - start_time)

## ======================================================================
## Save results
## ======================================================================
# save results
setwd('D:\\A工作\\2024_2\\CVMA\\simulation')
save(  mspe_by_rep,
       accuracy_by_rep,
       log_score_by_rep,
       mean_mspe,
       mean_accuracy,
       mean_log_score,
       weights_by_rep,
       mean_weights,
       cv_selection_rate,
       mt_selection_rate,
       successful_replications,
       scenario_list,
       n_set,
       method_names, file = "Scenario1_ml_v1.RData")


dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

save(
  mspe_by_rep,
  accuracy_by_rep,
  log_score_by_rep,
  mean_mspe,
  mean_accuracy,
  mean_log_score,
  weights_by_rep,
  mean_weights,
  cv_selection_rate,
  mt_selection_rate,
  successful_replications,
  scenario_list,
  n_set,
  method_names,
  file = file.path(output_dir, "Scenario1_random_coefficients_with_ML.RData")
)

## ======================================================================
## Plots
## ======================================================================
method_colors <- setNames(
  c(
    "#00A6A6", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3",
    "#A6D854", "#FFD92F", "#1F78B4", "#6A3D9A"
  ),
  method_names
)

mspe_data <- reshape2::melt(
  mspe_by_rep,
  varnames = c("SampleSize", "Scenario", "Replication", "Method"),
  value.name = "MSPE"
) %>%
  dplyr::filter(is.finite(MSPE)) %>%
  dplyr::mutate(
    SampleSize = factor(
      SampleSize,
      levels = as.character(n_set),
      labels = paste0("n = ", n_set)
    ),
    Scenario = factor(Scenario, levels = scenario_names),
    Method = factor(Method, levels = method_names)
  )

mspe_plot <- ggplot2::ggplot(
  mspe_data,
  ggplot2::aes(x = Method, y = MSPE, color = Method, fill = Method)
) +
  ggplot2::geom_boxplot(
    width = 0.8,
    alpha = 0.25,
    linewidth = 0.45,
    outlier.size = 0.25
  ) +
  ggplot2::facet_grid(SampleSize ~ Scenario, scales = "free_y") +
  ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
  ggplot2::scale_fill_manual(values = method_colors, drop = FALSE) +
  ggplot2::labs(x = NULL, y = "MSPE") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    legend.position = "none",
    axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 6),
    axis.text.y = ggplot2::element_text(size = 7),
    strip.text = ggplot2::element_text(size = 8)
  )

print(mspe_plot)
ggplot2::ggsave(
  file.path(figure_dir, "scenario1_mspe.pdf"),
  mspe_plot,
  width = 42,
  height = 25,
  units = "cm"
)

weight_data <- reshape2::melt(
  mean_weights,
  varnames = c("SampleSize", "Scenario", "Candidate"),
  value.name = "Weight"
) %>%
  dplyr::filter(is.finite(Weight)) %>%
  dplyr::mutate(
    SampleSize = factor(
      SampleSize,
      levels = as.character(n_set),
      labels = paste0("n = ", n_set)
    ),
    Scenario = factor(Scenario, levels = scenario_names),
    Candidate = factor(Candidate, levels = candidate_names),
    BarColor = ifelse(Weight > 1 / n_candidates, "#00A6A6", "#E78AC3")
  )

weight_plot <- ggplot2::ggplot(
  weight_data,
  ggplot2::aes(x = Candidate, y = Weight, fill = BarColor)
) +
  ggplot2::geom_col(width = 0.8) +
  ggplot2::facet_grid(SampleSize ~ Scenario) +
  ggplot2::scale_fill_identity() +
  ggplot2::scale_y_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    expand = ggplot2::expansion(mult = c(0, 0.03))
  ) +
  ggplot2::coord_flip() +
  ggplot2::labs(x = "Random-coefficient candidate", y = "CVMA-RC weight") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text = ggplot2::element_text(size = 7),
    strip.text = ggplot2::element_text(size = 8)
  )

print(weight_plot)
ggplot2::ggsave(
  file.path(figure_dir, "scenario1_weights.pdf"),
  weight_plot,
  width = 38,
  height = 25,
  units = "cm"
)

