# Helpers for the archived joint model (checks/joint_model/03_joint_model_{AD,MCI}.Rmd); not part of the current analysis.
#
# Used by 03_joint_model_{AD,MCI}.Rmd (JMbayes2 joint model).
#
# Inputs : results/{EP}/inputs/            (written by 01_prepare_inputs.ipynb)
# Outputs: results/{EP}/{model}/predictions.csv, plots/{EP}/FigureS*.png|pdf
#          (smoke mode: results/smoke/... and plots/smoke/...)
#          one row per participant x fit, columns ID, cohort, fit, set, model, history, S_24 ... S_60
#          S_m = P(event-free at landmark + m months | event-free at the landmark).
# Performance metrics are NOT computed here: 04_evaluate_and_plot.ipynb scores every model with the
# same IBS / AUC / Brier / calibration / bootstrap code.

suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})

project_dir <- function() {
  Sys.getenv("SURVIVAL_PROJECT_DIR", "/Users/yihan/Documents/Github/survival_Anlysis_Clustering")
}

read_ids <- function(path) read.csv(path, colClasses = c(ID = "character"), check.names = FALSE)

results_root <- function(mode) {
  if (mode == "standard") file.path(project_dir(), "results") else file.path(project_dir(), "results", "smoke")
}

plots_root <- function(mode) {
  if (mode == "standard") file.path(project_dir(), "plots") else file.path(project_dir(), "plots", "smoke")
}

load_shared_inputs <- function(endpoint, mode) {
  dir <- file.path(results_root(mode), endpoint, "inputs")
  if (!file.exists(file.path(dir, "config.json"))) {
    stop("Shared inputs not found in ", dir, ": run 01_prepare_inputs.ipynb first.")
  }
  cfg <- fromJSON(file.path(dir, "config.json"))
  sh <- list(
    dir = dir, config = cfg,
    cohort = read_ids(file.path(dir, "landmark_cohort.csv")),
    pre48 = read_ids(file.path(dir, "visits_pre48.csv")),
    dev = read_ids(file.path(dir, "development_history.csv")),
    folds = read_ids(file.path(dir, "outer_folds.csv"))
  )
  stopifnot(cfg$endpoint == endpoint, cfg$run_mode == mode)
  sh
}

# Outer folds 1..K on NACC plus the locked full-NACC fit evaluated in AIBL.
make_fits <- function(sh) {
  nacc <- sh$cohort$ID[sh$cohort$cohort == "NACC"]
  aibl <- sh$cohort$ID[sh$cohort$cohort == "AIBL"]
  fits <- list()
  for (k in sort(unique(sh$folds$outer_test_fold))) {
    test <- sh$folds$ID[sh$folds$outer_test_fold == k]
    fits[[paste0("outer_", k)]] <- list(train = setdiff(nacc, test), test = test, test_cohort = "NACC")
  }
  fits$full_NACC <- list(train = nacc, test = aibl, test_cohort = "AIBL")
  fits
}

# One row per participant with baseline features and latest pre-landmark values (the
# "latest value frozen" prediction input).
latest_rows <- function(sh, ids, cohort) {
  d <- sh$cohort[sh$cohort$cohort == cohort, ]
  d <- d[match(ids, d$ID), ]
  stopifnot(!anyNA(d$ID))
  d
}

# ---------------------------------------------------------------- preprocessing
# Fitted on training rows only: median imputation, then standardisation.
fit_preprocess <- function(train, cols) {
  med <- vapply(cols, function(c) median(train[[c]], na.rm = TRUE), numeric(1))
  med[is.na(med)] <- 0
  x <- train[, cols, drop = FALSE]
  for (c in cols) x[[c]][is.na(x[[c]])] <- med[[c]]
  center <- vapply(cols, function(c) mean(x[[c]]), numeric(1))
  scale <- vapply(cols, function(c) sd(x[[c]]), numeric(1))
  scale[is.na(scale) | scale == 0] <- 1
  list(cols = cols, median = med, center = center, scale = scale)
}

apply_preprocess <- function(prep, data) {
  for (c in prep$cols) {
    v <- data[[c]]
    v[is.na(v)] <- prep$median[[c]]
    data[[c]] <- (v - prep$center[[c]]) / prep$scale[[c]]
  }
  data
}

# Convert an n x length(grid) survival matrix to the shared prediction format.
prediction_frame <- function(ids, cohort, fit, set, model, history, surv, grid) {
  surv <- as.matrix(surv)
  stopifnot(nrow(surv) == length(ids), ncol(surv) == length(grid))
  colnames(surv) <- paste0("S_", as.integer(round(grid)))
  data.frame(ID = ids, cohort = cohort, fit = fit, set = set, model = model, history = history,
             surv, check.names = FALSE, stringsAsFactors = FALSE)
}

check_survival <- function(surv) {
  surv <- as.matrix(surv)
  if (any(!is.finite(surv)) || any(surv < -1e-8) || any(surv > 1 + 1e-8)) stop("Invalid survival probabilities.")
  if (ncol(surv) > 1 && any(t(apply(surv, 1, diff)) > 1e-6)) stop("Predicted survival must be non-increasing.")
  invisible(TRUE)
}

output_dirs <- function(endpoint, mode, model) {
  res <- file.path(results_root(mode), endpoint, model)
  plt <- file.path(plots_root(mode), endpoint)
  dir.create(res, recursive = TRUE, showWarnings = FALSE)
  dir.create(plt, recursive = TRUE, showWarnings = FALSE)
  list(results = res, plots = plt)
}

write_predictions <- function(pred, sh, endpoint, mode, model, history, settings = list()) {
  dirs <- output_dirs(endpoint, mode, model)
  # Every fit x history needs all test participants; training rows are either complete or absent
  # (the joint model writes training predictions only for the primary history).
  fits <- make_fits(sh)
  for (h in unique(pred$history)) for (name in names(fits)) {
    n_test <- sum(pred$history == h & pred$fit == name & pred$set == "test")
    n_train <- sum(pred$history == h & pred$fit == name & pred$set == "train")
    if (n_test != length(fits[[name]]$test) || !n_train %in% c(0, length(fits[[name]]$train)))
      stop("Prediction rows do not cover every participant: ", h, " / ", name)
  }
  check_survival(pred[, sh$config$survival_columns])
  write.csv(pred, file.path(dirs$results, "predictions.csv"), row.names = FALSE)
  pkgs <- c("survival", "JMbayes2", "nlme")
  info <- list(
    model = model, history = history, endpoint = endpoint, run_mode = mode,
    created = format(Sys.time(), tz = "UTC", usetz = TRUE),
    shared_config = sh$config, settings = settings, r_version = R.version.string,
    packages = lapply(setNames(pkgs, pkgs), function(p)
      if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA))
  writeLines(toJSON(info, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(dirs$results, "run_info.json"))
  invisible(dirs)
}

# ---------------------------------------------------------------- JMbayes2 joint model
# History conditions, identical to landmark_common.HISTORY_LENGTHS: visits at or before the landmark,
# selected backward from it (first visit = the earliest visit only).
HISTORIES <- c("baseline_visit", "1_visit", "2_visits", "3_visits", "all_visits")

select_history <- function(rows, history) {
  key <- if ("cohort" %in% names(rows)) paste(rows$cohort, rows$ID) else rows$ID  # IDs can repeat across cohorts
  o <- order(key, rows$visit_month)
  rows <- rows[o, ]; key <- key[o]
  pos <- ave(rows$visit_month, key, FUN = seq_along)        # 1 = earliest visit
  n <- ave(rows$visit_month, key, FUN = length)
  keep <- switch(history,
                 baseline_visit = pos == 1,
                 "1_visit" = pos > n - 1, "2_visits" = pos > n - 2, "3_visits" = pos > n - 3,
                 all_visits = rep(TRUE, nrow(rows)),
                 stop("Unknown history: ", history))
  rows[keep, ]
}

# MCMC settings. Standard = the analysis; smoke = execution check only.
jm_settings <- function(mode, penalty_alphas, seed) {  # seed: the project seed from config.json
  n_cores <- parallel::detectCores()
  base <- if (mode == "smoke") {
    list(n_chains = 2, n_iter = 600, n_burnin = 200, n_thin = 1, predict_samples = 30)
  } else {
    list(n_chains = 3, n_iter = 3500, n_burnin = 500, n_thin = 1, predict_samples = 200)
  }
  c(base, list(penalty_alphas = penalty_alphas, cores = max(1, min(3, n_cores - 1)),
               predict_cores = max(1, n_cores - 1), predict_parallel = "multicore", seed = seed,
               landmark_time_eps = 1e-4, predict_chunk = 250))
}

# Landmark time scale: t = (visit month - landmark) / 12 years (negative before the landmark);
# survival time = years after the landmark. Each marker: linear mixed model with random
# intercept and slope. Survival submodel: the four baseline features. Association: current
# value of every marker, with a shrinkage prior (default ridge) on the association parameters.
jm_long_data <- function(rows, month_col, prep, markers, landmark) {
  rows <- apply_long_scaling(rows, prep, markers)
  rows$t <- (rows[[month_col]] - landmark) / 12
  rows
}

# Markers are standardised with training means/SDs; missing marker values stay missing
# (each mixed model uses the measurements actually observed).
apply_long_scaling <- function(rows, prep, markers) {
  for (m in markers) rows[[m]] <- (rows[[m]] - prep$center[[m]]) / prep$scale[[m]]
  for (b in prep$basic) {
    v <- rows[[b]]
    v[is.na(v)] <- prep$median[[b]]
    rows[[b]] <- (v - prep$center[[b]]) / prep$scale[[b]]
  }
  rows
}

jm_fit <- function(sh, train_ids, markers, basic, settings) {
  suppressPackageStartupMessages({ library(JMbayes2); library(nlme) })
  L <- sh$config$landmark_month
  dev <- sh$dev[sh$dev$ID %in% train_ids, ]
  prep <- fit_preprocess(dev, c(markers, basic))
  prep$basic <- basic
  long <- jm_long_data(dev, "start", prep, markers, L)
  long <- long[order(long$ID, long$t), ]
  surv <- latest_rows(sh, train_ids, "NACC")
  surv <- apply_long_scaling(surv[, c("ID", "time", "event", basic)], prep, character(0))
  surv$Time <- surv$time / 12
  surv <- surv[order(surv$ID), ]
  # JMbayes2 requires every mixed model to use the same data frame; missing marker values
  # are excluded per model through na.exclude.
  long <- long[, c("ID", "t", markers)]
  lmes <- lapply(markers, function(m) {
    lme(as.formula(paste(m, "~ t")), random = ~ t | ID, data = long, na.action = na.exclude,
        control = lmeControl(opt = "optim", msMaxIter = 200))
  })
  cox <- coxph(as.formula(paste("Surv(Time, event) ~", paste(basic, collapse = " + "))), data = surv)
  priors <- list(penalty_alphas = settings$penalty_alphas)
  fit <- jm(cox, lmes, time_var = "t", priors = priors,
            n_chains = settings$n_chains, n_iter = settings$n_iter, n_burnin = settings$n_burnin,
            n_thin = settings$n_thin, cores = settings$cores, seed = settings$seed)
  list(model = fit, prep = prep, markers = markers, basic = basic)
}

# Conditional survival from the landmark using each participant's own pre-landmark visits of the
# chosen history condition (the model itself is always fitted on the full development history).
jm_predict <- function(fit, sh, ids, cohort, grid, settings, history = "all_visits") {
  L <- sh$config$landmark_month
  rows <- sh$pre48[sh$pre48$cohort == cohort & sh$pre48$ID %in% ids, ]
  rows <- select_history(rows, history)
  nd <- jm_long_data(rows, "visit_month", fit$prep, fit$markers, L)
  nd$Time <- settings$landmark_time_eps  # known event-free at the landmark
  nd$event <- 0
  nd <- nd[order(nd$ID, nd$t), ]
  times <- grid / 12
  chunks <- split(unique(nd$ID), ceiling(seq_along(unique(nd$ID)) / settings$predict_chunk))
  out <- lapply(chunks, function(chunk) {
    p <- predict(fit$model, newdata = nd[nd$ID %in% chunk, ], process = "event", times = times,
                 return_newdata = TRUE,
                 control = list(n_samples = settings$predict_samples,
                                 # "snow" silently runs on 1 core inside JMbayes2; "multicore" forks.
                                 parallel = if (is.null(settings$predict_parallel)) "snow" else settings$predict_parallel,
                                 cores = if (is.null(settings$predict_cores)) settings$cores else settings$predict_cores,
                                 seed = settings$seed))  # project seed, not JMbayes2's default (123)
    p <- p[p$ID %in% chunk, ]
    p$grid_index <- match(round(p$t * 12, 6), round(grid, 6))
    p <- p[!is.na(p$grid_index), c("ID", "grid_index", "pred_CIF")]
    w <- reshape(p, idvar = "ID", timevar = "grid_index", direction = "wide")
    w <- w[, c("ID", paste0("pred_CIF.", seq_along(grid)))]
    w
  })
  w <- do.call(rbind, out)
  w <- w[match(ids, w$ID), ]
  if (anyNA(w$ID)) stop("Joint-model predictions missing for some participants.")
  surv <- 1 - as.matrix(w[, -1])
  # Enforce exact monotonicity against Monte Carlo noise (differences are ~1e-4 or less).
  surv <- t(apply(pmin(pmax(surv, 0), 1), 1, cummin))
  surv
}

jm_association_table <- function(fit, fit_name) {
  s <- summary(fit$model)
  surv <- as.data.frame(s$Survival)
  surv$term <- rownames(s$Survival)
  surv$fit <- fit_name
  rownames(surv) <- NULL
  surv
}
