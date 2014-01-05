library(sparsediscrim)
library(sortinghat)
library(microbenchmark)
library(klaR)

# Simulation Parameters
K <- 4
n <- rep.int(x = 25, times = K)
p <- 250
means <- lapply(c(-3, -1, 1, 3), rep.int, times = p)
cov_mat <- diag(p)

# Generates multivariate normal data
set.seed(42)
data <- simdata_normal(n = n, mean = means, cov = cov_mat)

# klaR's Implementation of Friedman (1989)
rda_model_selection <- function(num_tuning = 3) {
  tuning_grid <- expand.grid(lambda = seq(0, 1, length = num_tuning),
                             gamma = seq(0, 1, length = num_tuning))

  rda_errors <- vector(mode = "list", length = nrow(tuning_grid))

  for (i in seq_len(nrow(tuning_grid))) {
    rda_errors[[i]] <- try({
      rda_out <- rda(x = data$x,
                     grouping = data$y,
                     estimate.error = FALSE)
      rda_out$error.rate['crossval']
    })
  }

  rda_errors[!sapply(rda_errors, is.numeric)] <- NA

  tuning_grid[which.min(unlist(rda_errors)), ]
}  

# HDRDA
hdrda_model_selection <- function(num_tuning = 3) {
  hdrda_cv(x = data$x,
           y = data$y,
           num_lambda = num_tuning,
           num_gamma = num_tuning,
           shrinkage_type = "convex")
}

timing_results <- list()

# 3 x 3 Grid of Tuning Parameters
timing_results[["3x3"]] <- microbenchmark(
  rda_model_selection(num_tuning = 3),
  hdrda_model_selection(num_tuning = 3),
  times = 50)

print(timing_results[["3x3"]])

# 4 x 4 Grid of Tuning Parameters
timing_results[["4x4"]] <- microbenchmark(
  rda_model_selection(num_tuning = 4),
  hdrda_model_selection(num_tuning = 4),
  times = 50)

print(timing_results[["4x4"]])

# 5 x 5 Grid of Tuning Parameters
timing_results[["5x5"]] <- microbenchmark(
  rda_model_selection(num_tuning = 5),
  hdrda_model_selection(num_tuning = 5),
  times = 50)

print(timing_results[["5x5"]])

save(timing_results, file = "data/timing-results.RData") 
