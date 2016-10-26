library(ProjectTemplate)
library(sortinghat)
library(microbenchmark)

load.project()

# Number of timing-comparison replications to perform
num_reps <- 5

# Grid Size (m x m)
grid_size <- 5

# Simulation Parameters
K <- 4
n <- rep.int(x=25, times=K)
#feature_dims <- c(
#  25, 50, 75,
#  seq(100, 1000, by=100),
# seq(1500, 3000, by=500)
#)
feature_dims <- c(25, 50)

# klaR's Implementation of Friedman (1989)
# Uses a grid model selection
rda_grid <- function(x, y, num_tuning=5) {
  tuning_grid <- expand.grid(lambda=seq(0, 1, length=num_tuning),
                             gamma=seq(0, 1, length=num_tuning))
   rda_errors <- vector(mode="list", length=nrow(tuning_grid))
   for (i in seq_len(nrow(tuning_grid))) {
     rda_errors[[i]] <- try({
       rda_out <- klaR::rda(x=x,
                            grouping=y,
                            lambda=tuning_grid$lambda[i],
                            gamma=tuning_grid$gamma[i],
                            estimate.error=FALSE)
       rda_out$error.rate['crossval']
     })
   }

   rda_errors[!sapply(rda_errors, is.numeric)] <- NA
   tuning_grid[which.min(unlist(rda_errors)), ]
}

# klaR's Implementation of Friedman (1989)
# Uses a Nelder-Mead model selection
rda_nelder_mead <- function(x, y) {
  rda_out <- klaR::rda(x=x,
                       grouping=y,
                       estimate.error=FALSE)
  rda_out$regularization
}

# HDRDA
hdrda_grid <- function(x, y, num_tuning=5) {
  hdrda_cv(x=x,
           y=y,
           num_lambda=num_tuning,
           num_gamma=num_tuning,
           shrinkage_type="convex")
}

# Timing comparison for various values of p
timing_results <- lapply(feature_dims, function(p) {
  # Creates a list of 4 mean vectors. Each vector has a length of `p`.
  means <- lapply(c(-3, -1, 1, 3), rep.int, times=p)
  cov_mat <- diag(p)

  # Generates multivariate normal data
  set.seed(42)
  data <- simdata_normal(n=n, mean=means, cov=cov_mat)

  mb_results <- microbenchmark(
    rda_grid(x=data$x, y=data$y, num_tuning=grid_size),
    rda_nelder_mead(x=data$x, y=data$y),
    hdrda_grid(x=data$x, y=data$y, num_tuning=grid_size),
    times=num_reps)
  cat("Timing results for p =", p, "\n")
  print(mb_results)

  mb_results
})
names(timing_results) <- feature_dims

save(timing_results, file="data/timing-results.RData")
