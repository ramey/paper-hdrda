library(ProjectTemplate)
library(sortinghat)
library(microbenchmark)

load.project()

# Number of timing-comparison replications to perform
# num_reps <- 50
num_reps <- 5

# Grid Size (m x m)
grid_size <- 5

# Simulation Parameters
K <- 4
n <- rep.int(x=25, times=K)
# p <- c(10, 100, 1000, 10000)
p <- c(10, 20)

# TODO: Update means
means <- lapply(c(-3, -1, 1, 3), rep.int, times=p)
cov_mat <- diag(p)

# Generates multivariate normal data
set.seed(42)
data <- simdata_normal(n=n, mean=means, cov=cov_mat)

# klaR's Implementation of Friedman (1989)
rda_model_selection <- function(x, y, num_tuning=5) {
  tuning_grid <- expand.grid(lambda=seq(0, 1, length=num_tuning),
                             gamma=seq(0, 1, length=num_tuning))

  rda_errors <- vector(mode="list", length=nrow(tuning_grid))

  for (i in seq_len(nrow(tuning_grid))) {
    rda_errors[[i]] <- try({
      rda_out <- rda(x=x,
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

# HDRDA
hdrda_model_selection <- function(x, y, num_tuning=5) {
  hdrda_cv(x=x,
           y=y,
           num_lambda=num_tuning,
           num_gamma=num_tuning,
           shrinkage_type="convex")
}

timing_results <- microbenchmark(
  rda_model_selection(x=data$x, y=data$y, num_tuning=grid_size),
  hdrda_model_selection(x=data$x, y=data$y, num_tuning=grid_size),
  times=num_reps)

save(timing_results, file="data/timing-results.RData") 

