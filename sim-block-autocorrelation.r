library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 2

num_iterations <- 1000

K <- 3
sample_sizes <- rep(25, K)
autocorrelations <- c(0.1, 0.5, 0.9)
feature_dim <- seq(50, 500, by=50)
block_size <- 25

# The number of observations to generate from each population as a test data set.
test_sizes <- rep(1000, K)

sim_config <- rep(feature_dim, each=num_iterations)

flog.logger("sim", INFO, appender=appender.file('sim-block-autocorrelation.log'))

results <- mclapply(seq_along(sim_config), function(i) {
  p <- sim_config[i]
  flog.info("Dimension: %s -- Sim:  %s of %s", p, i, length(sim_config), name="sim")
  set.seed(i)

  mu1 <- rep(0, p)
  mu2 <- c(rep(1/2, p), rep(0, p - 50))
  mu3 <- c(rep(-1/2, p), rep(0, p - 50))
  num_blocks <- p / block_size

  # Generates training data
  train_data <- generate_blockdiag(n=sample_sizes,
                                   mu=cbind(mu1, mu2, mu3),
                                   block_size=block_size,
                                   num_blocks=num_blocks,
                                   rho=autocorrelations)

  train_x <- train_data$x
  train_y <- train_data$y

  # Generates test data
  test_data <- generate_blockdiag(n=test_sizes,
                                  mu=cbind(mu1, mu2, mu3),
                                  block_size=block_size,
                                  num_blocks=num_blocks,
                                  rho=autocorrelations)

  test_x <- test_data$x
  test_y <- test_data$y

  # Removes the loaded data from memory and then garbage collects.
  rm(train_data)
  rm(test_data)
  gc(reset = TRUE)

  # Sets the prior probabilities as equal
  num_classes <- nlevels(train_y)
  prior_probs <- rep(1, num_classes) / num_classes

  # kNN
  knn_errors <- try({
      tune.knn_out <- tune.knn(x = train_x,
                               y = train_y,
                               k = 1:5,
                               tunecontrol = tune.control(sampling = "cross", cross = 10))
      knn_out <- knn(train = train_x,
                     cl = train_y,
                     test = test_x,
                     k = tune.knn_out$best.model$k)
      mean(knn_out != test_y)
  })

  # SVM with Radial Basis Functions
  ksvm_radial_errors <- try({
      ksvm_out <- ksvm(x = train_x,
                       y = train_y,
                       kernel = "rbfdot",
                       kpar = "automatic",
                       cross = 10)
      ksvm_radial <- predict(ksvm_out, test_x)
      mean(ksvm_radial != test_y)
  })

  # Witten and Tibshirani (2011) - JRSS B
  Witten_errors <- try({
      Witten_out <- Witten_Tibshirani(train_x = train_x,
                                      train_y = train_y,
                                      test_x = test_x)
      mean(Witten_out$predictions != test_y)
  })

  # Guo, Hastie, and Tibshirani (2007) - Biostatistics
  Guo_errors <- try({
      Guo_out <- Guo(train_x = train_x,
                     train_y = train_y,
                     test_x = test_x,
                     prior = prior_probs)
      mean(Guo_out != test_y)
  })

  # HDRDA - Ridge
  hdrda_ridge_errors <- try({
      cv_out <- hdrda_cv(x = train_x,
                         y = train_y,
                         prior = prior_probs)
      hdrda_ridge <- list(lambda = cv_out$lambda, gamma = cv_out$gamma)
      flog.info("HDRDA Ridge. Lambda: %s. Gamma: %s", cv_out$lambda, cv_out$gamma, name="sim")
      mean(predict(cv_out, test_x)$class != test_y)
  })

  # HDRDA - Convex
  hdrda_convex_errors <- try({
      cv_out <- hdrda_cv(x = train_x,
                         y = train_y,
                         prior = prior_probs,
                         num_gamma = 21,
                         shrinkage_type = "convex")
      hdrda_convex <- list(lambda = cv_out$lambda, gamma = cv_out$gamma)
      flog.info("HDRDA Convex. Lambda: %s. Gamma: %s", cv_out$lambda, cv_out$gamma, name="sim")
      mean(predict(cv_out, test_x)$class != test_y)
  })

  # Tong, Chen, and Zhao (2012) - Bioinformatics
  Tong_errors <- try({
      Tong_out <- dlda(x = train_x,
                       y = train_y,
                       est_mean = "tong",
                       prior = prior_probs)
      mean(predict(Tong_out, test_x)$class != test_y)
  })

  # Pang, Tong, and Zhao (2009) - Bioinformatics
  Pang_errors <- try({
      Pang_out <- sdlda(x = train_x,
                        y = train_y,
                        prior = prior_probs)
      mean(predict(Pang_out, test_x)$class != test_y)
  })

  error_rates <- list(
    Guo = Guo_errors,
    HDRDA_Ridge = hdrda_ridge_errors,
    HDRDA_Convex = hdrda_convex_errors,
    kNN = knn_errors,
    Pang = Pang_errors,
    SVM_Radial = ksvm_radial_errors,
    Tong = Tong_errors,
    Witten = Witten_errors
  )

  list(error_rates = error_rates,
       p = p,
       hdrda_ridge = hdrda_ridge,
       hdrda_convex = hdrda_convex
  )
}, mc.cores = num_cores)

results_block_autocorrelation <- list(sim_results=results,
                                      num_iterations=num_iterations)

save(results_block_autocorrelation,
     file='data/results-sim-block-autocorrelation.RData')
