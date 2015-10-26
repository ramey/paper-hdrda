library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 8

num_iterations <- 250

K <- 3
sample_sizes <- rep(25, K)
means <- c(1/2, 1)
autocorrelations <- c(0.1, 0.5, 0.9)
feature_dim <- seq(100, 500, by=100)
block_size <- 100
contamination_probs <- seq(0, 0.2, by=0.05)
uncontaminated_var <- rep(1, K)
contaminated_var <- rep(100, K)

# The number of observations to generate from each population as a test data set.
test_sizes <- rep(10000, K)

iter_product <- itertools2::iproduct(p=feature_dim,
                                     contamination_prob=contamination_probs,
                                     mu=means,
                                     rho=autocorrelations,
                                     times=num_iterations)
sim_config <- itertools2::ienum(iter_product)
num_iterations <- length(feature_dim) * length(contamination_probs) * length(means) * length(autocorrelations) * num_iterations

flog.logger("sim", INFO, appender=appender.file('sim-block-contaminated.log'))

results <- mclapply(sim_config, function(sim_i) {
  i <- sim_i$index
  p <- sim_i$value$p
  mu <- sim_i$value$mu
  rho <- rep(sim_i$value$rho, K)
  contamination_prob <- sim_i$value$contamination_prob
  flog.info("Dimension: %s. Contamination Prob: %s -- Sim: %s of %s",
            p, contamination_prob, i, num_iterations, name="sim")
  set.seed(i)

  # Difference in means in the first 100 variables (first block)
  mu1 <- rep(0, p)
  mu2 <- c(rep(mu, block_size), rep(0, p - block_size))
  mu3 <- -mu2
  num_blocks <- p / block_size

  # Generates training data
  train_data <- generate_contaminated(n=sample_sizes,
                                      mu=cbind(mu1, mu2, mu3),
                                      block_size=block_size,
                                      num_blocks=num_blocks,
                                      rho=rho,
                                      uncontaminated_var=uncontaminated_var,
                                      contaminated_var=contaminated_var,
                                      contamination_prob=contamination_prob)
  train_x <- train_data$x
  train_y <- train_data$y

  flog.info("Training Data: (Observations, Dimension): (%s, %s) -- Sim: %s of %s", nrow(train_x), ncol(train_x), i, num_iterations, name="sim")

  # Generates test data
  test_data <- generate_contaminated(n=test_sizes,
                                     mu=cbind(mu1, mu2, mu3),
                                     block_size=block_size,
                                     num_blocks=num_blocks,
                                     rho=autocorrelations,
                                     uncontaminated_var=uncontaminated_var,
                                     contaminated_var=contaminated_var,
                                     contamination_prob=contamination_prob)
  test_x <- test_data$x
  test_y <- test_data$y

  flog.info("Test Data: (Observations, Dimension): (%s, %s) -- Sim: %s of %s", nrow(test_x), ncol(test_x), i, num_iterations, name="sim")

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

  flog.info("kNN Error Rate: %s -- Sim: %s of %s", knn_errors, i, num_iterations, name="sim")

  # Random Forest
  rf_errors <- try({
    rf_out <- randomForest(x=train_x,
                           y=train_y,
                           ntree=250,
                           maxnodes=100)
    rf_predict <- predict(rf_out, test_x)
    mean(rf_predict != test_y)
  })

  flog.info("Random Forest Error Rate: %s -- Sim: %s of %s", rf_errors, i, num_iterations, name="sim")

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

  flog.info("SVM Error Rate: %s -- Sim: %s of %s", ksvm_radial_errors, i, num_iterations, name="sim")

  # Witten and Tibshirani (2011) - JRSS B
  Witten_errors <- try({
      Witten_out <- Witten_Tibshirani(train_x = train_x,
                                      train_y = train_y,
                                      test_x = test_x)
      mean(Witten_out$predictions != test_y)
  })

  flog.info("Witten Error Rate: %s -- Sim: %s of %s", Witten_errors, i, num_iterations, name="sim")

  # Guo, Hastie, and Tibshirani (2007) - Biostatistics
  Guo_errors <- try({
      Guo_out <- Guo(train_x = train_x,
                     train_y = train_y,
                     test_x = test_x,
                     prior = prior_probs)
      mean(Guo_out != test_y)
  })

  flog.info("Guo Error Rate: %s -- Sim: %s of %s", Guo_errors, i, num_iterations, name="sim")

  # HDRDA - Ridge
  hdrda_ridge_errors <- try({
      cv_out <- hdrda_cv(x = train_x,
                         y = train_y,
                         prior = prior_probs)
      hdrda_ridge <- list(lambda = cv_out$lambda, gamma = cv_out$gamma)
      flog.info("HDRDA Ridge. Lambda: %s. Gamma: %s", cv_out$lambda, cv_out$gamma, name="sim")
      mean(predict(cv_out, test_x)$class != test_y)
  })

  flog.info("HDRDA Ridge Error Rate: %s -- Sim: %s of %s", hdrda_ridge_errors, i, num_iterations, name="sim")

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

  flog.info("HDRDA Convex Error Rate: %s -- Sim: %s of %s", hdrda_convex_errors, i, num_iterations, name="sim")

  # Tong, Chen, and Zhao (2012) - Bioinformatics
  Tong_errors <- try({
      Tong_out <- dlda(x = train_x,
                       y = train_y,
                       est_mean = "tong",
                       prior = prior_probs)
      mean(predict(Tong_out, test_x)$class != test_y)
  })

  flog.info("Tong Error Rate: %s -- Sim: %s of %s", Tong_errors, i, num_iterations, name="sim")

  # Pang, Tong, and Zhao (2009) - Bioinformatics
  Pang_errors <- try({
      Pang_out <- sdlda(x = train_x,
                        y = train_y,
                        prior = prior_probs)
      mean(predict(Pang_out, test_x)$class != test_y)
  })

  flog.info("Pang Error Rate: %s -- Sim: %s of %s", Pang_errors, i, num_iterations, name="sim")

  error_rates <- list(
    Guo = Guo_errors,
    HDRDA_Ridge = hdrda_ridge_errors,
    HDRDA_Convex = hdrda_convex_errors,
    kNN = knn_errors,
    Pang = Pang_errors,
    Random_Forest = rf_errors,
    SVM_Radial = ksvm_radial_errors,
    Tong = Tong_errors,
    Witten = Witten_errors
  )

  list(error_rates = error_rates,
       p = p,
       mu=mu,
       rho=rho,
       contamination_prob = contamination_prob,
       hdrda_ridge = hdrda_ridge,
       hdrda_convex = hdrda_convex
  )
}, mc.cores = num_cores)

results_block_contaminated <- list(sim_results=results,
                                   num_iterations=num_iterations)

save(results_block_contaminated,
     file='data/results-sim-block-contaminated.RData')
