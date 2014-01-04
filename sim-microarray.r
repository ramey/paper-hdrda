library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 3

train_pct <- 2/3
num_iterations <- 5
#data_sets <- as.character(describe_data()$author)
#data_sets <- c("burczynski", "nakayama", "shipp", "singh")
data_sets <- c("sorlie", "christensen")

sim_config <- sim_factors(data_sets = data_sets, num_reps = num_iterations,
                          stringsAsFactors = FALSE)

results <- mclapply(seq_len(nrow(sim_config)), function(i) {
  message("Dataset: ", data_set, " -- Sim Config: ", i, " of ", nrow(sim_config))
  set.seed(i)
  
  # Load data and randomly partition into training/test data sets.
  data_set <- sim_config$data_set[i]
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)
  rand_split <- rand_partition(y = data$y, num_partitions = 1,
                               train_pct = train_pct)[[1]]
  train_x <- as.matrix(data$x[rand_split$training, ])
  train_y <- data$y[rand_split$training]
  test_x <- as.matrix(data$x[rand_split$test, ])
  test_y <- data$y[rand_split$test]

  # Removes the loaded data from memory and then garbage collects.
  rm(data)
  gc(reset = TRUE)

  # Sets the prior probabilities as equal
  num_classes <- nlevels(train_y)
  prior_probs <- rep(1, num_classes) / num_classes  

  # Clemmensen, Hastie, Witten and Ersbøll (2011) - Technometrics
  Clemmensen_errors <- try({
      Clemmensen_out <- Clemmensen(train_x = train_x,
                                   train_y = train_y,
                                   test_x = test_x,
                                   normalize_data = TRUE,
                                   lambda = 1e-3,
                                   verbose = TRUE)
      mean(Clemmensen_out != test_y)
  })

  # Friedman (1989)
  rda_errors <- try({
      rda_out <- klaR:::rda(x = train_x,
                            grouping = train_y,
                            prior = prior_probs,
                            output = TRUE)
      rda_class <- klaR:::predict.rda(rda_out, test_x)$class
      mean(rda_class != test_y)
  })

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
      Guo_out <- scrda_train(x = train_x, y = train_y, prior = prior_probs)
      Guo_pred <- scrda_predict(rda_out = Guo_out$rda_out,
                                train_x = train_x,
                                train_y = train_y,
                                test_x = test_x,
                                alpha = Guo_out$alpha,
                                delta = Guo_out$delta)
      mean(Guo_pred != test_y)
  })

  # HDRDA - Ridge
  hdrda_ridge_errors <- try({
      cv_out <- hdrda_cv(x = train_x,
                         y = train_y,
                         prior = prior_probs)
      hdrda_ridge <- list(lambda = cv_out$lambda, gamma = cv_out$gamma)
      hdrda_ridge_out <- hdrda(x = train_x,
                               y = train_y,
                               lambda = cv_out$lambda,
                               gamma = cv_out$gamma,
                               prior = prior_probs)
      mean(predict(hdrda_ridge_out, test_x)$class != test_y)
  })

  # HDRDA - Convex
  hdrda_convex_errors <- try({
      cv_out <- hdrda_cv(x = train_x,
                         y = train_y,
                         prior = prior_probs,
                         shrinkage_type = "convex")
      hdrda_convex <- list(lambda = cv_out$lambda, gamma = cv_out$gamma)
      hdrda_convex_out <- hdrda(x = train_x,
                                y = train_y,
                                lambda = cv_out$lambda,
                                gamma = cv_out$gamma,
                                prior = prior_probs,
                                shrinkage_type = "convex")
      mean(predict(hdrda_convex_out, test_x)$class != test_y)
  })

  # Tong, Chen, and Zhao (2012) - Bioinformatics
  Tong_errors <- try({
      Tong_out <- dlda(x = train_x,
                       y = train_y,
                       est_mean = "tong",
                       prior = prior_probs)
      mean(predict(Tong_out, test_x)$class != test_y)
  })

  list(errors = list(
         HDRDA_Ridge = hdrda_ridge_errors,
         HDRDA_Convex = hdrda_convex_errors,
         Guo = Guo_errors,
         Tong = Tong_errors,
         Witten = Witten_errors,
         RDA = rda_errors,
         kNN = knn_errors,
         SVM_Radial = ksvm_radial_errors
    ), data_set = data_set, hdrda_ridge = hdrda_ridge, hdrda_convex = hdrda_convex)
}, mc.cores = num_cores)

save(results, file = 'data/results-microarray.RData')

