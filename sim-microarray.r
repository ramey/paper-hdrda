library(ProjectTemplate)
load.project()

set.seed(42)
num_cores <- 4

train_pct <- 3/4
num_iterations <- 50

data_sets <- c("burczynski", "nakayama", "shipp", "singh")
sim_config <- rep(data_sets, each=num_iterations)

flog.logger("sim", INFO, appender=appender.file('sim-microarray.log'))

results <- mclapply(seq_along(sim_config), function(i) {
  # Load data and randomly partition into training/test data sets.
  data_set <- sim_config[i]
  flog.info("Dataset: %s -- Sim Config: %s of %s", data_set, i, length(sim_config), name="sim")
  set.seed(i)
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)
  rand_split <- rand_partition(y=data$y, num_partitions=1,
                               train_pct=train_pct)[[1]]
  train_x <- as.matrix(data$x[rand_split$training, ])
  train_y <- data$y[rand_split$training]
  test_x <- as.matrix(data$x[rand_split$test, ])
  test_y <- data$y[rand_split$test]

  # Removes the loaded data from memory and then garbage collects.
  rm(data)
  gc(reset=TRUE)

  # Sets the prior probabilities as equal
  num_classes <- nlevels(train_y)
  prior_probs <- rep(1, num_classes) / num_classes

  # kNN
  knn_errors <- try({
      tune.knn_out <- tune.knn(x=train_x,
                               y=train_y,
                               k=1:5,
                               tunecontrol=tune.control(sampling="cross", cross=10))
      knn_out <- knn(train=train_x,
                     cl=train_y,
                     test=test_x,
                     k=tune.knn_out$best.model$k)
      mean(knn_out != test_y)
  })

  # Random Forest
  rf_errors <- try({
    rf_out <- randomForest(x=train_x,
                           y=train_y,
                           ntree=250,
                           maxnodes=100)
    rf_predict <- predict(rf_out, test_x)
    mean(rf_predict != test_y)
  })

  # SVM with Radial Basis Functions
  ksvm_radial_errors <- try({
      ksvm_out <- ksvm(x=train_x,
                       y=train_y,
                       kernel="rbfdot",
                       kpar="automatic",
                       cross=10)
      ksvm_radial <- predict(ksvm_out, test_x)
      mean(ksvm_radial != test_y)
  })

  # Witten and Tibshirani (2011) - JRSS B
  Witten_errors <- try({
      Witten_out <- Witten_Tibshirani(train_x=train_x,
                                      train_y=train_y,
                                      test_x=test_x)
      mean(Witten_out$predictions != test_y)
  })

  # Guo, Hastie, and Tibshirani (2007) - Biostatistics
  Guo_errors <- try({
      Guo_out <- Guo(train_x=train_x,
                     train_y=train_y,
                     test_x=test_x,
                     prior=prior_probs)
      mean(Guo_out != test_y)
  })

  # HDRDA - Ridge
  hdrda_ridge_errors <- try({
      cv_out <- hdrda_cv(x=train_x,
                         y=train_y,
                         prior=prior_probs)
      hdrda_ridge <- list(lambda=cv_out$lambda, gamma=cv_out$gamma)
      flog.info("HDRDA Ridge. Lambda: %s. Gamma: %s", cv_out$lambda, cv_out$gamma, name="sim")
      mean(predict(cv_out, test_x)$class != test_y)
  })

  # HDRDA - Convex
  hdrda_convex_errors <- try({
      cv_out <- hdrda_cv(x=train_x,
                         y=train_y,
                         prior=prior_probs,
                         num_gamma=21,
                         shrinkage_type="convex")
      hdrda_convex <- list(lambda=cv_out$lambda, gamma=cv_out$gamma)
      flog.info("HDRDA Convex. Lambda: %s. Gamma: %s", cv_out$lambda, cv_out$gamma, name="sim")
      mean(predict(cv_out, test_x)$class != test_y)
  })

  # Tong, Chen, and Zhao (2012) - Bioinformatics
  Tong_errors <- try({
      Tong_out <- dlda(x=train_x,
                       y=train_y,
                       est_mean="tong",
                       prior=prior_probs)
      mean(predict(Tong_out, test_x)$class != test_y)
  })

  # Pang, Tong, and Zhao (2009) - Bioinformatics
  Pang_errors <- try({
      Pang_out <- sdlda(x=train_x,
                        y=train_y,
                        prior=prior_probs)
      mean(predict(Pang_out, test_x)$class != test_y)
  })

  error_rates <- list(
    Guo=Guo_errors,
    HDRDA_Ridge=hdrda_ridge_errors,
    HDRDA_Convex=hdrda_convex_errors,
    kNN=knn_errors,
    Pang=Pang_errors,
    Random_Forest=rf_errors,
    SVM_Radial=ksvm_radial_errors,
    Tong=Tong_errors,
    Witten=Witten_errors
  )

  list(error_rates=error_rates,
       data_set=data_set,
       hdrda_ridge=hdrda_ridge,
       hdrda_convex=hdrda_convex
  )
}, mc.cores=num_cores)

results <- list(sim_results=results,
                train_pct=train_pct,
                num_iterations=num_iterations)

save(results, file='data/results-microarray.RData')
