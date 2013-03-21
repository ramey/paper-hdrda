library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 12

data_sets <- c("alon", "burczynski", "chiaretti", "christensen", "gordon",
               "gravier", "khan", "nakayama", "shipp", "singh", "sorlie", "su",
               "yeoh")

results <- mclapply(data_sets, function(data_set) {
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)

  cv_folds <- cv_partition(y = data$y, num_folds = 10)

  cv_results <- lapply(seq_along(cv_folds), function(i) {
    message("Dataset: ", data_set, " -- Fold ", i)
    cv_fold <- cv_folds[[i]]
    train_x <- as.matrix(data$x[cv_fold$training, ])
    train_y <- data$y[cv_fold$training]
    test_x <- as.matrix(data$x[cv_fold$test, ])
    test_y <- data$y[cv_fold$test]

    # Sets the prior probabilities equal
    num_classes <- nlevels(train_y)
    prior_probs <- rep(1, num_classes) / num_classes

    # For now, I've commented out the Technometrics method because it's too slow.
    # As in, it took about 12 hours for a single fold in my 10-fold CV simulations
    # of several data sets.
    # Clemmensen, Hastie, Witten and Ersbøll (2012) - Technometrics
    # Clemmensen_out <- Clemmensen(train_x, train_y, test_x, cv_variables = TRUE,
    #                             normalize_data = TRUE)
    # Clemmensen_errors <- sum(Clemmensen_out != test_y)

    # Witten and Tibshirani (2011) - JRSS B
    Witten_out <- Witten_Tibshirani(train_x, train_y, test_x)
    Witten_errors <- sum(Witten_out$predictions != test_y)

    # Applies the variable selection from Witten and Tibshirani (2011)
    # train_x <- train_x[, Witten_out$variables]
    # test_x <- test_x[, Witten_out$variables]

    # Classifier from Guo, Hastie, and Tibshirani (2007) - Biostatistics
    Guo_out <- scrda_train(x = train_x, y = train_y)
    Guo_pred <- with(Guo_out, scrda_predict(rda_out, train_x, train_y, test_x,
                                            alpha, delta))
    Guo_errors <- sum(Guo_pred != test_y)

    # GRDA
    cv_out <- grda_cv(x = train_x, y = train_y, prior = prior_probs)
    grda_out <- with(cv_out, grda(x = train_x, y = train_y, lambda = lambda, gamma = gamma, prior = prior_probs))
    grda_errors <- sum(predict(grda_out, test_x)$class != test_y)

    # DLDA and DQDA
    dlda_errors <- sum(predict(dlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y)
    dqda_errors <- sum(predict(dqda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y)

    # Tong et al. (2012) - Bioinformatics
    Tong_out <- dlda(x = train_x, y = train_y, est_mean = "tong", prior = prior_probs)
    Tong_errors <- sum(predict(Tong_out, test_x)$class != test_y)

    # Cao, Boitard, and Besse (2011) - BMC Bioinformatics
    Cao_errors <- sum(Cao(train_x, train_y, test_x) != test_y)

    list(
      Cao = Cao_errors,
      dlda = dlda_errors,
      dqda = dqda_errors,
      grda = grda_errors,
      Guo = Guo_errors,
      Tong = Tong_errors,
      Witten = Witten_errors
      )
  })
  melt_results <- melt(cv_results)
  results_errors <- ddply(melt_results, .(L2), summarize, errors = sum(value))
  colnames(results_errors) <- c("Classifier", "errors")
  results_errors$error_rate <- results_errors$errors / length(data$y)
  cbind(data = data_set, results_errors)
}, mc.cores = num_cores)
results <- do.call(rbind, results)

save(results, file = 'results-microarray.RData')

