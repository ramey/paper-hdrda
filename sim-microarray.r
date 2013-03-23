library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 12

# Problem Data Sets
# burczynski
# chiaretti
data_sets <- c("alon", "burczynski", "chiaretti", "christensen", "gordon",
               "gravier", "khan", "nakayama", "shipp", "singh", "sorlie", "yeoh")

results <- mclapply(data_sets, function(data_set) {
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)

  set.seed(42)

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

    # Witten and Tibshirani (2011) - JRSS B
    Witten_out <- try_default(Witten_Tibshirani(train_x, train_y, test_x), NA)
    Witten_errors <- sum(Witten_out$predictions != test_y)

    # Guo, Hastie, and Tibshirani (2007) - Biostatistics
    Guo_out <- try_default(scrda_train(x = train_x, y = train_y), NA)
    Guo_pred <- try_default(with(Guo_out, scrda_predict(rda_out, train_x, train_y, test_x,
                                            alpha, delta)), NA)
    Guo_errors <- sum(Guo_pred != test_y)

    # GRDA
    cv_out <- try_default(grda_cv(x = train_x, y = train_y, prior = prior_probs), NA)
    grda_out <- try_default(with(cv_out, grda(x = train_x, y = train_y, lambda = lambda, gamma = gamma, prior = prior_probs)), NA)
    grda_errors <- try_default(sum(predict(grda_out, test_x)$class != test_y), NA)

    # DLDA and DQDA
    dlda_errors <- try_default(sum(predict(dlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)
    dqda_errors <- try_default(sum(predict(dqda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)

    # Pang, Tong, and Zhao (2009) - Biometrics
    sdlda_errors <- try_default(sum(predict(sdlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)
    sdqda_errors <- try_default(sum(predict(sdqda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)

    # Tong, Chen, and Zhao (2012) - Bioinformatics
    Tong_out <- try_default(dlda(x = train_x, y = train_y, est_mean = "tong", prior = prior_probs), NA)
    Tong_errors <- try_default(sum(predict(Tong_out, test_x)$class != test_y), NA)

    # Cao, Boitard, and Besse (2011) - BMC Bioinformatics
    Cao_errors <- try_default(sum(Cao(train_x, train_y, test_x) != test_y), NA)

    list(
      Cao = Cao_errors,
      dlda = dlda_errors,
      dqda = dqda_errors,
      grda = grda_errors,
      Guo = Guo_errors,
      sdlda = sdlda_errors,
      sdqda = sdqda_errors,
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

