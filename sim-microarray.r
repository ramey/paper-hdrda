library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 1

data_sets <- c("burczynski", "christensen", "khan", "nakayama", "sorlie", "su", "sun", "yeoh")


results <- mclapply(data_sets, function(data_set) {
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)

  cv_folds <- cv_partition(y = data$y, num_folds = 10)

  cv_results <- lapply(cv_folds, function(cv_fold) {
    train_x <- as.matrix(data$x[cv_fold$training, ])
    train_y <- data$y[cv_fold$training]
    test_x <- as.matrix(data$x[cv_fold$test, ])
    test_y <- data$y[cv_fold$test]

    # Sets the prior probabilities equal
    num_classes <- nlevels(train_y)
    prior_probs <- rep(1, num_classes) / num_classes

    # Witten and Tibshirani (2011) - JRSS B
    Witten_out <- Witten_Tibshirani(train_x, train_y, test_x)
    Witten_errors <- sum(Witten_out$predictions != test_y)

    # Applies the variable selection from Witten and Tibshirani
    train_x <- train_x[, Witten_out$variables]
    test_x <- test_x[, Witten_out$variables]

    # Hard-ETRDA
    cv_out <- etrda_cv(x = train_x, y = train_y, threshold = "hard", num_delta = 10, prior = prior_probs)
    etrda_out <- with(cv_out, etrda(x = train_x, y = train_y, threshold = "hard", lambda = lambda, delta = delta, prior = prior_probs))
    etrda_hard_errors <- sum(predict(etrda_out, test_x)$class != test_y)

    # RDA from Friedman (1989)
    cv_out <- rda_cv(x = train_x, y = train_y, prior = prior_probs)
    rda_out <- with(cv_out, rda(x = train_x, y = train_y, lambda = lambda, gamma = gamma, prior = prior_probs))
    rda_errors <- sum(predict(rda_out, test_x)$class != test_y)
    
    # DLDA and DQDA
    dlda_errors <- sum(predict(dlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y)
    dqda_errors <- sum(predict(dqda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y)

    list(
      dlda = dlda_errors,
      dqda = dqda_errors,
      etrda_hard = etrda_hard_errors,
      rda = rda_errors,
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



