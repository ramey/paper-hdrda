library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 12

train_pct <- 2/3
num_iterations <- 1000
data_sets <- c("burczynski", "nakayama", "shipp", "singh")
d_vals <- 100 * seq_len(10)

sim_config <- sim_factors(data_sets = data_sets, num_reps = num_iterations,
                          stringsAsFactors = FALSE)

results <- mclapply(seq_len(nrow(sim_config)), function(i) {
  data_set <- sim_config$data_set[i]
  
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)

  rand_split <- rand_partition(y = data$y, num_partitions = 1,
                               train_pct = train_pct)[[1]]

  message("Dataset: ", data_set, " -- Sim Config: ", i, " of ", nrow(sim_config))
  train_x <- as.matrix(data$x[rand_split$training, ])
  train_y <- data$y[rand_split$training]
  test_x <- as.matrix(data$x[rand_split$test, ])
  test_y <- data$y[rand_split$test]

  # Apply Dudoit variable selection to rank the genes in descending order
  dudoit_out <- dudoit(train_x, train_y)

  # For each value of 'd', we compute the error-rate estimates for each of the
  # classifiers.
  errors_d <- lapply(d_vals, function(d) {
    trn_x <- train_x[, head(dudoit_out, d)]
    tst_x <- test_x[, head(dudoit_out, d)]

    # Sets the prior probabilities equal
    num_classes <- nlevels(train_y)
    prior_probs <- rep(1, num_classes) / num_classes

    # Witten and Tibshirani (2011) - JRSS B
    Witten_out <- try_default(Witten_Tibshirani(trn_x, train_y, tst_x), NA)
    Witten_errors <- mean(Witten_out$predictions != test_y)

    # TODO: Double-check cross-validation method to determine why method is classifying
    # into one class and therefore having terrible classification accuracy
    # Guo, Hastie, and Tibshirani (2007) - Biostatistics
    Guo_out <- try_default(scrda_train(x = trn_x, y = train_y), NA)
    Guo_pred <- try_default(with(Guo_out, scrda_predict(rda_out, trn_x, train_y, tst_x,
                                             alpha, delta)), NA)
    Guo_errors <- mean(Guo_pred != test_y)

    # HDRDA
    cv_out <- try_default(hdrda_cv(x = trn_x, y = train_y, prior = prior_probs), NA)
    hdrda_lambda <- cv_out$lambda
    hdrda_gamma <- cv_out$gamma
    hdrda_out <- try_default(with(cv_out, hdrda(x = trn_x, y = train_y, lambda = lambda, gamma = gamma, prior = prior_probs)), NA)
    hdrda_errors <- try_default(mean(predict(hdrda_out, tst_x)$class != test_y), NA)

    # DLDA
    dlda_errors <- try_default(mean(predict(dlda(x = trn_x, y = train_y, prior = prior_probs), tst_x)$class != test_y), NA)

    # Pang, Tong, and Zhao (2009) - Biometrics
    sdlda_errors <- try_default(mean(predict(sdlda(x = trn_x, y = train_y, prior = prior_probs), tst_x)$class != test_y), NA)

    # Tong, Chen, and Zhao (2012) - Bioinformatics
    Tong_out <- try_default(dlda(x = trn_x, y = train_y, est_mean = "tong", prior = prior_probs), NA)
    Tong_errors <- try_default(mean(predict(Tong_out, tst_x)$class != test_y), NA)

    list(errors = list(
           Dudoit = dlda_errors,
           HDRDA = hdrda_errors,
           Guo = Guo_errors,
           Pang = sdlda_errors,
           Tong = Tong_errors,
           Witten = Witten_errors
      ), data_set = data_set, d = d, lambda = hdrda_lambda, gamma = hdrda_gamma)
  })
}, mc.cores = num_cores)

save(results, file = 'data/results-microarray.RData')

