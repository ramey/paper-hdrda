library('ProjectTemplate')
load.project()

set.seed(42)
num_cores <- 12

train_pct <- 2/3
num_iterations <- 1000
data_sets <- c("burczynski", "nakayama", "shipp", "singh")
d <- c(100, 250, 500, 1000)

sim_config <- expand.grid(data_sets = data_sets, d = d, stringsAsFactors = FALSE)

results <- mclapply(seq_len(nrow(sim_config)), function(i) {
  data_set <- sim_config$data_set[i]
  d <- sim_config$d[i]
  
  data <- load_microarray(data_set)
  data$x <- as.matrix(data$x)

  rand_splits <- rand_partition(y = data$y, num_partitions = num_iterations,
                                train_pct = train_pct)

  split_results <- lapply(seq_along(rand_splits), function(split) {
    message("Dataset: ", data_set, " -- Split ", split)
    rand_split <- rand_splits[[split]]
    train_x <- as.matrix(data$x[rand_split$training, ])
    train_y <- data$y[rand_split$training]
    test_x <- as.matrix(data$x[rand_split$test, ])
    test_y <- data$y[rand_split$test]

    # Apply Dudoit variable selection
    dudoit_out <- dudoit(train_x, train_y, test_x, q = d)
    train_x <- dudoit_out$train_x
    test_x <- dudoit_out$test_x

    # Sets the prior probabilities equal
    num_classes <- nlevels(train_y)
    prior_probs <- rep(1, num_classes) / num_classes

    # Witten and Tibshirani (2011) - JRSS B
    Witten_out <- try_default(Witten_Tibshirani(train_x, train_y, test_x), NA)
    Witten_errors <- mean(Witten_out$predictions != test_y)

    # Guo, Hastie, and Tibshirani (2007) - Biostatistics
    Guo_out <- try_default(scrda_train(x = train_x, y = train_y), NA)
    Guo_pred <- try_default(with(Guo_out, scrda_predict(rda_out, train_x, train_y, test_x,
                                            alpha, delta)), NA)
    Guo_errors <- mean(Guo_pred != test_y)

    # GRDA
    cv_out <- try_default(grda_cv(x = train_x, y = train_y, prior = prior_probs), NA)
    grda_out <- try_default(with(cv_out, grda(x = train_x, y = train_y, lambda = lambda, gamma = gamma, prior = prior_probs)), NA)
    grda_errors <- try_default(mean(predict(grda_out, test_x)$class != test_y), NA)

    # DLDA
    dlda_errors <- try_default(mean(predict(dlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)

    # Pang, Tong, and Zhao (2009) - Biometrics
    sdlda_errors <- try_default(mean(predict(sdlda(x = train_x, y = train_y, prior = prior_probs), test_x)$class != test_y), NA)

    # Tong, Chen, and Zhao (2012) - Bioinformatics
    Tong_out <- try_default(dlda(x = train_x, y = train_y, est_mean = "tong", prior = prior_probs), NA)
    Tong_errors <- try_default(mean(predict(Tong_out, test_x)$class != test_y), NA)

    list(
      Dudoit = dlda_errors,
      GRDA = grda_errors,
      Guo = Guo_errors,
      Pang = sdlda_errors,
      Tong = Tong_errors,
      Witten = Witten_errors
      )
  })

  melt_results <- melt(split_results)
  melt_results <- subset(melt_results, select = -L1)
  colnames(melt_results) <- c("Error", "Classifier")
  cbind(data_set, d, melt_results)
}, mc.cores = num_cores)
results <- do.call(rbind, results)

save(results, file = 'data/results-microarray.RData')

