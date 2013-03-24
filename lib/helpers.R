#' Generates the Cartesian product of all simulation experiment factors.
#' Then, replicates the Cartesian product a specified number of times.
#' This function is especially useful for simulations run with parallel
#' processing, where we need to run the same set of simulation configurations
#' multiples times.
#'
#' If num_reps = 1, the output is nearly identical to that of the 'expand.grid' function.
#' In fact, much of the documentation for this function has been copied
#' directly from the 'expand.grid' function.
#'
#' @param ... vectors, factors or a list containing these.
#' @param num_reps the number of times to replicate each simulation configuration
#' @param stringsAsFactors logical specifying if character vectors are converted
#'   to factors.
#' @return A data frame containing D rows for each combination of the supplied factors.
#'  The first factors vary fastest. The columns are labelled by the factors if these
#'  are supplied as named arguments or named components of a list. The row names are
#'  'automatic.'
sim_factors <- function(..., num_reps = 1, stringsAsFactors = FALSE) {
  grid <- expand.grid(..., stringsAsFactors = stringsAsFactors)
  do.call(rbind, replicate(num_reps, grid, simplify = FALSE))
}

#' Random partition of classification data.
#'
#' For a vector of training labels, we return a list of random partitions,
#' where each partition has the indices of both the training and test
#' observations.
#'
#' We partition the vector \code{y} based on its length, which we treat as the
#' sample size, \code{n}. If an object other than a vector is used in \code{y},
#' its length can yield unexpected results. For example, the output of
#' \code{length(diag(3))} is 9.
#'
#' @export
#' @param y a vector of the labels of the training data
#' @return list the indices of the training and test observations for each fold.
#' @examples
#' library(MASS)
#' # The following three calls to \code{rand_partition} yield the same partitions.
#' set.seed(42)
#' rand_partition(iris$Species)
#' rand_partition(iris$Species, num_partitions = 5)
#' rand_partitions(iris$Species, train_pct = 2/3)
rand_partition <- function(y, num_partitions = 10, train_pct = 1/2) {
  n <- length(y)
  seq_n <- seq_len(n)

  splits <- replicate(num_partitions, expr = {
    training <- sample(seq_n, train_pct * n)
    test <- which(!seq_n %in% training)
    list(
      training = training,
      test = test
    )
  }, simplify = FALSE)

  names(splits) <- paste0("Split", seq_len(num_partitions))
  splits
}

# This function rounds to the nearest digits and does not truncated leading
# zeroes, unlike the 'round' function.
round_digits <- function(x, num_digits = 3) {
  sprintf(paste("%.", num_digits, "f", sep = ""), round(x, num_digits))
}

# Classifier from Witten and Tibshirani (2011) - JRSS B
Witten_Tibshirani <- function(train_x, train_y, test_x) {
  # Uses cross-validation to select the best lambda and K
  cv_out <- PenalizedLDA.cv(x = train_x, y = as.integer(train_y))

  # Classifies the test data set.
  PLDA_out <- with(cv_out, PenalizedLDA(x = train_x, y = as.integer(train_y),
                                        xte = test_x, lambda = bestlambda,
                                        K = bestK))

  # The variables selected are those with values of 'beta' that are nonzero.
  variables_selected <- which(with(PLDA_out, discrim[, K]) != 0)

  # Obtains the classifier's predictions on the test data
  # Note that the test predictions are integers and not factors. We convert the
  # test predictions back to the factors using the levels of the training data set.
  test_predictions <- factor(with(PLDA_out, ypred[, K]),
                             levels = seq_len(nlevels(train_y)),
                             labels = levels(train_y))

  list(predictions = test_predictions, variables = variables_selected)
}

# Helper function to load the appropriate microarray data set.
# If the Nakayama data set is specified, we restrict it to the classes with
# at least 15 samples.
load_microarray <- function(dataset) {
  # Reduces the 'nakayama' data set to classes with sample sizes >= 15
  # Witten and Tibshirani (2011) do this as well
  if (dataset == 'nakayama') {
    data('nakayama', package = 'datamicroarray')
    classes <- with(nakayama, levels(y)[table(y) >= 15])
    idx <- which(nakayama$y %in% classes)
    nakayama$x <- nakayama$x[idx, ]
    nakayama$y <- factor(nakayama$y[idx])
    data_out <- nakayama
  } else if (dataset == 'su') {
    # The Witten and Tibshirani classifier results in an error from features that
    # have within-class standard deviation of 0. I double-checked the data, and it
    # turns out there are features within a class that have the same value observed,
    # hence the error. I resolve the issue by removing all features that have 0
    # variance before I begin the simulation. Here is the error from Witten's code:
    # Error in PenalizedLDA(x = xtr, y = ytr, xte = xte, lambda = lambdas[i],  :
    #   Some features have 0 within-class standard deviation.
    data('su', package = 'datamicroarray')
    class_var <- tapply(seq_along(su$y), su$y, function(i) apply(su$x[i,], 2, var))
    class_var <- do.call(rbind, class_var)
    features_kept <- which(apply(class_var, 2, function(x) all(x != 0)))

    su$x <- su$x[, features_kept]
    data_out <- su

  } else {
    data(list = dataset, package = "datamicroarray")
    data_out <- get(dataset)
  }
  data_out
}

# Variable selection method from Dudoit et al.
dudoit <- function(train_x, train_y, test_x = NULL, q) {
  library('multtest')

  F_stat <- mt.teststat(t(train_x), train_y, test = "f")
  top_q_genes <- rev(order(F_stat))[seq_len(q)]

  out <- list()
  out$train_x <- train_x[, top_q_genes]

  if (!is.null(test_x)) {
    out$test_x <- test_x[, top_q_genes]
  }

  out
}

# Classifier from Clemmensen, Hastie, Witten and Ersbøll (2012) - Technometrics
Clemmensen <- function(train_x, train_y, test_x, normalize_data = FALSE, cv_variables = FALSE, num_folds = 10, ...) {
  if (normalize_data) {
    normalize_out <- normalize(train_x)
    train_x <- normalize_out$Xc
    test_x <- normalizetest(test_x, normalize_out)
  }

  # Selects the number of variables by cross-validation.
  # The candidate values are 10 to p in increments of 10.
  # If there is a tie, we choose the smallest value.
  # The code often results in the following error:
  # Error in lda.default(x, grouping, ...)
  #  variable 1 appears to be constant within groups
  # We ignore any cases where this occurs.
  if (cv_variables) {
    num_vars <- seq.int(10, ncol(train_x), by = 10)
    cv_folds <- cv_partition(y = train_y, num_folds = num_folds)

    # Traverse through each cross-validation fold and compute the number of
    # cross-validation errors for each reduced dimension
    cv_errors <- lapply(cv_folds, function(cv_fold) {
      trn_x <- train_x[cv_fold$training, ]
      trn_y <- train_y[cv_fold$training]
      tst_x <- train_x[cv_fold$test, ]
      tst_y <- train_y[cv_fold$test]

      # For each reduced dimension considered, we calculate the number of test errors
      # resulting for the current cross-validation fold.
      cv_num_vars <- sapply(num_vars, function(q) {
        sda_predict <- try_default(predict(sda(x = trn_x, y = trn_y, stop = -q), tst_x)$class, NA, quiet = TRUE)
        sum(sda_predict != tst_y)
      })
    })
    cv_errors <- colSums(do.call(rbind, cv_errors))

    # Determines the optimal value of 'num_vars' to be the one that yielded the minimized
    # the cross-validation error rate. If there is a tie, we break the tie by
    # choosing the smallest value of 'num_vars' for parsimony.
    num_vars <- num_vars[which.min(cv_errors)]
  } else {
    num_vars <- ncol(train_x)
  }
     
  predict(sda(x = train_x, y = train_y, stop = -num_vars), test_x)$class
}

# Classifier from Cao, Boitard, and Besse (2011) - BMC Bioinformatics
# We use the 'max.dist' 
Cao <- function(train_x, train_y, test_x) {
  plsda_out <- plsda(X = train_x, Y = train_y)
  plsda_predict <- predict(plsda_out, test_x, method = "max.dist")
  factor(plsda_predict$class$max.dist[, 2], levels = seq_len(nlevels(train_y)),
         labels = levels(train_y))
}


# Classifier from Guo, Hastie, and Tibshirani (2007) - Biostatistics
scrda_train <- function(x, y) {
  x <- t(x)
  rda_out <- rda(x = x, y = y)
  rda_cv_out <- rda.cv(rda_out, x = x, y = y)

  # Which alpha and delta give min cv error?
  min_cv_error <- with(rda_cv_out, which(cv.err == min(cv.err), arr.ind = TRUE))

  # Out of the minima, what are the gene dimensions?
  num_genes <- rda_cv_out$ngene[min_cv_error]

  # Which alpha-delta has smallest dimension? If there's still a tie, choose one
  # at random.
  min_alpha_delta <- sample(which(num_genes == min(num_genes)), size = 1)

  # Now we determine the corresponding alpha-delta pair.
  alpha_delta <- min_cv_error[min_alpha_delta,]
  alpha <- as.numeric(rownames(rda_cv_out$cv.err)[alpha_delta[1]])
  delta <- as.numeric(colnames(rda_cv_out$cv.err)[alpha_delta[2]])

  list(rda_out = rda_out, rda_cv_out = rda_cv_out, alpha = alpha, delta = delta)
}

# Classifier from Guo, Hastie, and Tibshirani (2007) - Biostatistics
scrda_predict <- function(rda_out, train_x, train_y, test_x, alpha, delta) {
  group_names <- levels(train_y)
  factor(
    predict(rda_out, x = t(train_x), y = train_y, xnew = t(test_x), alpha = alpha, delta = delta),
    levels = seq_along(group_names), labels = group_names)
}

