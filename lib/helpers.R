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

#' Helper function to return predictions from the list element 'class'
#'
#' A common convention for classifiers in R is to return a list from the
#' 'predict' with the classifications (class labels) of each test observation
#' along with a variety of information about each specific test
#' classification, such as discriminant scores and posterior probabilities.
#' For our 'classify' function, we need only the class labels. This helper
#' function is a wrapper function that returns only the class labels when
#' the aforementioned convention is followed.
predict_helper <- function(obj, newdata) {
  predict(obj, newdata)$class
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

