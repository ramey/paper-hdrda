library('ProjectTemplate')
load.project()

data(alon)
q_vals <- seq.int(50, 2000, by = 50)

log_det <- lapply(q_vals, function(q) {
  x <- with(alon, dudoit(x, y, q = q))$train_x
  Sigma <- cov_pool(x = x, y = alon$y)
  Sigma_eigenvalues <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values

  shrinkage_rda <- mean(Sigma_eigenvalues)
  shrinkage_mdeb <- sum(Sigma_eigenvalues) / length(alon$y)
  
  eigen_rda <- Sigma_eigenvalues + shrinkage_rda
  eigen_mdeb <- Sigma_eigenvalues + shrinkage_mdeb
  eigen_nlda <- replace(Sigma_eigenvalues, Sigma_eigenvalues < shrinkage_rda, shrinkage_rda)
  
  log_determinants <- list(MDEB = sum(log(eigen_mdeb)),
                           NLDA = sum(log(eigen_nlda)),
                           RDA = sum(log(eigen_rda)),
                           q = q)
  log_determinants
})
log_det <- do.call(rbind.data.frame, log_det)
m_log_det <- melt(log_det, id.vars = 'q', variable.name = 'Classifier',
                  value.name = 'Log_Determinant')

# Log-determinants of the covariance matrices for each of the three shrinkage methods
# We used the pooled sample covariance matrix
# In the figure caption, say that this corresponds to $\lambda = 1$
p <- ggplot(m_log_det, aes(x = q, y = Log_Determinant, group = Classifier, color = Classifier))
p + geom_line(aes(linetype = Classifier))
