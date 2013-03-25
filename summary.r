library('ProjectTemplate')
load.project()

# Remove the following data sets because there is too much separation for a good comparison
results <- na.omit(results)

results$d <- as.factor(results$d)
results$data_set <- as.character(results$data_set)

p <- ggplot(results, aes(x = d, y = Error))
p <- p + geom_boxplot(aes(fill = Classifier))
p + facet_wrap(~ data_set)

results_summary <- ddply(results, .(data_set, d, Classifier), summarize,
                         Avg_Error = mean(Error),
                         SE_Error = sqrt(mean(Error) * (1 - mean(Error)) / length(Error)))

# Table containing a summary of the data sets
data_summary <- subset(describe_data(), author %in% unique(results$data_set))
data_summary$K[data_summary$author == "nakayama"] <- 5
rownames(data_summary) <- NULL
colnames(data_summary)[1:3] <- c("Author", "Year", "N")
data_summary$Author <- Hmisc:::capitalize(as.character(data_summary$Author))
data_summary$Author <- with(data_summary, paste0(Author, " et al. (", Year, ")"))
data_summary <- subset(data_summary, select = -Year)
print(xtable(data_summary, caption = "Summary of high-dimensional microarray data sets.",
             label = "tab:data-summary", digits = 0), include.rownames = FALSE)


# Tables containing the classification results for each value of 'q'
results_summary$data_set <- Hmisc:::capitalize(results_summary$data_set)
results_summary$d <- as.factor(results_summary$d)
results_summary$Avg_Error <- round_digits(results_summary$Avg_Error, 3)
results_summary$SE_Error <- round_digits(results_summary$SE_Error, 3)
results_summary$Error <- with(results_summary, paste0(Avg_Error, " (", SE_Error, ")"))
results_summary <- subset(results_summary, select = -c(Avg_Error, SE_Error))
colnames(results_summary) <- c("Dataset", "d", "Classifier", "Error")

results_tables <- tapply(seq_along(results_summary$d), results_summary$d, function(i) {
  dcast(results_summary[i, ], Dataset ~ Classifier, value.var = "Error")
})

results_xtables <- lapply(names(results_tables), function(d) {
  tab_caption <- "The average of the test error rates for the microarray data sets with $d = "
  tab_caption <- paste0(tab_caption, d, "$. Approximate standard errors are given in parentheses.")
  xtable(results_tables[[d]], caption = tab_caption, label = paste0("tab:results-d", d))
})
