source("http://bioconductor.org/biocLite.R")
biocLite("multtest", ask=FALSE)

cran_packages <- c('futile.logger', 'reshape2', 'plyr', 'ggplot2', 'stringr',
                   'rda', 'penalizedLDA', 'xtable', 'Hmisc', 'mixtools',
                   'gridExtra', 'sparseLDA', 'e1071', 'klaR', 'kernlab',
                   'devtools', 'randomForest', 'ProjectTemplate', 'sortinghat')
install.packages(cran_packages, dep=TRUE)

github_packages <- c('ramhiser/datamicroarray', 'ramhiser/sparsediscrim',
                     'ramhiser/itertools2')
library(devtools)
devtools::install_github(github_packages)
