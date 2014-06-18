install.packages(c("CircStats", "lattice", "matrixStats", "plyr",
                   "Rcpp", "R2admb", "secr", "testthat", "downloader"))
library(downloader)
download("https://github.com/b-steve/admbsecr/releases/download/v1.0.6/admbsecr_1.0.6-win.zip", destfile = "admbsecr_1.0.6.zip")
install.packages("admbsecr_1.0.6.zip", repos = NULL)
unlink("admbsecr_1.0.6.zip")
