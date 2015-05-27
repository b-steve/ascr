pkgs <- c("CircStats", "fastGHQuad", "knitr", "matrixStats", "mvtnorm", "optimx", "plyr",
          "Rcpp", "R2admb", "secr", "testthat", "truncnorm", "xtable", "downloader")
options(warn = -1)
for (i in pkgs){
    if (!require(i, quietly = TRUE, character.only = TRUE)){
        install.packages(i)
    }
}
options(warn = 0)
if (.Platform$OS == "windows"){
    bin.name <- "https://github.com/b-steve/admbsecr/releases/download/v1.1.2/admbsecr_1.1.2.zip"
    ext <- ".zip"
    type <- "win.binary"
} else if (.Platform$OS == "unix"){
    bin.name <- "https://github.com/b-steve/admbsecr/archive/v1.1.2.tar.gz"
    ext <- ".tar.gz"
    type <- "source"
} else {
    stop("Unknown OS type.")
}
dest <- paste("admbsecr_1.1.2", ext, sep = "")
library(downloader)
download(bin.name, destfile = dest)
install.packages(dest, repos = NULL, type = type)
unlink(dest)



