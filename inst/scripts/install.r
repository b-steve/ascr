pkgs <- c("CircStats", "fastGHQuad", "fields", "knitr", "matrixStats", "mgcv", "mvtnorm", "optimx", "plyr",
          "Rcpp", "R2admb", "secr", "testthat", "truncnorm", "viridis", "xtable", "downloader")
options(warn = -1)
for (i in pkgs){
    if (!require(i, quietly = TRUE, character.only = TRUE)){
        install.packages(i)
    }
}
options(warn = 0)
if (.Platform$OS == "windows"){
    bin.name <- "https://github.com/b-steve/ascr/releases/download/v2.2.3/ascr_2.2.3.zip"
    ext <- ".zip"
    type <- "win.binary"
} else if (.Platform$OS == "unix"){
    bin.name <- "https://github.com/b-steve/ascr/archive/v2.2.3.tar.gz"
    ext <- ".tar.gz"
    type <- "source"
} else {
    stop("Unknown OS type.")
}
dest <- paste("ascr_2.2.3", ext, sep = "")
library(downloader)
download(bin.name, destfile = dest)
install.packages(dest, repos = NULL, type = type)
unlink(dest)



