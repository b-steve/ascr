pkgs <- c("CircStats", "lattice", "matrixStats", "plyr",
          "Rcpp", "R2admb", "secr", "testthat", "downloader")
for (i in pkgs){
    if (!require(i, quietly = TRUE, character.only = TRUE)){
        install.packages(i)
    }
}
if (.Platform$OS == "windows"){
    bin.name <- "https://github.com/b-steve/admbsecr/releases/download/v1.0.6/admbsecr_1.0.6.zip"
    ext <- ".zip"
} else if (.Platform$OS == "unix"){
    bin.name <- "https://github.com/b-steve/admbsecr/archive/v1.0.6.tar.gz"
    ext <- ".tar.gz"
} else {
    stop("Unknown OS type.")
}
dest <- paste("admbsecr_1.0.6", ext, sep = "")
library(downloader)
download(bin.name, destfile = dest)
install.packages(dest, repos = NULL)
unlink(dest)
