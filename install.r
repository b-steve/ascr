install.packages(c("CircStats", "lattice", "matrixStats", "plyr",
                   "Rcpp", "R2admb", "secr", "testthat", "downloader"))
if (.Platform$OS == "windows"){
    os <- "win"
    ext <- ".zip"
} else if (.Platform$OS == "unix"){
    if (Sys.info()["sysname"] == "Linux"){
        os <- "linux"
        ext <- ".tar.gz"
    } else if (Sys.info()["sysname"] == "Darwin"){
        os <- "osx"
        ext <- ".tgz"
    } else {
        stop("Unknown OS type.")
    }
} else {
    stop("Unknown OS type.")
}
bin.name <- paste("https://github.com/b-steve/admbsecr/releases/download/v1.0.6/admbsecr_1.0.6-",
                  os, ext, sep = "")
dest <- paste("admbsecr_1.0.6", ext, sep = "")
library(downloader)
download(bin.name, destfile = dest)
install.packages(dest, repos = NULL)
unlink(dest)
