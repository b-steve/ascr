#' Run admbsecr tests
#'
#' Runs tests for admbsecr.
#'
#' This function was written to allow users to test their admbsecr
#' installation. Users cannot use the \link[testthat]{test_package}
#' function from the testthat package as the tests are not
#' installed. This is because tests fail on the R-forge servers due to
#' the AD model builder executable.
#'
#' @export
test.admbsecr <- function(){
    if (require(testthat)){
        dir <- paste(system.file(package = "admbsecr"), "testthat", sep = "/")
        test_dir(dir)
    } else {
        stop("Please install package 'testthat' from CRAN.")
    }
}
