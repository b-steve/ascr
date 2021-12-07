#' @export
compile.ascr <- function(dev = FALSE){
  if(dev){
    dir = paste0(getwd(), '/inst/TMB/')
  } else {
    dir = paste0(system.file(package = "ascr"), "/TMB/")
  }

  
  TMB::compile(paste0(dir, "ascrTmb.cpp"))

}
