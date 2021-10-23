#' @export
compile.ascr <- function(){

  dir = paste0(system.file(package = "ascr"), "/TMB/")
  
  TMB::compile(paste0(dir, "ascrTmb.cpp"))

}
