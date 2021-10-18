compile.ascr <- function(template = 'ascrTmb.cpp'){
  wd <- getwd()
  dir <- paste0(
    #system.file(package = "ascr"), 
    
    #for dev only
    wd,
    "/src")
  setwd(dir)
  TMB::compile(template)
  setwd(wd)
}