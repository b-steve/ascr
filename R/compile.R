#' @export
compile.ascr <- function(dev = FALSE, debug.mode = FALSE){
  if(dev){
    dir = paste0(getwd(), '/inst/TMB/')
  } else {
    dir = paste0(system.file(package = "ascr"), "/TMB/")
  }
  
  if(file.exists(paste0(dir, "ascrTmb.o"))) unlink(paste0(dir, "ascrTmb.o"))
  if(file.exists(paste0(dir, "ascrTmb.dll"))) unlink(paste0(dir, 'ascrTmb.dll'))

  if(!debug.mode){
    TMB::compile(paste0(dir, "ascrTmb.cpp"), framework = "CppAD")
  } else {
    if(Sys.info()['sysname'] == 'Windows'){
      TMB::compile(paste0(dir, "ascrTmb.cpp"), "-O1 -g", DLLFLAGS="")
    } else {
      TMB::compile(paste0(dir, "ascrTmb.cpp"), "-O0 -g")
    }
    
  }

}
