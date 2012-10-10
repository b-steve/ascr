## Funcitons for automatic generation of a .tpl file from admbsecr().
make.top.of.main <- function(memory){
  if (!is.null(memory)){
    cat("TOP_OF_MAIN_SECTION\n  arrmblsize=", memory, ";", file = "secr.tpl", sep = "") 
  }
}

make.common.variable.init <- function(){
  string <- "\n\nPROCEDURE_SECTION
  // Setting up variables
  int i,j;
  dvariable p,lambda,L1,L2,L3;
  dvar_matrix p1(1,ntraps,1,nmask);
  dvar_matrix p2(1,ntraps,1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);" 
  cat(string, file = "secr.tpl", sep = "", append = TRUE)
}

make.extra.variable.init <- function(methods){
  toastring <- "\n  dvariable nzz;
  dvar_vector toall(1,nmask);"["toa" %in% methods]
  angstring <- "\n  const double pi=3.14159265359;
  dvar_vector angll(1,nmask);"["ang" %in% methods]
  ssstring <- "\n  const double pi=3.14159265359;
  dvar_vector ssll(1,nmask);
  dvar_vector ess(1,nmask);"["ss" %in% methods]
  cat(toastring, angstring, ssstring, file = "secr.tpl", sep = "", append = TRUE)
}

make.probabilities <- function(){
  string <- "\n  // Probabilities of caputure at each location for each trap.
  // Add a small amount to prevent zeros.
  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;
  p2=1-p1;
  logp1=log(p1);
  logp2=log(p2);
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    for(j=1; j<=ntraps; j++){
      p*=p2(j)(i);
    }
    pm(i)=1-p;
  }"
  cat(string, file = "secr.tpl", append = TRUE)
}

make.L1.likelihood <- function(methods){
  common.start <- "\n  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
    wi1=capt(i)(1,ntraps);
    wi2=1-wi1;"
  toastring <-
    "\n    nzz=sum(wi1);
    toall=(1-nzz)*log(sigmatoa)-((row(toassq, i))/(2*square(sigmatoa)));"["toa" %in% methods]
  angstring <-
    "\n    angll=0;
    // Likelihood due to angles.
    for(j=1; j<=ntraps; j++){
      // Von-Mises density contribution for each trap.
      if(capt(i)(j)==1){
        angll+=kappa*cos(angcapt(i)(j)-row(ang,j));
      }
    }
    // Term in Von-Mises density not dependent on data.
    angll-=sum(wi1)*log(2*pi*bessi0(kappa));"["ang" %in% methods]
  ssstring <-
    "\n    ssll=0;
    // Likelihood due to signal strengths.
    for(j=1; j<=ntraps; j++){
      if (capt(i)(j)==1){
        ess=ssb0+ssb1*row(dist,j);
        ssll+=-log(sigmass)-(square(sscapt(i)(j)-ess)/(2*square(sigmass)));
      }
    }"["ss" %in% methods]
  common.end <- paste("\n    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)",
                      "+toall"["toa" %in% methods], "+angll"["ang" %in% methods],
                      "+ssll"["ss" %in% methods], "))+DBL_MIN);\n  }", sep = "")
  cat(common.start, toastring, angstring, ssstring, common.end, file = "secr.tpl",
      sep = "", append = TRUE)
}

make.together.likelihood <- function(){
  string <- "\n  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);"
  cat(string, file = "secr.tpl", append = TRUE)
}

make.globals <- function(methods){
  common.string <- "\n\nGLOBALS_SECTION
  #include <float.h>"
  ss.string <- "\n  #include <bessel.cxx>"["ang" %in% methods]
  cat(common.string, ss.string, "\n\n", file = "secr.tpl", sep = "", append = TRUE)
}

make.all.tpl <- function(memory, methods){
  if (file.exists("secr.tpl")){
    file.remove("secr.tpl")
  }
  file.create("secr.tpl")
  make.top.of.main(memory)
  make.common.variable.init()
  make.extra.variable.init(methods)
  make.probabilities()
  make.L1.likelihood(methods)
  make.together.likelihood()
  make.globals(methods)
}

