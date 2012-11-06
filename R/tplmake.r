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


## Above functions will allow you to combine methods. For now, let's just import the
## .tpl files as they are in the ADMB folder.
make.all.tpl.easy <- function(memory, methods){
  make.top.of.main(memory)
  simpletext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  int i,j;\n  dvariable p,lambda,L1,L2,L3;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector wi2(1,ntraps);\n  // Probabilities of caputure at each location for each trap.\n  // Add a small amount to prevent zeros.\n  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;\n  p2=1-p1;\n  logp1=log(p1);\n  logp2=log(p2);\n  // Probability of detection at any trap for each location.\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p*=p2(j)(i);\n    }\n    pm(i)=1-p;\n  }\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    wi1=row(capt,i);\n    wi2=1-wi1;\n    L1+=log(D*sum(mfexp(wi1*logp1+wi2*logp2)));\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm);\n  L2=-n*log(D*sum(pm));\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace == 1){\n    cout << \"D: \" << D << \", g0: \" << g0 << \", sigma: \" << sigma << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n"
  toatext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  int i,j;\n  dvariable p,lambda,L1,L2,L3,nzz;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector wi2(1,ntraps);\n  dvar_vector toall(1,nmask);\n  // Probabilities of caputure at each location for each trap.\n  // Add a small amount to prevent zeros.\n  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;\n  p2=1-p1;\n  logp1=log(p1);\n  logp2=log(p2);\n  // Probability of detection at any trap for each location.\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p*=p2(j)(i);\n    }\n    pm(i)=1-p;\n  }\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    wi1=row(capt,i);\n    wi2=1-wi1;\n    nzz=sum(wi1);\n    toall=(1-nzz)*log(sigmatoa)-((row(toassq, i))/(2*square(sigmatoa)));\n    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+toall)));\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm);\n  L2=-n*log(D*sum(pm));\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace == 1){\n    cout << \"D: \" << D << \", g0: \" << g0 << \", sigma: \" << sigma << \", sigmatoa: \" << sigmatoa << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n"
  angtext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i,j;\n  dvariable p,lambda,L1,L2,L3;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector wi2(1,ntraps);\n  dvar_vector angll(1,nmask);\n  // Probabilities of caputure at each location for each trap.\n  // Add a small amount to prevent zeros.\n  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;\n  p2=1-p1;\n  logp1=log(p1);\n  logp2=log(p2);\n  // Probability of detection at any trap for each location.\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p*=p2(j)(i);\n    }\n    pm(i)=1-p;\n  }\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    wi1=capt(i)(1,ntraps);\n    wi2=1-wi1;\n    angll=0;\n    // Likelihood due to angles.\n    for(j=1; j<=ntraps; j++){\n      // Von-Mises density contribution for each trap.\n      if(capt(i)(j)==1){\n        angll+=kappa*cos(angcapt(i)(j)-row(ang,j));\n      }\n    }\n    // Term in Von-Mises density not dependent on data.\n    angll-=sum(wi1)*log(2*pi*bessi0(kappa));\n    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+angll))+DBL_MIN);\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm);\n  L2=-n*log(D*sum(pm));\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace == 1){\n  cout << \"D: \" << D << \", g0: \" << g0 << \", sigma: \" << sigma << \", kappa: \" << kappa << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n  #include <bessel.cxx>\n\nREPORT_SECTION\n"
  sstext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  const double dmin=DBL_MIN;\n  int i,j,k;\n  dvariable p,lambda,L1,L2,L3;\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector ci1(1,ntraps);\n  dvar_vector ssll(1,nmask);\n  dvar_matrix muss(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  muss=mfexp(ssb0+ssb1*dist);\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p2(j,i)=cumd_norm((c-muss(j,i))/sigmass);\n      p*=p2(j,i);\n    }\n    pm(i)=1-p;\n  }\n  logp2=log(p2+dmin);\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    logp1=0;\n    ssll=0;\n    wi1=row(sscapt,i);\n    ci1=row(capt,i);\n    for(j=1; j<=ntraps; j++){\n      if (ci1(j)==1){\n        logp1(j)(1,nmask)=-log(sigmass)-log(sqrt(2*pi))+(square(wi1(j)-row(muss,j))/(-2*square(sigmass)));\n      }\n    }\n    L1+=log(D*sum(mfexp(ci1*logp1+(1-ci1)*logp2)+dmin));\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm)+dmin;\n  L2=-n*log(D*sum(pm)+dmin);\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace==1){\n    cout << \"D: \" << D << \", ssb0: \" << ssb0 << \", ssb1: \" << ssb1 << \", sigmass: \" << sigmass << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n"
  sstoatext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i,j,k;\n  dvariable p,lambda,L1,L2,L3,nzz;\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector ci1(1,ntraps);\n  dvar_vector ssll(1,nmask);\n  dvar_vector toall(1,nmask);\n  dvar_matrix muss(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  muss=mfexp(ssb0+ssb1*dist);\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p2(j,i)=cumd_norm((c-muss(j,i))/sigmass);\n      p*=p2(j,i);\n    }\n    pm(i)=1-p;\n  }\n  logp2=log(p2+DBL_MIN);\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    logp1=0;\n    ssll=0;\n    wi1=row(sscapt,i);\n    ci1=row(capt,i);\n    nzz=sum(ci1);\n    for(j=1; j<=ntraps; j++){\n      if (ci1(j)==1){\n        logp1(j)(1,nmask)=-log(sigmass)-log(sqrt(2*pi))+(square(wi1(j)-row(muss,j))/(-2*square(sigmass)));\n      }\n    }\n    toall=(1-nzz)*log(sigmatoa)-((row(toassq, i))/(2*square(sigmatoa)));\n    L1+=log(D*sum(mfexp((ci1*logp1+(1-ci1)*logp2)+toall)+DBL_MIN));\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm)+DBL_MIN;\n  L2=-n*log(D*sum(pm)+DBL_MIN);\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace == 1){\n    cout << \"D: \" << D << \", sigmatoa: \" << sigmatoa << \", ssb0: \" << ssb0 << \", ssb1: \" << ssb1 << \", sigmass: \" << sigmass << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n"
  disttext <- "\n\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i,j;\n  dvariable p,lambda,L1,L2,L3;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector wi2(1,ntraps);\n  dvar_vector distll(1,nmask);\n  dvar_vector beta(1,nmask);\n  // Probabilities of caputure at each location for each trap.\n  // Add a small amount to prevent zeros.\n  p1=g0*mfexp(-square(dist)/(2*square(sigma)))+DBL_MIN;\n  p2=1-p1;\n  logp1=log(p1);\n  logp2=log(p2);\n  // Probability of detection at any trap for each location.\n  for(i=1; i<=nmask; i++){\n    p=1;\n    for(j=1; j<=ntraps; j++){\n      p*=p2(j)(i);\n    }\n    pm(i)=1-p;\n  }\n  L1=0;\n  // Probability of capture histories for each animal.\n  for(i=1; i<=n; i++){\n    wi1=capt(i)(1,ntraps);\n    wi2=1-wi1;\n    distll=0;\n    // Likelihood due to distances.\n    for(j=1; j<=ntraps; j++){\n      // Gamma density contribution for each trap.\n      if(capt(i)(j)==1){\n\tbeta=alpha/row(dist,j);\n\tdistll+=alpha*log(beta)+(alpha-1)*log(distcapt(i)(j))-(beta*distcapt(i)(j))-gammln(alpha);\n      }\n    }\n    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+distll))+DBL_MIN);\n  }\n  // Putting log-likelihood together.\n  lambda=A*D*sum(pm);\n  L2=-n*log(D*sum(pm));\n  L3=log_density_poisson(n,lambda);\n  f=-(L1+L2+L3);\n  if (trace == 1){\n  cout << \"D: \" << D << \", g0: \" << g0 << \", sigma: \" << sigma << \", alpha: \" << alpha << \", loglik: \" << -f << endl;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n  #include <bessel.cxx>\n\nREPORT_SECTION\n"
  alltext <- c(simpletext, toatext, angtext, sstext, sstoatext, disttext)
  names(alltext) <- c("simple", "toa", "ang", "ss", "sstoa", "dist")
  cat(alltext[methods], file = "secr.tpl", append = !is.null(memory))
}

get.tpl.text <- function(methods = c("simple", "toa", "ang", "ss", "sstoa", "dist"),
                         admbwd = "/home/ben/admbsecr/ADMB/"){
  t <- character(length(methods))
  names(t) <- methods
  for (i in methods){
    file <- paste(admbwd, i, "secr.tpl", sep = "")
    t[i] <- readChar(file, file.info(file)$size)
  }
  t
}
