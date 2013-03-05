## Functions for automatic generation of a .tpl file from admbsecr().

make.all.tpl.easy <- function(memory, method, detfn, parnames){
  make.top.of.main(memory)
  simpletext <- "\nPROCEDURE_SECTION\n  // Setting up variables\n  int i, j;\n  dvariable p, d;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  // Probability of detection at any trap for each location.\n  for (i = 1; i <= nmask; i++){\n    p = 1;\n    for (j = 1; j <= ntraps; j++){\n      d = dist(j,i);\n      // Flag for detection function insertion.\n      p1(j,i) = //@DETFN;\n      p2(j,i) = 1 - p1(j,i);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp1 = log(p1 + DBL_MIN);\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1 = 0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    wi1 = row(capt,i);\n    L1 += log(D*sum(mfexp(wi1*logp1 + (1 - wi1)*logp2)) + DBL_MIN);\n  }\n\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1 + L2 + L3);\n  if (trace == 1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n" 
  toatext <- "\nPROCEDURE_SECTION\n  // Setting up variables\n  int i, j;\n  dvariable p, d, nzz;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector toall(1,nmask);\n  // Probability of detection at any trap for each location.\n  for (i = 1; i <= nmask; i++){\n    p = 1;\n    for (j = 1; j <= ntraps; j++){\n      d = dist(j,i);\n      // Flag for detection function insertion.\n      p1(j,i) = //@DETFN;\n      p2(j,i) = 1 - p1(j,i);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp1 = log(p1 + DBL_MIN);\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1=0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    wi1 = row(capt,i);\n    nzz = sum(wi1);\n    toall = (1-nzz)*log(sigmatoa) - ((row(toassq,i))/(2*square(sigmatoa)));\n    L1 += log(sum(mfexp(log(D) + (wi1*logp1 + (1-wi1)*logp2) + toall)) + DBL_MIN);\n  }\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1 + L2 + L3);\n  if (trace == 1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n" 
  angtext <- "\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i, j;\n  dvariable p, d;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector angll(1,nmask);\n  // Probability of detection at any trap for each location.\n  for (i = 1; i <= nmask; i++){\n    p = 1;\n    for (j = 1; j <= ntraps; j++){\n      d = dist(j,i);\n      // Flag for detection function insertion.\n      p1(j,i) = //@DETFN;\n      p2(j,i) = 1 - p1(j,i);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp1 = log(p1 + DBL_MIN);\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1 = 0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    wi1 = row(capt,i);\n    angll = 0;\n    // Likelihood due to angles.\n    for (j = 1; j <= ntraps; j++){\n      // Von-Mises density contribution for each trap.\n      if (capt(i,j) == 1){\n        angll += kappa*cos(angcapt(i,j)-row(ang,j));\n      }\n    }\n    // Term in Von-Mises density not dependent on data.\n    angll -= sum(wi1)*log(2*pi*bessi0(kappa));\n    L1 += log(sum(mfexp(log(D) + (wi1*logp1 + (1 - wi1)*logp2) + angll)) + DBL_MIN);\n  }\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1 + L2 + L3);\n  if (trace == 1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n  #include <bessel.cxx>\n\nREPORT_SECTION\n" 
  sstext <- "\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i, j;\n  dvariable p, d;\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector ci1(1,ntraps);\n  dvar_matrix muss(1,ntraps,1,nmask);\n\n  muss = //@LINKFN(ssb0 + ssb1*dist);\n  for (i = 1; i <= nmask; i++){\n    p=1;\n    for (j = 1; j <= ntraps; j++){\n      p2(j,i) = cumd_norm((c - muss(j,i))/sigmass);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1=0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    logp1 = 0;\n    wi1 = row(sscapt,i);\n    ci1 = row(capt,i);\n    for(j = 1; j <= ntraps; j++){\n      if (ci1(j) == 1){\n        logp1(j)(1,nmask) = -log(sigmass) - log(sqrt(2*pi)) + (square(wi1(j)-row(muss,j))/(-2*square(sigmass)));\n      }\n    }\n    L1 += log(D*sum(mfexp(ci1*logp1 + (1-ci1)*logp2) + DBL_MIN));\n  }\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1+L2+L3);\n  if (trace==1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n" 
  sstoatext <- "\nPROCEDURE_SECTION\n  // Setting up variables\n  const double pi=3.14159265359;\n  int i, j;\n  dvariable p, d, nzz;\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector ci1(1,ntraps);\n  dvar_vector toall(1,nmask);\n  dvar_matrix muss(1,ntraps,1,nmask);\n\n  muss = //@LINKFN(ssb0 + ssb1*dist);\n  for (i = 1; i <= nmask; i++){\n    p=1;\n    for (j = 1; j <= ntraps; j++){\n      p2(j,i) = cumd_norm((c - muss(j,i))/sigmass);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1=0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    logp1 = 0;\n    wi1 = row(sscapt,i);\n    ci1 = row(capt,i);\n    nzz = sum(ci1);\n    for (j = 1; j <= ntraps; j++){\n      if (ci1(j) == 1){\n        logp1(j)(1,nmask)= -log(sigmass) - log(sqrt(2*pi)) + (square(wi1(j)-row(muss,j))/(-2*square(sigmass)));\n      }\n    }\n    toall = (1 - nzz)*log(sigmatoa) - ((row(toassq,i))/(2*square(sigmatoa)));\n    L1 += log(D*sum(mfexp((ci1*logp1 + (1 - ci1)*logp2) + toall) + DBL_MIN));\n  }\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1+L2+L3);\n  if (trace==1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n" 
  disttext <- "\nPROCEDURE_SECTION\n  // Setting up variables.\n  const double pi=3.14159265359;\n  int i, j;\n  dvariable p, d;\n  dvar_matrix p1(1,ntraps,1,nmask);\n  dvar_matrix p2(1,ntraps,1,nmask);\n  dvar_matrix logp1(1,ntraps,1,nmask);\n  dvar_matrix logp2(1,ntraps,1,nmask);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector distll(1,nmask);\n  dvar_vector beta(1,nmask);\n  // Probability of detection at any trap for each location.\n  for (i = 1; i <= nmask; i++){\n    p = 1;\n    for (j = 1; j <= ntraps; j++){\n      d = dist(j,i);\n      // Flag for detection function insertion.\n      p1(j,i) = //@DETFN;\n      p2(j,i) = 1 - p1(j,i);\n      p *= p2(j,i);\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  logp1 = log(p1 + DBL_MIN);\n  logp2 = log(p2 + DBL_MIN);\n  dvariable L1 = 0;\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    wi1 = row(capt,i);\n    distll = 0;\n    // Likelihood due to distances.\n    for(j = 1; j <= ntraps; j++){\n      // Gamma density contribution for each trap.\n      if(capt(i,j)==1){\n\tbeta = alpha/row(dist,j);\n\tdistll += alpha*log(beta) + (alpha - 1)*log(distcapt(i)(j)) - (beta*distcapt(i)(j)) - gammln(alpha);\n      }\n    }\n    L1 += log(sum(mfexp(log(D) + (wi1*logp1 + (1 - wi1)*logp2) + distll)) + DBL_MIN);\n  }\n  // Putting log-likelihood together.\n  dvariable lambda = A*D*sum(pm);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1 + L2 + L3);\n  if (trace == 1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n"
  mrdstext <- "\nPROCEDURE_SECTION\n  // Setting up variables.\n  int i,j;\n  dvariable p, p1, d;\n  dvar_matrix indivp1(1,n,1,ntraps);\n  dvar_matrix indivp2(1,n,1,ntraps);\n  dvar_matrix logindivp1(1,n,1,ntraps);\n  dvar_matrix logindivp2(1,n,1,ntraps);\n  dvar_matrix logprobs(1,n,1,ntraps);\n  dvar_vector pm(1,nmask);\n  dvar_vector wi1(1,ntraps);\n  dvar_vector wi2(1,ntraps);\n  // Probability of detection at any trap for each location.\n  for (i = 1; i <= nmask; i++){\n    p = 1;\n    for(j = 1; j <= ntraps; j++){\n      d = dist(j,i);\n      // Flag for detection function insertion.\n      p1 = //@DETFN;\n      p *= 1 - p1;\n    }\n    pm(i) = 1 - p + DBL_MIN;\n  }\n  // Probability of capture histories for each animal.\n  for (i = 1; i <= n; i++){\n    for (j = 1; j <= ntraps; j++){\n      d = indivdist(i,j);\n      indivp1(i,j) = //@DETFN;\n      indivp2(i,j) = 1 - indivp1(i,j);\n    }\n  }\n  logindivp1 = log(indivp1 + DBL_MIN);\n  logindivp2 = log(indivp2 + DBL_MIN);\n  logprobs = elem_prod(logindivp1,capt) + elem_prod(logindivp2,1 - capt);\n  dvariable L1 = sum(logprobs) + n*log(D);\n  dvariable L2 = -n*log(D*sum(pm));\n  dvariable lambda = A*D*sum(pm);\n  dvariable L3 = log_density_poisson(n,lambda);\n  f = -(L1+L2+L3);\n  if (trace == 1){\n    //@TRACE;\n  }\n\nGLOBALS_SECTION\n  #include <float.h>\n\nREPORT_SECTION\n\n" 
  alltext <- c(simpletext, toatext, angtext, sstext, sstoatext, disttext, mrdstext)
  names(alltext) <- c("simple", "toa", "ang", "ss", "sstoa", "dist", "mrds")
  codestring <- alltext[method]
  if (!(method == "ss" | method == "sstoa")){
    codesplit <- strsplit(codestring, "//@DETFN")[[1]]
    if (detfn == "hn"){
      dfstr <- "g0*mfexp(-square(d)/(2*square(sigma)))"
    } else if (detfn == "th"){
      dfstr <- "0.5 - 0.5*(2*cumd_norm((shape - scale*d)*pow(2,0.5)) - 1)"
    } else if (detfn == "logth"){
      dfstr <- "0.5 - 0.5*(2*cumd_norm((shape1 - mfexp(shape2 + scale*d))*pow(2,0.5)) - 1)"
    } else if (detfn == "hr"){
      dfstr <- "g0*(1 - mfexp(-pow(d/sigma,-z)))"
    } else {
      stop("detection function not supported.")
    }
  } else {
    codesplit <- strsplit(codestring, "//@LINKFN")[[1]]
    if (detfn == "identity"){
      dfstr <- ""
    } else if (detfn == "log"){
      dfstr <- "mfexp"
    } else {
      stop("detection function not supported.")
    }
  }
  clen <- length(codesplit)
  codestring <- splitinsert(codesplit, dfstr)
  codesplit <- strsplit(codestring, "//@TRACE")[[1]]
  trstr <- paste("cout << \"", parnames[1], ": \" << ", parnames[1], " <<",
                 paste(" \", ", parnames[-1], ": \" << ", parnames[-1], " <<",
                       sep = "", collapse = ""), " \", loglik: \" << -f << endl",
                 sep = "")
  codestring <- splitinsert(codesplit, trstr)
  cat(codestring, file = "secr.tpl", append = !is.null(memory))
}

make.top.of.main <- function(memory){
  if (!is.null(memory)){
    cat("TOP_OF_MAIN_SECTION\n  arrmblsize=", memory, ";", file = "secr.tpl", sep = "")
  }
}

splitinsert <- function(codesplit, str){
  clen <- length(codesplit)
  codestring <- codesplit[1]
  if (clen > 1){
    for (i in 2:clen){
      codestring <- paste(codestring, str, codesplit[i], sep = "")
    }
  }
  codestring
}

make.bessel <- function(){
  besseltext <- "#include <admodel.h>\n\ndvariable bessi0(dvariable x)\n{\n  dvariable ax,ans;\n  dvariable y;\n  \n  if ((ax=fabs(x)) < 3.75) {\n    y=x/3.75;\n    y*=y;\n    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492\n\t\t\t\t\t +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));\n  } else {\n    y=3.75/ax;\n    ans=(exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1\n\t\t\t\t\t  +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2\n\t\t\t\t\t\t\t\t\t     +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t  +y*0.392377e-2))))))));\n  }\n  return ans;\n}\n"
  cat(besseltext, file = "bessel.cxx")
}

get.tpl.text <- function(methods = c("simple", "toa", "ang", "ss", "sstoa", "dist", "mrds"),
                         admbwd = "~/admbsecr/ADMB/"){
  t <- character(length(methods))
  names(t) <- methods
  for (i in methods){
    file <- paste(admbwd, i, "secr.tpl", sep = "")
    t[i] <- readChar(file, file.info(file)$size)
  }
  t
}
