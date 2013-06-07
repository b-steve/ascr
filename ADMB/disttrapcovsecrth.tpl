
PROCEDURE_SECTION
  // Setting up variables
  const double DBL_MIN = 1e-150;
  int i,j;
  dvariable d,p,lambda,L1,L2,L3;
  dvar_vector p11(1,nmask);
  dvar_vector p12(1,nmask);
  dvar_vector p21(1,nmask);
  dvar_vector p22(1,nmask);
  dvar_matrix logp1(1,ntraps,1,nmask);
  dvar_matrix logp2(1,ntraps,1,nmask);
  dvar_vector pm(1,nmask);
  dvar_vector wi1(1,ntraps);
  dvar_vector wi2(1,ntraps);
  dvar_vector distll(1,nmask);
  dvar_vector beta(1,nmask);
  // Probabilities of caputure at each location for each trap.
  // Add a small amount to prevent zeros.
  for(i=1; i<=nmask; i++){
    d = dist(1,i);
    p11(i) = 0.5 - 0.5*(2*cumd_norm((shape1 - scale1*d)*pow(2,0.5)) - 1) + DBL_MIN;
    p21(i) = 0.5 - 0.5*(2*cumd_norm((shape2 - scale2*d)*pow(2,0.5)) - 1) + DBL_MIN;
    p12(i) = 1 - p11(i);
    p22(i) = 1 - p21(i);
    logp1(1,i)=log(p11(i));
    logp1(2,i)=log(p21(i));
    logp2(1,i)=log(p12(i));
    logp2(2,i)=log(p22(i));
  }
  // Probability of detection at any trap for each location.
  for(i=1; i<=nmask; i++){
    p=1;
    p*=p12(i);
    p*=p22(i);
    pm(i)=1-p;
  }
  L1=0;
  // Probability of capture histories for each animal.
  for(i=1; i<=n; i++){
    wi1=capt(i)(1,ntraps);
    wi2=1-wi1;
    distll=0;
    // Likelihood due to distances.
    for(j=1; j<=ntraps; j++){
      // Gamma density contribution for each trap.
      if(capt(i)(j)==1){
	beta=alpha/row(dist,j);
	distll+=alpha*log(beta)+(alpha-1)*log(distcapt(i)(j))-(beta*distcapt(i)(j))-gammln(alpha);
      }
    }
    L1+=log(sum(mfexp(log(D)+(wi1*logp1+wi2*logp2)+distll))+DBL_MIN);
  }
  // Putting log-likelihood together.
  lambda=A*D*sum(pm);
  L2=-n*log(D*sum(pm));
  L3=log_density_poisson(n,lambda);
  f=-(L1+L2+L3);
  if (trace == 1){
  cout << "D: " << D << ", shape1: " << shape1 << ", scale1: " << scale1 << ", shape2: " << shape2 << ", scale2: " << scale2 << ", alpha: " << alpha << ", loglik: " << -f << endl;
  }



