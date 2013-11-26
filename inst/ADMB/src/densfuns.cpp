#include <admodel.h>

// Helper prototypes:
dvariable bessi0 (dvariable x);

// Gamma distribution functions:

dvariable log_dgamma (double x, const prevariable& alpha, const prevariable& beta)
{
  return alpha*log(beta) + (alpha - 1)*log(x) - (beta*x) - gammln(alpha);
}

dvar_vector log_dgamma (double x, const prevariable& alpha, const dvar_vector& beta)
{
  return alpha*log(beta) + (alpha - 1)*log(x) - (beta*x) - gammln(alpha);
}


// Normal distribution functions:

dvariable log_dnorm (double x, const prevariable& mu, const prevariable& sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - square(x - mu)/(2*square(sigma));
}

dvar_vector log_dnorm (double x, const dvar_vector& mu, const prevariable& sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - square(x - mu)/(2*square(sigma));
}

// Poisson distribution functions:

dvariable log_dpois (double x, const prevariable& mu)
{
  return log_density_poisson(x, mu);
}

// Von-mises distribution functions:

dvariable log_dvm (double x, double mu, const prevariable& kappa)
{
  return kappa*cos(x - mu) - log(2*M_PI*bessi0(kappa));
}

dvar_vector log_dvm (double x, dvector mu, const prevariable& kappa)
{
  return kappa*cos(x - mu) - log(2*M_PI*bessi0(kappa));
}

//Helpers:
dvariable bessi0 (dvariable x)
{
  dvariable ax,ans;
  dvariable y;
    if ((ax = fabs(x)) < 3.75){
    y = x/3.75;
    y *= y;
    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 + y*(0.2659732 + y*(0.360768e-1 + y*0.45813e-2)))));
  } else {
    y = 3.75/ax;
    ans = (exp(ax)/sqrt(ax))*(0.39894228 + y*(0.1328592e-1 + y*(0.225319e-2 + y*(-0.157565e-2 + y*(0.916281e-2 + y*(-0.2057706e-1 + y*(0.2635537e-1 + y*(-0.1647633e-1 + y*0.392377e-2))))))));
  }
  return ans;
}
