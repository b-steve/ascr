#include <admodel.h>

// Helper prototypes:
dvariable bessi0 (dvariable x);

// Gamma distribution density functions:

dvariable log_dgamma (double x, const prevariable& alpha, const prevariable& beta)
{
  return alpha*log(beta) + (alpha - 1)*log(x) - (beta*x) - gammln(alpha);
}

dvar_vector log_dgamma (double x, const prevariable& alpha, const dvar_vector& beta)
{
  return alpha*log(beta) + (alpha - 1)*log(x) - (beta*x) - gammln(alpha);
}


// Normal distribution density functions:

dvariable log_dnorm (double x, const prevariable& mu, const prevariable& sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - square(x - mu)/(2*square(sigma));
}

dvar_vector log_dnorm (double x, const dvar_vector& mu, const prevariable& sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - square(x - mu)/(2*square(sigma));
}

double log_dnorm (double x, double mu, const double sigma)
{
  return -0.5*log(2*M_PI) - log(sigma) - square(x - mu)/(2*square(sigma));
}


// Poisson distribution density functions:

dvariable log_dpois (double x, const prevariable& mu)
{
  return log_density_poisson(x, mu);
}

// Von-mises distribution density functions:

dvariable log_dvm (double x, double mu, const prevariable& kappa)
{
  return kappa*cos(x - mu) - log(2*M_PI*bessi0(kappa));
}

dvar_vector log_dvm (double x, dvector mu, const prevariable& kappa)
{
  return kappa*cos(x - mu) - log(2*M_PI*bessi0(kappa));
}

// Multivariate normal distribution density functions:

dvariable log_dmvn (dvector x, const dvar_vector& mu, const dvar_matrix& sigma)
{
  double k = x.size();
  dvar_matrix diff(1,1,1,k);
  dvar_matrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  dvariable e = (diff*inv(sigma)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det(sigma)) + e);
}

dvariable log_dmvn (dvector x, dvector mu, const dvar_matrix& sigma)
{
  double k = x.size();
  dmatrix diff(1,1,1,k);
  dmatrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  dvariable e = (diff*inv(sigma)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det(sigma)) + e);
}

double log_dmvn (dvector x, dvector mu, dmatrix sigma)
{
  double k = x.size();
  dmatrix diff(1,1,1,k);
  dmatrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  double e = (diff*inv(sigma)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det(sigma)) + e);
}

// Multivariate normal cumulative distribution functions:

// Must be normalised, i.e., mu = 0, diagonal elements of sigma are all 1.
dvariable log_pmvn (const dvar_vector& x, const prevariable& corr, double lower = -5, double upper = 5, double n_quadpoints = 10)
{
  int k = x.size();
  int i, j;
  double bin_width = (upper - lower)/n_quadpoints;
  double bin_mid;
  dvariable out = 0;
  dvariable summand;
  dvariable z;
  for (i = 1; i <= n_quadpoints; i++){
    bin_mid = lower + bin_width*(i - 0.5);
    summand = mfexp(log_dnorm(bin_mid, 0, 1));
    for (j = 1; j <= k; j++){
      summand *= cumd_norm((x(j) - pow(corr, 0.5)*bin_mid)/pow(1 - corr, 0.5));
    }
    summand *= bin_width;
    out += summand;
  }
  return out;
}

double log_pmvn (dvector x, double corr, double lower = -5, double upper = 5, double n_quadpoints = 10)
{
  int k = x.size();
  int i, j;
  double bin_width = (upper - lower)/n_quadpoints;
  double bin_mid;
  double out = 0;
  double summand;
  double z;
  for (i = 1; i <= n_quadpoints; i++){
    bin_mid = lower + bin_width*(i - 0.5);
    summand = mfexp(log_dnorm(bin_mid, 0, 1));
    for (j = 1; j <= k; j++){
      summand *= cumd_norm((x(j) - pow(corr, 0.5)*bin_mid)/pow(1 - corr, 0.5));
    }
    summand *= bin_width;
    out += summand;
  }
  return out;
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
