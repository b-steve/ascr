#include <admodel.h>

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

dvariable log_dmvn_diag (dvector x, const dvar_vector& mu, const prevariable& diag, const prevariable& offdiag, double dbl_min = 1e-150)
{
  double k = x.size();
  dvar_matrix diff(1,1,1,k);
  dvar_matrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  dvariable e = (diff*inv_diag(k, diag, offdiag)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det_diag(k, diag, offdiag) + dbl_min) + e);
}

dvariable log_dmvn_diag (dvector x, dvector mu, const prevariable& diag, const prevariable& offdiag, double dbl_min = 1e-150)
{
  double k = x.size();
  dmatrix diff(1,1,1,k);
  dmatrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  dvariable e = (diff*inv_diag(k, diag, offdiag)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det_diag(k, diag, offdiag) + dbl_min) + e);
}

double log_dmvn_diag (dvector x, dvector mu, double diag, double offdiag, double dbl_min = 1e-150)
{
  double k = x.size();
  dmatrix diff(1,1,1,k);
  dmatrix tdiff(1,k,1,1);
  for (int i = 1; i <= k; i++){
    diff(1,i) = x(i) - mu(i);
    tdiff(i,1) = diff(1,i);
  }
  double e = (diff*inv_diag(k, diag, offdiag)*tdiff)(1,1);
  return -0.5*(k*log(2*M_PI) + log(det_diag(k, diag, offdiag) + dbl_min) + e);
}

// Multivariate normal cumulative distribution functions. 

// Uses either Gauss-Hermite quadrature (gh = true) or the rectangle
// rule (gh = false). Implements the simplification of the MVN CDF
// from Dunnett and Sobel (1955), see Kotz, Balakrishnan and Johnson
// (2000) pp. 134.

// Vector x must be normalised, i.e., mu = 0, diagonal elements of
// sigma are all 1.
dvariable pmvn (const dvar_vector& x, const dvariable& corr, bool gh, dvector weights, dvector nodes, double n_quadpoints, double lower, double upper)
{
  int k = x.size();
  int i, j;
  dvariable out = 0;
  dvariable summand;
  if (gh){
    for (i = 1; i <= n_quadpoints; i++){
      summand = weights(i)/pow(M_PI, 0.5);
      for (j = 1; j <= k; j++){
        // Note nodes altered due to change of variable technique.
	summand *= cumd_norm((x(j) - pow(corr, 0.5)*nodes(i)*pow(2, 0.5))/pow(1 - corr, 0.5));
      }
      out += summand;
    }
  } else {
    double bin_width = (upper - lower)/n_quadpoints;
    double bin_mid;
    for (i = 1; i <= n_quadpoints; i++){
      bin_mid = lower + bin_width*(i - 0.5);
      summand = mfexp(log_dnorm(bin_mid, 0, 1));
      for (j = 1; j <= k; j++){
	summand *= cumd_norm((x(j) - pow(corr, 0.5)*bin_mid)/pow(1 - corr, 0.5));
      }
      summand *= bin_width;
      out += summand;
    }
  }
  return out;
}

double pmvn (dvector x, double corr, bool gh, dvector weights, dvector nodes, double n_quadpoints, double lower, double upper)
{
  int k = x.size();
  int i, j;
  double out = 0;
  double summand;
  if (gh){
    cout << "gh" << endl;
    // Altering nodes due to change of variable technique.
    for (i = 1; i <= n_quadpoints; i++){
      summand = weights(i)/pow(M_PI, 0.5);
      for (j = 1; j <= k; j++){
	summand *= cumd_norm((x(j) - pow(corr, 0.5)*nodes(i)*pow(2, 0.5))/pow(1 - corr, 0.5));
      }
      out += summand;
    }
  } else {
    cout << "rect" << endl;
    double bin_width = (upper - lower)/n_quadpoints;
    double bin_mid;
    for (i = 1; i <= n_quadpoints; i++){
      bin_mid = lower + bin_width*(i - 0.5);
      summand = mfexp(log_dnorm(bin_mid, 0, 1));
      for (j = 1; j <= k; j++){
	summand *= cumd_norm((x(j) - pow(corr, 0.5)*bin_mid)/pow(1 - corr, 0.5));
      }
      summand *= bin_width;
      out += summand;
    }
  }
  return out;
}

