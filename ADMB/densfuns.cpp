#include <admodel.h>

// Normal distribution functions:

dvariable log_dnorm (double x, const prevariable mu, const prevariable sigma)
{
  const double pi = 3.141592653589793238463;
  return -0.5*log(2*pi) - log(sigma) - square(x - mu)/(2*square(sigma));
}

dvariable log_dpois (double x, const prevariable mu)
{
  return log_density_poisson(x, mu);
}
