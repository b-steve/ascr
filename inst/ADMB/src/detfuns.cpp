#include <admodel.h>
typedef dvariable (*detfn_pointer)(double, const dvar_vector&, double, double);

// Half normal.
// Order of detpars: g0, sigma.
dvariable detfn_hn (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return detpars(1)*mfexp(-square(x)/(2*square(detpars(2))));
}

// Hazard rate.
// Order of detpars: g0, sigma, z.
dvariable detfn_hr (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return detpars(1)*(1 - mfexp(-pow(x/detpars(2), -detpars(3))));
}

// Threshold.
// Order of detpars: shape, scale.
dvariable detfn_th (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  dvariable z = (x/detpars(2) - detpars(1))*pow(2,0.5);
  return 0.5 - 0.5*(2*cumd_norm(z) - 1);
}

// Log-link threshold.
// Order of detpars: shape1, shape2, scale.
dvariable detfn_logth (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return 0.5 - 0.5*(2*cumd_norm((detpars(1) - mfexp(detpars(2) - detpars(3)*x))*pow(2,0.5)) - 1);
}

// Signal strength.
// Order of detpars: b0ss, b1ss, b2ss, sigmab0ss, sigmass.
dvariable detfn_ss (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return 1 - cumd_norm((cutoff - (detpars(1) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)/2))*x))/detpars(5));
}

// Log-link signal strength.
// Order of detpars: b0ss, b1ss, b2ss, sigmab0ss, sigmass.
dvariable detfn_logss (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return 1 - cumd_norm((cutoff - mfexp(detpars(1) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)/2))*x))/detpars(5));
}

// Spherical spreading signal strength.
// Order of detpars: b0ss, b1ss, b2ss, sigmab0ss, sigmass.
dvariable detfn_sphericalss (double x, const dvar_vector &detpars, double cutoff, double orientation)
{
  return 1 - cumd_norm((cutoff - (detpars(1) - 10*log10(square(x)) - (detpars(2) - (detpars(3)*(cos(orientation) - 1.0)/2))*(x - 1.0)))/detpars(5));
}

detfn_pointer get_detfn(int detfn_id)
{
  detfn_pointer detfn;
  switch(detfn_id){
  case 1: detfn = detfn_hn; break;
  case 2: detfn = detfn_hr; break;
  case 3: detfn = detfn_th; break;
  case 4: detfn = detfn_logth; break;
  case 5: detfn = detfn_ss; break;
  case 6: detfn = detfn_logss; break;
  case 7: detfn = detfn_sphericalss; break;
  }
  return(detfn) ;
}
