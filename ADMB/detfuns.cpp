#include <admodel.h>
typedef dvariable (*detfn_pointer)(double, const dvar_vector&, double); // Mirrors what goes in detfn functions

// Half normal.
// Order of detpars: g0, sigma.
dvariable detfn_hn (double x, const dvar_vector &detpars, double cutoff)
{
  return detpars(1)*mfexp(-square(x)/(2*square(detpars(2))));
}

// Hazard rate.
// Order of detpars: g0, sigma, z.
dvariable detfn_hr (double x, const dvar_vector &detpars, double cutoff)
{
  return detpars(1)*(1 - mfexp(-pow(x/detpars(2),-detpars(3))));
}

// Threshold.
// Order of detpars: shape, scale.
dvariable detfn_th (double x, const dvar_vector &detpars, double cutoff)
{
  return 0.5 - 0.5*(2*cumd_norm((x/detpars(2) - detpars(1))*pow(2,0.5)) - 1);
}

// Log-link threshold.
// Order of detpars: shape1, shape2, scale.
dvariable detfn_logth (double x, const dvar_vector &detpars, double cutoff)
{
  return 0.5 - 0.5*(2*cumd_norm((detpars(1) - mfexp(detpars(2) + detpars(3)*x))*pow(2,0.5)) - 1);
}

// Signal strength.
// Order of detpars: b0ss, b1ss, sigmass.
dvariable detfn_ss (double x, const dvar_vector &detpars, double cutoff)
{
  return 1 - cumd_norm((cutoff - (detpars(1) - detpars(2)*x))/detpars(3));
}

// Log-link signal strength.
// Order of detpars: b0ss, b1ss, sigmass.
dvariable detfn_logss (double x, const dvar_vector &detpars, double cutoff)
{
  return 1 - cumd_norm((cutoff - mfexp(detpars(1) - detpars(2)*x))/detpars(3));
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
  }
  return(detfn) ;
}

typedef dvar_vector (*linkfn_pointer)(dvector, const dvar_vector&);

// Identity link.
dvar_vector linkfn_identity (dvector x, const dvar_vector &detpars)
{
  return detpars(1) - detpars(2)*x;
}

// Log link.
dvar_vector linkfn_log (dvector x, const dvar_vector &detpars)
{
  return mfexp(detpars(1) - detpars(2)*x);
}

// Set to zero.
dvar_vector linkfn_zero (dvector x, const dvar_vector &detpars)
{
  return 0;
}

linkfn_pointer get_linkfn(int linkfn_id)
{
  linkfn_pointer linkfn;
  switch(linkfn_id){
    case 1: linkfn = linkfn_identity; break;
    case 2: linkfn = linkfn_log; break;
    case 3: linkfn = linkfn_zero; break;
  }
}
