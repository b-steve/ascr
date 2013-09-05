#include <admodel.h>
typedef dvariable (*detfn_pointer)(double, const dvar_vector&); // Mirrors what goes in detfn functions

// Half normal.
// Order of detpars: g0, sigma.
dvariable detfn_hn (double x, const dvar_vector &detpars)
{
  return detpars[1]*mfexp(-square(x)/(2*square(detpars[2])));
}

// Hazard rate.
// Order of detpars: g0, sigma, z.
dvariable detfn_hr (double x, const dvar_vector &detpars)
{
  return detpars[1]*(1 - mfexp(-pow(x/detpars[2],-detpars[3])));
}

// Log-link threshold.
// Order of detpars: shape1, shape2, scale.
dvariable detfn_logth (double x, const dvar_vector &detpars)
{
  return 0.5 - 0.5*(2*cumd_norm((detpars[1] - mfexp(detpars[2] + detpars[3]*x))*pow(2,0.5)) - 1);
}

// Threshold.
// Order of detpars: shape, scale.
dvariable detfn_th (double x, const dvar_vector &detpars)
{
  return 0.5 - 0.5*(2*cumd_norm((x/detpars[2] - detpars[1])*pow(2,0.5)) - 1);
}

detfn_pointer get_detfn(int detfn_id)
{     
  detfn_pointer detfn;    
  switch(detfn_id){
    case 1: detfn = detfn_hn; break;
    case 2: detfn = detfn_hr; break;
    case 3: detfn = detfn_th; break;
    case 4: detfn = detfn_logth; break;
  }
  return(detfn) ;
}
