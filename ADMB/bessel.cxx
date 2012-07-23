#include <admodel.h>

dvariable bessi0(dvariable x)
{
  // Returns the modified Bessel function of order zero for any real x
  // Numerical Recipes in C, 1992  
  dvariable ax,ans;
  dvariable y;
  if ((ax=fabs(x)) < 3.75) {
    y=x/3.75;
    y*=y;
    
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768-1+y*0.45813e-2)))));
  } else {
    y=3.75/ax;
    
    ans=(mfexp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.196281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
  }
  return ans;
}
