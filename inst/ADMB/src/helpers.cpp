#include <admodel.h>

// Bessel function.
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

// Filling a diagonal matrix.
dvar_matrix fill_diag (int dim, const prevariable& diag, const prevariable& offdiag)
{
  int i, j;
  dvar_matrix out(1, dim, 1, dim);
  for (i = 1; i <= dim; i++){
    for (j = i; j <= dim; j++){
      if (i == j){
        out(i, j) = diag;
      } else {
        out(i, j) = offdiag;
        out(j, i) = offdiag;
      }
    }
  }
  return out;
}

dmatrix fill_diag (int dim, double diag, double offdiag)
{
  int i, j;
  dmatrix out(1, dim, 1, dim);
  for (i = 1; i <= dim; i++){
    for (j = i; j <= dim; j++){
      if (i == j){
        out(i, j) = diag;
      } else {
        out(i, j) = offdiag;
        out(j, i) = offdiag;
      }
    }
  }
  return out;
}


// Inverting a matrix where all diagonal elements are equal to diag,
// and all off-diagonal elements are equal to offdiag.
dvar_matrix inv_diag (int dim, const prevariable& diag, const prevariable& offdiag)
{
  double k = dim;
  dvar_matrix out(1, k, 1, k);
  if (k == 1){
    out(1, 1) = 1/diag;
  } else {
    dvariable mult = 1/(square(diag) + (k - 2)*diag*offdiag - (k - 1)*square(offdiag));
    dvariable a = diag + (k - 2)*offdiag;
    dvariable b = -offdiag;
    out = mult*fill_diag(k, a, b);
  }
  return out;
}

dmatrix inv_diag (int dim, double diag, double offdiag)
{
  double k = dim;
  dmatrix out(1, k, 1, k);
  if (k == 1){
    out(1, 1) = 1/diag;
  } else {
    double mult = 1/(square(diag) + (k - 2)*diag*offdiag - (k - 1)*square(offdiag));
    double a = diag + (k - 2)*offdiag;
    double b = -offdiag;
    out = mult*fill_diag(k, a, b);
  }
  return out;
}

// Discriminant of a matrix where all diagonal elements are equal to
// diag, and all off-diagonal elements are equal to offdiag.
dvariable det_diag (int dim, const prevariable& diag, const prevariable& offdiag)
{
  double k = dim;
  return pow(diag - offdiag, k - 1)*(diag + (k - 1)*offdiag);
}

double det_diag (int dim, double diag, double offdiag)
{
  double k = dim;
  return pow(diag - offdiag, k - 1)*(diag + (k - 1)*offdiag);
}
