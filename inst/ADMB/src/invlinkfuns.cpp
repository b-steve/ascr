#include <admodel.h>
typedef dvariable (*invlinkfn_pointer)(const dvariable&);

// Inverse identity (identity).
dvariable invlinkfn_identity (const dvariable &x)
{
  return x;
}

// Inverse log (exponential).
dvariable invlinkfn_log (const dvariable &x)
{
  return mfexp(x);
}

// Inverse logit.
dvariable invlinkfn_logit (const dvariable &x)
{
  return mfexp(x)/(mfexp(x) + 1);
}

invlinkfn_pointer get_invlinkfn(int linkfn_id)
{     
  invlinkfn_pointer invlinkfn;    
  switch(linkfn_id){
    case 1: invlinkfn = invlinkfn_identity; break;
    case 2: invlinkfn = invlinkfn_log; break;
    case 3: invlinkfn = invlinkfn_logit; break;
  }
  return(invlinkfn) ;
}
