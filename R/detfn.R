#dx and each component of det_par should be vectors with the same length
#of course, if any one of them is only a scaler, there will be no problem

#argument orientation is used for ss.dir, not supported yet

det_prob = function(det_fn, det_par, dx, ss.link = NULL, orientation = NULL){
  
  for(i in names(det_par)) assign(i, det_par[[i]])
  
  if(det_fn == 'hn'){
    out = g0 * exp((-1 * dx^2) / (2 * sigma^2))
  } else if(det_fn == 'hhn'){
    out = 1 - exp((-1) * lambda0 * exp(dx^2 * (-0.5) / sigma^2))
  } else if(det_fn == 'hr'){
    out = g0 * (1 - exp((-1) * (dx/sigma)^(-z)))
  } else if(det_fn == 'th'){
    out = 0.5 - 0.5 * erf(dx / scale - shape)
  } else if(det_fn == 'lth'){
    out = 0.5 - 0.5 * erf(shape.1 - exp(shape.2 - scale * dx))
  } else if(det_fn == 'ss'){
    #ss is a little different, if the simulated ss > cut off, the detection probability is 1
    #and 0 otherwise, so we return E(SS) instead of detected probability
    if(ss.link == 'identity'){
      out = b0.ss - b1.ss * dx
    } else if(ss.link == 'log'){
      out = exp(b0.ss - b1.ss * dx)
    } else if(ss.link == 'spherical'){
      out = ifelse(dx > 1, b0.ss - 20 * log10(dx) - b1.ss * (dx - 1), b0.ss)
    }
  }
  
  return(out)
  
}



