natural_number_check = function(session, key){
  output = TRUE
  for(s in unique(session)){
    tem_key = key[which(session == s)]
    seq_key = seq(length(unique(tem_key)))
    if(any(!tem_key %in% seq_key)) output = FALSE
  }
  
  return(output)
}

convert_natural_number = function(dat, is.animalID, which.convert){
  
  dat.na = subset(dat, is.na(ID))
  dat = subset(dat, !is.na(ID))
  
  if(which.convert == 'ID'){
    if(is.animalID){
      keys = paste(dat$session, dat$animal_ID, sep = '_')
    } else {
      keys = dat$session
    }
    uni_keys = unique(keys)
    output = vector('list', length(uni_keys))
    for(i in 1:length(uni_keys)){
      output[[i]] = dat[keys == uni_keys[i],,drop = FALSE]
      output[[i]]$ID = as.numeric(as.factor(output[[i]]$ID))
    }
    output = do.call('rbind', output)
  } else if(which.convert == 'animal_ID'){
    stopifnot(is.animalID)
    keys = dat$session
    uni_keys = unique(keys)
    output = vector('list', length(uni_keys))
    for(i in 1:length(uni_keys)){
      output[[i]] = dat[keys == uni_keys[i],,drop = FALSE]
      output[[i]]$animal_ID = as.numeric(as.factor(output[[i]]$animal_ID))
    }
    output = do.call('rbind', output)
  } else if(which.convert == 'both'){
    if(is.animalID){
      dat = convert_natural_number(dat, is.animalID, 'animal_ID')
    }
    dat = convert_natural_number(dat, is.animalID, 'ID')
    output = dat
  } else {
    stop('Invalid input.')
  }
  
  output = rbind(output, dat.na)
  return(output)

}


agg_sort = function(dat, obj, lst, f){
  if(length(obj) == 1){
    x = dat[[obj]]
  } else {
    x = dat[[obj[1]]]
    for(i in 2:length(obj)){
      x = paste(x, dat[[obj[i]]], sep = '---')
    }
  }
  
  y = vector('list', length(lst))
  names(y) = lst
  for(i in lst) y[[i]] = dat[[i]]
  output = aggregate(x, y, function(x) f(x))
  
  sort_by = vector('list', length(lst))
  for(i in 1:length(lst)) sort_by[[i]] = output[[lst[i]]]
  
  o = do.call('order', sort_by)
  output = output[o,]
  return(output)
}


extend_dat_check = function(dat, check_var, ori_dat, n.sessions, n.var, identical_flag){
  stopifnot(is(dat, 'data.frame'))
  if(!(check_var %in% colnames(dat))){
    stopifnot(nrow(dat) == n.var[1])
    dat[[check_var]] = seq(n.var[1])
  }
  
  if('session' %in% colnames(dat)){
    stopifnot(length(colnames(dat)) > 2)
  } else {
    #if all sessions share the same trap/mask, the user can skip 'session' in the traps/masks' covariates
    #or there is only 1 session
    stopifnot(length(colnames(dat)) > 1)
    if(n.sessions > 1) stopifnot(identical_flag)
    tem = vector('list', n.sessions)
    for(s in 1:n.sessions){
      tem[[s]] = dat
      tem[[s]][['session']] = s
    }
    dat = do.call('rbind', tem)
  }
  
  if(any(duplicated(dat[, c('session', check_var)]))) stop(paste0('duplicated ', check_var, ' input.'))
  
  stopifnot(all(paste(dat$session, dat[[check_var]], sep = '-') %in%
                  unique(paste(ori_dat$session, ori_dat[[check_var]], sep = '-'))))
  
  return(dat)
}



covariates_mask_check = function(dat, n.sessions, n.masks, identical_flag){
  stopifnot(is(dat, 'list') | is(dat, 'data.frame'))
  
  if(is(dat, 'list')){
    if(length(dat) > 1){
      #if more than one component in this list, each component represent one session
      stopifnot(length(dat) == n.sessions)
      stopifnot(all(sapply(dat, nrow) == n.masks))
    } else {
      #if only one component in this list, then if n.sessions > 1 is a special scenario
      stopifnot(nrow(dat[[1]]) == n.masks[1])
      
      if(n.sessions > 1){
        stopifnot(identical_flag)
        tem = dat[[1]]
        dat = vector('list', n.sessions)
        for(s in 1:n.sessions) dat[[s]] = tem
      }
      
    }
    
  } else {
    stopifnot(nrow(dat) == n.masks[1])
    if(n.sessions > 1){
      stopifnot(identical_flag)
    }
    
    tem = dat
    dat = vector('list', n.sessions)
    for(s in 1:n.sessions) dat[[s]] = tem
  }
  
  for(s in 1:n.sessions){
    dat[[s]][['session']] = s
    dat[[s]][['mask']] = 1:n.masks[s]
  }
  dat = do.call('rbind', dat)
  
  return(dat)
}


default.link = function(i){
  if(i == "mu") return('log')
  if(i == "D") return('log')
  if(i == "g0") return('logit')
  if(i == "sigma") return('log')
  if(i == "lambda0") return('log')
  if(i == "z") return('log')
  if(i == "shape.1") return('log')
  if(i == "shape.2") return('identity')
  if(i == "shape") return('identity')
  if(i == "scale") return('log')
  if(i == "b0.ss") return('log')
  if(i == "b1.ss") return('log')
  if(i == "b2.ss") return('log')
  if(i == "sigma.ss") return('log')
  if(i == "kappa") return('log')
  if(i == "alpha") return('log')
  if(i == "sigma.toa") return('log')
  if(i == "sigma.b0.ss") return('log')
}

detfn.params = function(detfn){
  if(detfn == 'hn'){
    param.og = c('g0', 'sigma')
  }
  if(detfn == 'hhn'){
    param.og = c('lambda0', 'sigma')
  }
  if(detfn == 'hr'){
    param.og = c('g0', 'sigma', 'z')
  }
  if(detfn == 'th'){
    param.og = c('shape', 'scale')
  }
  if(detfn == 'lth'){
    param.og = c('shape.1', 'shape.2', 'scale')
  }
  if(detfn == 'ss'){
    param.og = c('b0.ss', 'b1.ss', 'sigma.ss')
  }
  if(detfn == 'ss_dir'){
    param.og = c('b0.ss', 'b1.ss', 'sigma.ss', 'b2.ss')
  }
  if(detfn == 'ss_het'){
    param.og = c('b0.ss', 'b1.ss', 'sigma.ss', 'sigma.b0.ss')
  }
  return(param.og)
}

p.dot.defaultD = function(points = NULL, traps = NULL, detfn = NULL, ss_dir = NULL,
                          ss.link = NULL, pars = NULL, A, n.quadpoints = 8){
  
  dists <- distances(traps, points)
  
  if(ss_dir){
    n.traps <- nrow(traps)
    n.points <- nrow(points)
    dirs <- (0:(n.quadpoints - 1))*2*pi/n.quadpoints
    probs <- numeric(n.points)
    ## Integrating over all possible directions.
    ## in the future: Write all this in C++.
    for (i in 1:n.quadpoints){
      dir <- dirs[i]
      bearings <- bearings(traps, points)
      orientations <- abs(dir - bearings)
      for (j in 1:n.points){
        ## Probabilities of detection given orientation.
        o.prob <- numeric(n.traps)
        for (k in 1:n.traps){
          o.prob[k] <- det_prob(detfn, pars, dists[k, j], ss.link
            #direction ss is not supported yet, hide this part
                                  #      , orientations[k, j]
                                  )
        }
        probs[j] <- probs[j] + (1/n.quadpoints)*(1 - prod(1 - o.prob))
      }
    }
    out <- probs
  } else {
    probs <- det_prob(detfn, pars, dists, ss.link)
    if(detfn == 'ss'){
      sigma.ss <- pars$sigma.ss
      cutoff <- pars$cutoff
      probs = 1 - pnorm(cutoff, mean = probs, sd = sigma.ss)
    }
    out <- plyr::aaply(probs, 2, function(x) 1 - prod(1 - x))
  }
  
  out <- A*sum(out)
  
  return(out)
}


formula_separate = function(foo, var.m){
  foo_vars = all.vars(foo[[3]])
  response_var = all.vars(foo[[2]])
  index = which(foo_vars %in% var.m)
  if(length(index) == 0){
    full_vars = foo_vars
    mask_vars = character(0)
  } else {
    full_vars = foo_vars[-index]
    mask_vars = foo_vars[index]
  }
  
  foo_terms = attr(stats::terms(foo), 'factors')[-1,,drop = FALSE]
  
  row_names = as.list(rownames(foo_terms))
  col_names = as.list(colnames(foo_terms))
  #use col names to check whether there is interaction term
  #between variables from mask level and other levels
  for(i in 1:length(col_names)){
    tem.name = col_names[[i]]
    tem.foo = stats::as.formula(paste0('gam.resp~',tem.name))
    tem.foo_vars = all.vars(tem.foo[[3]])
    if(!(sum(tem.foo_vars %in% var.m) %in% c(0, length(tem.foo_vars)))){
      stop("interaction between variables from mask level and other level is not supported.")
    }
  }
  #use row names to separate this foo_terms matrix into two matrix
  #one for mask level and other one for other levels
  for(i in 1:length(row_names)){
    tem.name = row_names[[i]]
    tem.foo = stats::as.formula(paste0('gam.resp~',tem.name))
    tem.foo_vars = all.vars(tem.foo[[3]])
    row_names[[i]] = tem.foo_vars
  }
  
  len = lapply(row_names, length)
  foo_terms = foo_terms[rep(1:nrow(foo_terms), len),,drop = FALSE]
  row_names = do.call('c', row_names)
  rownames(foo_terms) = row_names
  
  n.var.full = sum(!(foo_vars %in% var.m))
  n.var.mask = sum((foo_vars %in% var.m))
  
  if(n.var.full > 0) {
    tem_full = vector('list', n.var.full)
    names(tem_full) = full_vars
  }
  
  if(n.var.mask > 0) {
    tem_mask = vector('list', n.var.mask)
    names(tem_mask) = mask_vars
  }
  for(i in foo_vars){
    if(i %in% full_vars){
      tem_full[[i]] = foo_terms[which(row_names == i),,drop = FALSE]
      tem_full[[i]] = apply(tem_full[[i]], 2, function(x) sum(x) > 0)
      tem_full[[i]] = matrix(tem_full[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_full[[i]]) = colnames(foo_terms)
      rownames(tem_full[[i]]) = i
    } else {
      tem_mask[[i]] = foo_terms[which(row_names == i),,drop = FALSE]
      tem_mask[[i]] = apply(tem_mask[[i]], 2, function(x) sum(x) > 0)
      tem_mask[[i]] = matrix(tem_mask[[i]], ncol = ncol(foo_terms),
                             byrow = T)
      colnames(tem_mask[[i]]) = colnames(foo_terms)
      rownames(tem_mask[[i]]) = i
    }
  }
  
  if(n.var.full > 0){
    tem = apply(do.call('rbind', tem_full), 2, any)
    tem.terms = names(tem)[which(tem)]
    foo_full = stats::as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_full = stats::as.formula(paste0(response_var, "~ 1"))
  }
  
  if(n.var.mask > 0){
    tem = apply(do.call('rbind', tem_mask), 2, any)
    tem.terms = names(tem)[which(tem)]
    foo_mask = stats::as.formula(paste0(response_var, "~", paste(tem.terms, collapse = "+")))
  } else {
    foo_mask = stats::as.formula(paste0(response_var, "~ 1"))
  }
  
  return(list(foo_full = foo_full, foo_mask = foo_mask))
}


#if the value is valid, return TRUE

param.range.validate = function(param, value){
  if(param %in% c("sigma", "lambda0", "z", "shape.1","scale", "sigma.b0.ss", "kappa", "alpha", "sigma.toa",
                  "b0.ss", "b1.ss", "b2.ss", "sigma.ss", 'D', 'mu')){
    if(value >= 0){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else if(param == 'g0'){
    if(value >= 0 & value <= 1){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}


#default sv, this part needs read the original function carefully

default.sv = function(param, info){
  data.full = info[['data.full']]
  data.mask = info[['data.mask']]
  buffer = info[['buffer']]
  A = info[['A']]
  ss.opts = info[['ss.opts']]
  cutoff = ss.opts[['cutoff']]
  ss.link = ss.opts[["ss.link"]]
  detfn = info[['detfn']]
  dims = info[["dims"]]
  survey.length = info[["survey.length"]]
  param.og = info[['param.og']]
  sv = info[['sv']]
  
  
  
  if(param == 'sigma'){
    easy.out = mean(buffer) / 4
    is.animal_ID = "animal_ID" %in% colnames(data.full)
    same.traplocs = (nrow(data.full[!duplicated(data.full[,c('trap_x', 'trap_y')]), ]) == 1)
    if(same.traplocs){
      return(easy.out)
    } else {
      tem = data.full[!duplicated(data.full[,c("session", "animal_ID"[is.animal_ID], "ID", "trap")]), ]
      tem = subset(tem, !is.na(bincapt) & bincapt == 1)
      u.id = unique(tem[, c('session', "animal_ID"[is.animal_ID], "ID")])
      u.id = paste(u.id$session, u.id$animal_ID[is.animal_ID], u.id$ID, sep = "_")
      
      mean.dists = data.frame(avg_dist = numeric(0), weight = numeric(0))
      j = 1
      for(i in 1:length(u.id)){
        tem1 = subset(tem, paste(tem$session, tem$animal_ID[is.animal_ID], tem$ID, sep = "_") == u.id[i],
                      select = c("trap_x", "trap_y"))
        
        if(nrow(tem1) > 1){
          rc.dists <- distances(as.matrix(tem1), as.matrix(tem1))
          mean.dists[j, 1] = mean(rc.dists[rc.dists > 0])
          mean.dists[j, 2] = length(rc.dists[rc.dists > 0])
          j = j + 1
        }
      }
      
      if(nrow(mean.dists) > 0){
        out <- sum(mean.dists$avg_dist * mean.dists$weight)/sum(mean.dists$weight)
      } else {
        out = easy.out
      }
      return(out)
    }
  } else if(param == "D"){
    session_to_use = which(dims$n.IDs > 0)[1]
    tem.mask = subset(data.mask, session == session_to_use)
    mask = as.matrix(unique(tem.mask[, c('x', 'y')]))
    tem.trap = subset(data.full, session == session_to_use)
    traps = tem.trap[, c('trap', 'trap_x', 'trap_y')]
    traps = traps[!duplicated(traps[['trap']]), ]
    traps = as.matrix(traps[order(traps$trap), c('trap_x', 'trap_y')])
    survey.length = survey.length[session_to_use]
    A = A[session_to_use]
    
    detpar.names = param.og[which(!(param.og %in% c('kappa', 'alpha', 'sigma.toa', 'D')))]
    pars <- sv[detpar.names]
    
    if(detfn == 'ss_het') ss_het = TRUE else ss_het = FALSE
    if(detfn == 'ss_dir') ss_dir = TRUE else ss_dir = FALSE
    
    if(ss_het | ss_dir | detfn == 'ss'){
      detfn = 'ss'
      pars[["sigma.b0.ss"]] <- 0
      if(!any(detpar.names == 'b2.ss')) pars[['b2.ss']] = 0
    }
    
    if (!is.null(cutoff)){
      pars$cutoff <- cutoff
    }
    esa <- p.dot.defaultD(points = mask, traps = traps, detfn = detfn, ss_dir = ss_dir,
                          ss.link = ss.link, pars = pars, A = A, n.quadpoints = 8)
    
    return(dims$n.IDs[session_to_use]/(esa*survey.length))
    #this part in the original function is difficult, better to ask for help
  } else if(param == "g0"){
    return(0.95)
  } else if(param == "lambda0"){
    return(2)
  } else if(param == "z"){
    return(1)
  } else if(param == "sigma.toa"){
    return(0.0025)
  } else if(param == "kappa"){
    return(10)
  } else if(param == "b0.ss"){
    if(ss.link == 'log'){
      return(log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(max(data.full[['ss']], na.rm = TRUE))
    }
  } else if(param == "sigma.b0.ss"){
    if(ss.link == 'log'){
      return(0.01 * log(max(data.full[['ss']], na.rm = TRUE)))
    } else {
      return(0.01 * max(data.full[['ss']], na.rm = TRUE))
    }
  } else if(param == "b1.ss"){
    max.ss = max(data.full[["ss"]], na.rm = TRUE)
    out = (max.ss - cutoff)/(mean(buffer)/2)
    if(ss.link == "log"){
      out = out/max.ss
    }
    return(out)
  } else if(param == "b2.ss"){
    return(0.1)
  } else if(param == "sigma.ss"){
    ss_valid = data.full[['ss']][data.full[['ss']] >= cutoff & !is.na(data.full[['ss']])]
    return(stats::sd(ss_valid))
  } else if(param == "alpha"){
    return(2)
  } else if(param == "shape"){
    return(2)
  } else if(param == "scale"){
    if(detfn == "th"){
      out = default.sv(param = "sigma", info = info)
    } else if(detfn == "lth"){
      sigma = default.sv(param = "sigma", info = info)
      shape.1 = default.sv(param = "shape.1", info = info)
      out = (log(shape.1) - log(shape.1 - 0.5)) / sigma
    } else {
      stop("Detection function not recognised.")
    }
    return(out)
  } else if(param == "shape.1"){
    return(0.809017)
  } else if(param == "shape.2"){
    sigma = default.sv(param = "sigma", info = info)
    shape.1 = default.sv(param = "shape.1", info = info)
    scale = default.sv(param = "scale", info = info)
    return(shape.1 + scale * sigma)
  } else if(param == 'mu'){
    n.animals = dims$n.animals
    n.animal.call = dims$n.animal.call
    n.sessions = dims$n.sessions
    
    i = 1
    for(s in 1:n.sessions){
      if(n.animals[s] > 0){
        for(a in 1:n.animals[s]){
          n.animal.call[i] = n.animal.call[i]/survey.length[s]
          i = i + 1
        }
      }
    }
    
    return(mean(n.animal.call))
  }
}


#default bounds
#the default bounds of "D" is not specified in the origin
default.bounds = function(param){
  if(param == 'g0') {
    return(c(0, 1))
  } else if(param %in% c('sigma', 'lambda0', 'shape.1', 'scale', 'b0.ss',
                         'b1.ss', 'b2.ss', 'sigma.b0.ss', 'sigma.ss',
                         'z', 'sigma.toa', 'alpha', 'mu')){
    return(c(0, 1e8))
  } else if(param == 'kappa'){
    return(c(0, 700))
  } else if(param %in% c('shape', 'shape.2')){
    return(c(-100, 100))
  } else if(param == "D"){
    return(c(0, 1e20))
  }
}

#link function
link.fun = function(link, value){
  if(!link %in% c("identity", "log", "logit")){
    stop('Not a valid link function is assigned.')
  }
  
  if(link == "identity") return(cus_identity(value))
  if(link == "log") return(cus_log(value))
  if(link == "logit") return(cus_logit(value))
  
}

unlink.fun = function(link, value){
  if(!link %in% c("identity", "log", "logit")){
    stop('Not a valid link function is assigned.')
  }
  
  if(link == "identity") return(value)
  if(link == "log") return(cus_log_unlink(value))
  if(link == "logit") return(cus_logit_unlink(value))
}

delta.fun = function(link, sd, est){
  if(link == "identity") return(sd)
  
  #for link == 'log', we are calculating std(exp(x)) with x_hat = est and x_sd = sd
  if(link == "log") return(sd * exp(est))
  
  #for link == 'logit' we are calculating std(exp(x)/(exp(x) + 1))
  if(link == 'logit'){
    fx_grad = exp(est)/(exp(est) + 1) ^ 2
    return(sd * fx_grad)
  }
}

sort.data = function(dat, name){
  is.animal_ID = "animal_ID" %in% colnames(dat)
  
  if(name == 'data.full'){
    if(is.animal_ID){
      dat = dat[order(dat$session, dat$animal_ID,
                      dat$ID, dat$trap), ]
    } else {
      dat = dat[order(dat$session, dat$ID, dat$trap), ]
    }
  }
  
  if(name == "data.dists.thetas"){
    dat = dat[order(dat$session, dat$mask, dat$trap), ]
  }
  
  if(name == "data.ID_mask"){
    if(is.animal_ID){
      #the data.ID_mask in the c++ template will be sorted as
      #session - animal_ID - mask - ID, which differs from this
      #but in the data preparation process, we use this order
      dat = dat[order(dat$session, dat$animal_ID, dat$ID, dat$mask),]
    } else {
      dat = dat[order(dat$session, dat$ID, dat$mask), ]
    }
  }
  
  if(name == "data.mask"){
    dat = dat[order(dat$session, dat$mask),]
  }
  
  if(name == "u_id_match"){
    dat = dat[order(dat$session, dat$u_id, dat$ID),]
  }
  
  if(name == "index_traps_uid"){
    dat = dat[order(dat$session, dat$u_id),]
  }
  
  if(name == "data_u_bin"){
    dat = dat[order(dat$session, dat$u_id, dat$trap),]
  }
  
  return(dat)
}

numeric_link = function(link){
  if(is.null(link)){
    ans = 0
  } else {
    #if no match, return 0
    ans = numeric(length(link))
    for(i in 1:length(link)){
      if(link[i] == 'identity') ans[i] = 1
      if(link[i] == 'log') ans[i] = 2
      if(link[i] == 'logit') ans[i] = 3
      if(link[i] == 'spherical') ans[i] = 4
    }
  }
  return(ans)
}

numeric_detfn = function(detfn){
  if(detfn == 'hn') return(1)
  if(detfn == 'hhn') return(2)
  if(detfn == 'hr') return(3)
  if(detfn == 'th') return(4)
  if(detfn == 'lth') return(5)
  if(detfn == 'ss') return(6)
  if(detfn == 'ss_dir') return(7)
  if(detfn == 'ss_het') return(8)
}

numeric_het_method = function(het_method){
  if(is.null(het_method)) {
    return(1)
  } else {
    if(het_method == "GH") return(2)
    if(het_method == "rect") return(3)
  }
}

cal_n_det = function(data, is_uid = FALSE){
  if(is_uid){
    data = subset(data, !is.na(u_id))
  } else {
    data = subset(data, !is.na(ID))
  }

  if("animal_ID" %in% colnames(data)){
    tem = aggregate(data$bincapt, list(session = data$session,
                                       animal_ID = data$animal_ID,
                                       ID = data$ID), sum)
    tem = tem[order(tem$session, tem$animal_ID, tem$ID),]
    return(tem$x)
  } else {
    if(is_uid){
      tem = aggregate(data$bincapt, list(session = data$session,
                                         u_id = data$u_id), sum)
      tem = tem[order(tem$session, tem$u_id),]
    } else {
      tem = aggregate(data$bincapt, list(session = data$session,
                                         ID = data$ID), sum)
      tem = tem[order(tem$session, tem$ID),]
    }
    return(tem$x)
  }
  
}


extract_unique_id = function(data.full, dims){
  dat = subset(data.full,!is.na(ID))
  data.u.id.match = data.frame(session = numeric(0), ID = numeric(0), u_id = numeric(0))
  n.u.id = numeric(0)
  n.u.id.session = numeric(dims$n.sessions)
  for(s in 1:dims$n.sessions){
    tem = subset(dat, session == s)
    if(nrow(tem) > 0){
      tem = aggregate(tem$bincapt, list(ID = tem$ID), function(x) t(x))
      u.bin = unique(tem$x)
      n.u.id.session[s] = nrow(u.bin)
      tem$u_id = 0
      for(i in 1:nrow(tem)){
        captbin = as.vector(tem[i,'x'])
        tem$u_id[i] = which(apply(u.bin, 1, function(x) all(as.vector(x) == captbin)))
      }
      tem = tem[, c('ID', 'u_id')]
      tem$session = s
      n.u.id = c(n.u.id, as.vector(table(tem$u_id)))
      data.u.id.match = rbind(data.u.id.match, tem)
      u.bin = data.frame(session = s, u_id = rep(1:nrow(u.bin), each = dims$n.traps[s]),
                         trap = rep(1:dims$n.traps[s], nrow(u.bin)), bincapt = as.vector(t(u.bin)))
    }
  }
  
  data.u.bin = merge(data.full, data.u.id.match, by = c('session', "ID"))
  #firstly sort it as data.full and record the duplicated indices based on "session-uid-trap"
  #this bool vector will be used outside of this function to create design matrices from the one we have
  #already built from data.full
  data.u.bin = sort.data(data.u.bin, "data.full")
  dup = duplicated(data.u.bin[,c('session', 'u_id', 'trap')])
  
  data.u.bin = data.u.bin[!dup, -which(colnames(data.u.bin) %in% c('bearing','dist','ss',
                                                                   'toa', "ID", 'mrds_x', 'mrds_y'))]
  
  return(list(data_u_bin = data.u.bin, uid_dup_data_full = dup, u_id_match = data.u.id.match,
              n_id_uid = n.u.id, n_uids = n.u.id.session))
}

cus_log = function(x){
  x <- pmax(x, .Machine$double.eps)
  return(log(x))
}

cus_logit = function(x){
  x <- pmax(x, .Machine$double.eps)
  x <- pmin(x, 1 - .Machine$double.eps)
  return(log(x / (1 - x)))
}

cus_identity = function(x){
  return(x)
}

cus_log_unlink = .Primitive('exp')

cus_logit_unlink = function(x){
  return(exp(x)/(exp(x) + 1))
}

scale.closure = function(var.ex.info){
  #although named as "numeric.cov" here, but for numeric covariates under scale = FALSE
  #this "numeric.cov" will still be FALSE
  numeric.cov = !sapply(var.ex.info, is.null)
  cov.names = names(var.ex.info)
  
  out.fun = function(covariates_df){
    cov.names.new <- colnames(covariates_df)
    out <- covariates_df
    for (i in cov.names){
      if (numeric.cov[i] & (i %in% cov.names.new)){
        out[, i] <- (out[, i] - var.ex.info[[i]][1])/var.ex.info[[i]][2]
      }
    }
    
    return(out)
  }
  
  return(out.fun)
}

#extract values without name from a vector

val = function(vec){
  names(vec) = NULL
  return(vec)
}

#calculate bearings
bearings_by_vec = function(trap_x, trap_y, mask_x, mask_y){
  n = length(trap_x)
  if(length(trap_y)!=n | length(mask_x)!=n | length(mask_y)!=n){
    stop("different length coordinates.")
  }
  
  x_diff = mask_x - trap_x
  y_diff = mask_y - trap_y
  
  output = atan(x_diff / y_diff)
  
  output = ifelse(y_diff < 0, output + pi, output)
  output = ifelse(y_diff >=0 & x_diff < 0, output + 2 * pi, output)
  
  return(output)
}


#duplicate the rows with k times if a specific column is k
#in order to split multiple animals/calls (item) from one row to multiple rows
#used in simulation function "sim_det_history"
split_item = function(dat, item){
  #for each n_call/n_animal > 1, we need to split it to several identical individual call
  tem_gt1 = subset(dat, dat[[item]] > 1)
  
  if(nrow(tem_gt1) > 0){
    tem_eq1 = subset(dat, dat[[item]] == 1)
    
    u_n_items = unique(tem_gt1[[item]])
    
    tem_gt1_list = vector('list', length(u_n_items))
    
    for(i in 1:length(u_n_items)){
      #extract the part with n_items equals this unique number of calls
      tem = subset(tem_gt1, tem_gt1[[item]] == u_n_items[i])
      #repeat this part u_n_items[i] times
      tem_list = vector('list', u_n_items[i])
      for(j in 1:u_n_items[i]) tem_list[[j]] = tem
      tem_gt1_list[[i]] = do.call('rbind', tem_list)
    }
    
    tem_gt1 = do.call('rbind', tem_gt1_list)
    
    output = rbind(tem_gt1, tem_eq1)
    #since all calls become individual calls, assign the item to be 1
    output[[item]] = 1
  } else {
    output = dat
  }
  
  return(output)
}


delta_method_ascr_tmb = function(cov_linked, param_values, link_funs = NULL, new_covariates = NULL,
                                 pars = NULL, name_og = NULL, name_extend = NULL, df_param = NULL,
                                 gam.model.full = NULL, gam.output = NULL){
  #delta method, G(x) is a vector of n functions, where x is a vector of m variables
  #then VAR(G(x)) = G'(x) %*% VAR(x) %*% t[G'(x)]
  #VAR(x) is m*m matrix, G'(x) is n*m matrix, G'(x)[i, j] = d(g[i])/d(x[j])
  
  #In this case, m is always the n_rows or n_cols of the cov_linked, which is equal to length of param_values
  m = length(param_values)
  
  #for simple case, no new covariates involved, in which case, "link_funs" should be provided
  #otherwise, "new_covariates", "name_og", "object" will be used together
  if(!is.null(link_funs)){
    #in this case, G'(x) is a n*n diagonal matrix
    n = length(link_funs)
    G_grad = matrix(0, nrow = n, ncol = n)
    for(i in 1:n){
      
      link = link_funs[i]
      x = param_values[i]
      if(link == 'identity'){
        G_grad[i, i] = 1
      } else if(link == 'log'){
        G_grad[i, i] = exp(x)
      } else if(link == 'logit'){
        G_grad[i, i] = exp(x)/(exp(x) + 1)^2
      }
      
    }
    
  } else {
    if(nrow(new_covariates) != 1) stop('Only 1 row in new_covariates is accepted.')
    #if "new_covariates" is provided, then we are going to calculate the covariance matrix
    #for all original parameters used in this model, let "n" be length(unique(name_og))
    #and the G'(x) will be a n * m matrix
    unique_name_og = unique(name_og)
    n = length(unique_name_og)
    
    G_grad = matrix(0, nrow = n, ncol = m)
    
    
    for(i in 1:n){
      #for each row, the funciont is G_i(x) = h^(-1)(t(w)x), where w is the covariates values
      #x is parameters estimated values, h^(-1)(.) is reverse function of the link function
      par_name = unique_name_og[i]
      #extract the indices of this parameter in name_og, which is the
      #original parameter names of each column of cov_linked
      index_p = which(name_og == par_name)
      x = param_values[index_p]
      
      
      tem = subset(df_param, par == par_name)
      link = tem[, 'link']
      n_col_full = tem[, 'n_col_full']
      n_col_mask = tem[, 'n_col_mask']
      
      #let w be the covariates value for each parameter
      w = numeric(n_col_full + n_col_mask)
      
      if(length(w) != length(x)){
        stop('Critical Error.')
      }
      
      if(par_name %in% name_extend & par_name %in% pars){
        
        if(n_col_full > 1){
          gam.model = gam.output[[par_name]][["gam_non_mask"]]
          w[1:n_col_full] = as.vector(get_DX_new_gam(gam.model, new_covariates))
        } else {
          w[1] = 1
        }
        
        if(n_col_mask > 0){
          gam.model = gam.output[[par_name]][["gam_mask"]]
          tem = get_DX_new_gam(gam.model, new_covariates)
          #get rid of the first column since the intercept is not here
          tem = tem[, -1, drop = FALSE]
          w[(n_col_full + 1):(n_col_full + n_col_mask)] = as.vector(tem)
        }
        
      } else {
        w = 1
      }
      
      len_p = length(index_p)
      
      if(link == 'identity'){
        for(j in 1:len_p){
          G_grad[i, index_p[j]] = 1
        }
      } else if(link == 'log'){
        common_component = exp(sum(w * x))
        for(j in 1:len_p){
          G_grad[i, index_p[j]] = common_component * w[j]
        }
      } else if(link == 'logit'){
        common_component = exp(sum(w * x))
        for(j in 1:len_p){
          G_grad[i, index_p[j]] = common_component * w[j] / (common_component + 1) ^ 2
        }
      }
      
    }
    
  }
  
  #for test purpose
  #print('for test use only, display the G\'(x) matrix')
  #print(G_grad)
  #print('end of test display')
  
  return(G_grad %*% cov_linked %*% t(G_grad))
  
}

vector_to_df = function(vec){
  name = names(vec)
  if(is.null(name)) name = seq(length(vec))
  output = data.frame(name = name, value = as.vector(vec))
  return(output)
}

confint_gaussian_cal = function(object, types, pars, new_covariates, q_lower, q_upper){
  df_est = vector_to_df(coef(object, types = types, pars = pars, new_covariates = new_covariates))
  df_std = vector_to_df(stdEr(object, types = types, pars = pars, new_covariates = new_covariates))
  
  output = merge(df_est, df_std, by = 'name', all.x = TRUE)
  colnames(output) = c('par', 'est', 'std')
  output$std = ifelse(is.na(output$std), 0, output$std)
  output$lower = output$est + q_lower * output$std
  output$upper = output$est + q_upper * output$std
  return(output)
}

ori_name = function(char){
  #this order matters, we must try to match 'sigma.b0.ss' and 'sigma.toa' first
  #otherwise "sigma.b0.ss.(Intercept)" will be recognized as 'sigma'
  #so we must put 'sigma' behind any 'sigma.xx' parameters
  fulllist.par = c('g0', 'lambda0', 'z', 'shape.1',
                   'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
                   'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
                   "sigma.b0.ss", 'D', 'sigma', 'mu')
  #the pattern is the original parameter names followed by a dot
  pattern = paste(fulllist.par, ".", sep = "")
  n = length(char)
  output = character(n)
  for(i in 1:n){
    pattern_matched = FALSE
    
    #firstly try to match the pattern like "kappa."
    for(j in 1:length(pattern)){
      tem = regexpr(pattern[j], char[i])
      if(tem == 1){
        pattern_matched = TRUE
        output[i] = fulllist.par[j]
        break
      }
    }
    
    #if we cannot find the pattern like "kappa.", then we check whether the input is just the original parameter name
    
    if(!pattern_matched){
      for(j in 1:length(fulllist.par)){
        if(char[i] == fulllist.par[j]){
          pattern_matched = TRUE
          output[i] = fulllist.par[j]
          break
        }
      }
    }
    
    if(!pattern_matched){
      msg = paste0("Cannot find original parameter's name for the input: ", char[i])
      stop(msg)
    }
    
    #end of loop for char[i]
  }
  
  return(output)
}

## Finding nearest mask point to an observed MRDS location.
find.nearest.mask <- function(locs, mask){
  dists <- distances(locs, mask)
  as.list(apply(dists, 1, function(x) which(x == min(x))[1]))
}

mean_diy = function(vec){
  if(length(vec) > 0){
    if(is(vec, 'numeric')){
      return(mean(vec, na.rm = TRUE))
    } else {
      return(vec[1])
    }
  } else {
    if(is(vec, 'character')){
      return(character(0))
    } else {
      return(numeric(0))
    }
  }
}
