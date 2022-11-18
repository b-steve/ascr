fulllist.par.generator = function(){
  #the order matters, cannot be changed.
  return(c('g0', 'sigma', 'lambda0', 'z', 'shape.1',
            'shape.2', 'shape', 'scale', 'b0.ss', 'b1.ss',
            'b2.ss', 'sigma.ss', 'kappa', 'alpha', 'sigma.toa',
            "sigma.b0.ss", 'D', 'mu'))
}

natural_number_check = function(session, key){
  output = TRUE
  for(s in unique(session)){
    tem_key = key[which(session == s)]
    seq_key = seq(length(unique(tem_key)))
    if(any(!tem_key %in% seq_key)) output = FALSE
  }
  
  return(output)
}

convert_natural_number = function(dat, is.animalID, which.convert = 'both'){
  
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


df_to_list = function(df, n.sessions){
  if(is(df, 'list')){
    stopifnot(length(df) == n.sessions)
    return(df)
  } else {
    stopifnot(any(is(df, 'data.frame'), is(df, 'matrix')))
    output = vector('list', n.sessions)
    for(s in 1:n.sessions) output[[s]] = df
    return(output)
  }
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
  animal.model = info[['animal.model']]
  
  
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
    if(animal.model){
      return(dims$n.animals[session_to_use]/(esa*survey.length))
    } else {
      return(dims$n.IDs[session_to_use]/(esa*survey.length))
    }
    
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

sort.data = function(dat,  name){
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
  x = pmax(x, cus_logit(0))
  x = pmin(x, cus_logit(1))
  return(exp(x)/(exp(x) + 1))
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



list_2vector_4value = function(param){
  tem = sapply(param, length)
  output = do.call('c', param)
  names(output) = rep(names(tem), tem)
  return(output)
}


ext_par_in_new_df = function(par_ext, new_covariates, extend_par_covariates){
  output = NULL
  if(!is.null(par_ext)){
    for(i in par_ext){
      if(all(extend_par_covariates[[i]] %in% colnames(new_covariates))){
        output = c(output, i)
      }
    }
  }

  return(output)
}




delta_method_ascr_tmb = function(cov_linked, param_values, link_funs = NULL, new_covariates = NULL,
                                 name_og = NULL, name_extend = NULL, par_ext_cov_provided = NULL, df_param = NULL,
                                 gam.model.full = NULL, gam.output = NULL, back_trans = TRUE){
  #delta method, G(x) is a vector of n functions, where x is a vector of m variables
  #then VAR(G(x)) = G'(x) %*% VAR(x) %*% t[G'(x)]
  #VAR(x) is m*m matrix, G'(x) is n*m matrix, G'(x)[i, j] = d(g[i])/d(x[j])
  
  #In this case, m is always the n_rows or n_cols of the cov_linked, which is equal to length of param_values
  m = length(param_values)
  
  #for simple case, no new covariates involved, in which case, "link_funs" should be provided
  #otherwise, "new_covariates", "name_og", "df_param" will be used together, and "link_funs" should be NULL
  if(!is.null(link_funs)){
    #in this case, G'(x) is a n*n diagonal matrix
    n = length(link_funs)
    G_grad = matrix(0, nrow = n, ncol = n)
    for(i in 1:n){
      if(back_trans){
        link = link_funs[i]
      } else {
        link = 'identity'
      }
      
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
    #if(nrow(new_covariates) != 1) stop('Only 1 row in new_covariates is accepted.')
    #if "new_covariates" is provided, then we are going to calculate the covariance matrix
    #for all original parameters used in this model, let "n" be length(unique(name_og))
    #and the G'(x) will be a n * m matrix
    unique_name_og = unique(name_og)
    n = length(unique_name_og)
    G_grad = vector('list', n)
    n_new_df = nrow(new_covariates)
    #G_grad = matrix(0, nrow = n, ncol = m)
    
    
    for(i in 1:n){
      #for each row, the funciont is G_i(x) = h^(-1)(t(w)x), where w is the covariates values
      #x is parameters estimated values, h^(-1)(.) is reverse function of the link function
      par_name = unique_name_og[i]
      #extract the indices of this parameter in name_og, which is the
      #original parameter names of each column of cov_linked
      index_p = which(name_og == par_name)
      x = param_values[index_p]
      len_p = length(index_p)
      
      tem = subset(df_param, par == par_name)
      if(back_trans){
        #when a par is extended but not all covariates provided, assign the link as 'identity' as well
        #becuase in this case, we will show such as "D_int_link", "D_beta1_link"
        if((par_name %in% name_extend) & !(par_name %in% par_ext_cov_provided)){
          link = 'identity'
        } else {
          link = tem[, 'link']
        }
        
      } else {
        #if we do not do back transformation, it means, for example "D" has "D_int" and "D_beta1",
        #then we just use delta method to calculate var(D_int + D_beta1 * x1), instead of
        #var(exp(D_int + D_beta1 * x1))
        link = 'identity'
      }
      
      n_col_full = tem[, 'n_col_full']
      n_col_mask = tem[, 'n_col_mask']
      
      #let w be the covariates value for each parameter
      #w = numeric(n_col_full + n_col_mask)
      
      if((par_name %in% name_extend) & (par_name %in% par_ext_cov_provided)){
        G_grad[[i]] = matrix(0, nrow = n_new_df, ncol = m)
        w = matrix(0, nrow = n_new_df, ncol = n_col_full + n_col_mask)
        if(n_col_full > 1){
          gam.model = gam.output[[par_name]][["gam_non_mask"]]
          #browser()
          w[, 1:n_col_full] = get_DX_new_gam(gam.model, new_covariates)
        } else {
          w[, 1] = matrix(1, nrow = n_new_df, ncol = 1)
        }
        
        if(n_col_mask > 0){
          gam.model = gam.output[[par_name]][["gam_mask"]]
          tem = get_DX_new_gam(gam.model, new_covariates)
          #get rid of the first column since the intercept is not here
          w[, (n_col_full + 1):(n_col_full + n_col_mask)] = tem[, -1]
        }
        
      } else {
        G_grad[[i]] = matrix(0, nrow = len_p, ncol = m)
        w = diag(len_p)
      }

      
      n_row = nrow(w)
      
      if(link == 'identity'){
        for(k in 1:n_row){
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] = w[k, j]
          }
        }
      } else if(link == 'log'){
        for(k in 1:n_row){
          common_component = exp(sum(w[k,] * x))
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] = common_component * w[k, j]
          }
        }

      } else if(link == 'logit'){
        for(k in 1:n_row){
          common_component = exp(sum(w[k,] * x))
          for(j in 1:len_p){
            G_grad[[i]][k, index_p[j]] = common_component * w[k, j] / (common_component + 1) ^ 2
          }
        }
      }
      
    }
    #browser()
    
    G_grad = do.call('rbind', G_grad)
  }
  
  #for test purpose
  #print('for test use only, display the G\'(x) matrix')
  #print(G_grad)
  #print('end of test display')
  #browser()
  return(G_grad %*% cov_linked %*% t(G_grad))
  
}


#since DX could be very large, we only calculate and return the diagonal element
delta_for_pred = function(DX, par_value_linked, vcov_matrix, log_scale){
  #only use for D, so only consider the situation that link function is log
  n = nrow(DX)
  p = ncol(DX)
  
  if(!log_scale){
    G = exp(DX %*% par_value_linked)
    G = matrix(rep(G, p), ncol = p)
    G = G * DX
  } else {
    G = DX
  }

  
  output = numeric(n)
  for(i in 1:n){
    tem_G = G[i, , drop = FALSE]
    output[i] = (tem_G %*% vcov_matrix %*% t(tem_G))[1,1]
  }
  
  return(output)
  
}


vcov_fixed_par_add = function(m, p, type){
  #when p is NULL, n_added will be zero
  n_added = length(p)
  if(n_added == 0){
    output = m
  } else {
    if(type == 'linked') p = paste(p, "link", sep = "_")
    
    m_new = matrix(0, nrow = nrow(m) + n_added, ncol = ncol(m) + n_added)
    
    new_name = c(rownames(m), p)
    rownames(m_new) = new_name
    colnames(m_new) = new_name
    m_new[1:nrow(m), 1:ncol(m)] = m
    
    dim_m = data.frame(par = ori_name(new_name), dim_name = new_name, order_og = 1:nrow(m_new))
    full_par = fulllist.par.generator()
    full_par = data.frame(par = full_par, o = 1:length(full_par))
    dim_m = merge(dim_m, full_par, by = 'par')
    dim_m = dim_m[order(dim_m$o, dim_m$order_og),,drop = FALSE]
    
    output = matrix(0, nrow = nrow(m_new), ncol = ncol(m_new))
    rownames(output) = dim_m$dim_name
    colnames(output) = dim_m$dim_name
    
    tem = matrix(0, nrow = nrow(m_new), ncol = ncol(m_new))
    for(i in 1:nrow(m_new)){
      i_m_new = dim_m$order_og[i]
      tem[i, ] = m_new[i_m_new, ]
    }
    
    for(i in 1:ncol(m_new)){
      i_m_new = dim_m$order_og[i]
      output[, i] = tem[, i_m_new]
    }
    
  }
  
  return(output)
  
}







vector_to_df = function(vec){
  name = names(vec)
  if(is.null(name)) name = seq(length(vec))
  output = data.frame(name = name, value = as.vector(vec))
  return(output)
}

confint_gaussian_cal = function(object, types, pars, new_covariates, q_lower, q_upper){
  
  output = vector('list', length(types))
  names(output) = types
  for(i in types){
    #if(i == 'linked'){
    #  df_est = vector_to_df(coef.ascr_tmb(object, types = 'linked', pars = pars))
    #  df_std = vector_to_df(stdEr.ascr_tmb(object, types = 'linked', pars = pars))
    #}
    
    if(i == 'fitted' | i == 'linked'){
      #browser()
      #"fitted" is just back transformed from "linked" confidence interval
      df_est = vector_to_df(coef.ascr_tmb(object, types = 'linked', pars = pars, new_covariates = new_covariates))
      df_std = vector_to_df(stdEr.ascr_tmb(object, types = 'linked', pars = pars, new_covariates = new_covariates, show_fixed_par = FALSE))
    }
    
    if(i == 'derived'){
      df_est = vector_to_df(coef.ascr_tmb(object, types = 'derived'))
      df_std = vector_to_df(stdEr.ascr_tmb(object, types = 'derived', show_fixed_par = FALSE))
    }
    #browser()
    colnames(df_est) = c('par', 'est')
    df_est$std = 0
    u_name = unique(df_est$par)
    tem = vector('list', length(u_name))
    names(tem) = u_name
    for(k in u_name){
      tem_est = subset(df_est, df_est$par == k)
      tem_std = subset(df_std, df_std$name == k)
      if(nrow(tem_std) > 0) tem_est$std = tem_std$value
      tem[[k]] = tem_est
    }
    #browser()
    tem = do.call('rbind', tem)
    tem$lower = tem$est + q_lower * tem$std
    tem$upper = tem$est + q_upper * tem$std
    
    
    if(i == 'fitted') tem$par = gsub("_link", "", tem$par)
    
    output[[i]] = tem
    
  }
  
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
  #the pattern is the original parameter names at start
  fulllist_reg_exp = gsub("\\.", "\\\\\\.", fulllist.par)
  pattern = paste("^", fulllist_reg_exp, sep = "")
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

#weight method refers to Shepard and Modified Shepard
#mask: the "mask" object generated by create.mask()

#loc_cov: a list or a data frame or a matrix. If it is a list, each element should be a data frame or a matrix
#contains at least 3 columns, 2 of them are "x" and "y" for location, and the rest of them are the location related
#covariates. It could be a list means it could include different covariates from different sets of locations.

#control_nn2: a list contains the arugments for the function RANN::nn2(), please refer to the help document in that
#package for further details

#control_weight: a list contains 3 elements control the weighted method. The function RANN:nn2() gives us a M * k matrix,
#where M is the number of mask points, and each row contains k points provided by "loc_cov" which are nearest to this mask point.
#for this k nearest points, we use their values of each covariates to do the imputation to obtain the estimated value of this mask point.
#considering the distances between this mask point and these k nearest points are different, we need a method to determine their weights.
#The argument "control_weight" is used to do so. Its 3 elements are "method", which could be "Shepard" or "Modified"; "q", which is used
#to control both methods; "r", which is used to control method "Modified".
#"Shepard" method: w ~ (1/d)^q
#"Modified" method: modified shepard method, w ~ [max(0, r - d)/(r * d)]^q


#control_char: a list contains 2 elements (weighted_nn, k), and be used to control character type covariates.
#1st element - "weighted_nn": a logical value with FALSE by default. If TRUE, then used the value which has the highest weight as the value for each mask,
#for example, by RANN::nn2() and control_wight argument, a mask has 3 nearest points with w = c(0.4, 0.3, 0.3), and the values of these 3 points
#are ("a", "b", "b"), then "a" has weight of 0.4 and "b" has 0.6 in total, so the value assigned to this mask will be "b".
#If FALSE, then KNN method will be used.

#2nd element - k: it is used for the KNN method when weighted_nn is FALSE.

location_cov_to_mask = function(mask, loc_cov, control_nn2 = NULL, control_weight = NULL, control_char = NULL){
  #declare function nn2 from RANN
  f = RANN::nn2
  
  stopifnot(is(mask, 'list'))
  n.sessions = length(mask)
  n.masks = sapply(mask, nrow)
  
  if(is(loc_cov, 'data.frame') | is(loc_cov, 'matrix')){
    stopifnot(all(c('x', 'y') %in% colnames(loc_cov)))
    stopifnot(ncol(loc_cov) > 2)
    n_loc_cov = 1
    name_cov = colnames(loc_cov)[-which(colnames(loc_cov) %in% c('x', 'y'))]
    #this is the index to indicate which covariate is extracted from which element of loc_cov
    index_name_cov = rep(1, length(name_cov))
    loc_cov = list(loc_cov)
  } else {
    stopifnot(is(loc_cov, 'list'))
    n_loc_cov = length(loc_cov)
    name_cov = vector('list', n_loc_cov)
    #this is the index to indicate which covariate is extracted from which element of loc_cov
    index_name_cov = vector('list', n_loc_cov)
    for(i in 1:n_loc_cov){
      stopifnot(all(c('x', 'y') %in% colnames(loc_cov[[i]])))
      stopifnot(ncol(loc_cov[[i]]) > 2)
      name_cov[[i]] = colnames(loc_cov[[i]])[-which(colnames(loc_cov[[i]]) %in% c('x', 'y'))]
      index_name_cov[[i]] = rep(i, length(name_cov[[i]]))
    }

    name_cov = do.call('c', name_cov)
    index_name_cov = do.call('c', index_name_cov)
  }

  if(any(duplicated(name_cov))) stop('duplicated covariates found.')
  

  #set default values
  #when user does not assign any value to these arguments
  if(is.null(control_nn2)) control_nn2 = list(k = 10)
  if(is.null(control_weight)) control_weight = list(q = 2, method = 'Shepard', r = NULL)
  if(is.null(control_char)) control_char = list(weighted_nn = FALSE, k = 1)

  #when user only assign part of components in these arguments
  if(is.null(control_nn2$k)) control_nn2$k = 10
  if(is.null(control_char$weighted_nn)) control_char$weighted_nn = FALSE
  if(is.null(control_char$k)) control_char$k = 1
  if(is.null(control_weight$method)) control_weight$method = 'Shepard'
  if(is.null(control_weight$q)) control_weight$q = 2
  
  output = vector('list', n.sessions)
  
  for(s in 1:n.sessions){
    tem_mask = mask[[s]]
    output[[s]] = data.frame(session = s, mask = 1:n.masks[s])

    o_nn2 = vector('list', n_loc_cov)
    w = vector('list', n_loc_cov)
    
    #each i denote one loc_cov data frame provided, it might contain multiple covariates
    for(i in 1:n_loc_cov){
      tem_cov = loc_cov[[i]]
      tem_nn2_par = control_nn2
      tem_nn2_par$k = min(tem_nn2_par$k, nrow(tem_cov))
      tem_nn2_par$data = tem_cov[, c('x', 'y')]
      tem_nn2_par$query = tem_mask

      #the range denoted by tem_nn2_par$data should cover all tem_nn2_par$query ideally, or
      #at least the boundary of tem_nn2_par$query should not exceed tem_nn2_par$data too far
      #here the term of "too far" means 20% of the range of xlim of data or ylim of data
      xlim_data = range(tem_nn2_par$data[, 'x'])
      xscale_data = diff(xlim_data)
      #when x cordinates are the same for data, such as the data is provided in a line
      #i'm not sure how to determine the "scale", just take the abs(x) temporarily,
      if(xscale_data == 0) xscale_data = abs(xlim_data[1])
      xlim_data_max = c(xlim_data[1] - 0.7 * xscale_data, xlim_data[2] + 0.7 * xscale_data)
      xlim_query = range(tem_nn2_par$query[,'x'])

      ylim_data = range(tem_nn2_par$data[, 'y'])
      yscale_data = diff(ylim_data)
      if(yscale_data == 0) yscale_data = abs(ylim_data[1])
      ylim_data_max = c(ylim_data[1] - 0.7 * yscale_data, ylim_data[2] + 0.7 * yscale_data)
      ylim_query = range(tem_nn2_par$query[,'y'])

      if(any(xlim_query[1] < xlim_data_max[1], xlim_query[2] > xlim_data_max[2],
            ylim_query[1] < ylim_data_max[1], ylim_query[2] > ylim_data_max[2])){
        
        warning(paste0("for ", i, "'th element of 'loc_cov', the masks generated by the traps ",
                        "of ", s, "'th session exceeds the area boundary provided too much, ",
                        "it is recommened to provide location related covariates with a larger range."))
        
      }

      o_nn2[[i]] = do.call('f', tem_nn2_par)

      #remove zero index which may happen when we do "radius" search
      
      #nn.idx contains the k nearest points indices in the loc_cov for each mask point
      #nn.dists contains the distances between each mask point to those k nearest points
      #both of them are M * k matrix, M is the number of mask points
      
      
      zero_index = which(apply(o_nn2[[i]]$nn.idx, 2, sum) == 0)
      if(length(zero_index) > 0){
        o_nn2[[i]]$nn.idx = o_nn2[[i]]$nn.idx[, -zero_index, drop = FALSE]
        o_nn2[[i]]$nn.dists = o_nn2[[i]]$nn.dists[,-zero_index, drop = FALSE]
        if(ncol(o_nn2[[i]]$nn.idx) == 0){
          stop('Please specify a larger radius for radius searching.')
        }
      }
      
      
      #for each mask point, the value of the covariates on each mask point is a weighted average of all covariates values provided
      #there are two ways to determine the weights
      #Shepard: w ~ (1/d)^q, where q is a provided constant
      #Modified Shepard: w ~ [max(0, r - d) / (r * d)] ^ q, where r and q are provided constants
      
      #calculate weight matrix (for any row of w, sum(w[i,]) === 1)
      if(control_weight$method == 'Shepard'){
        w[[i]] = (1/o_nn2[[i]]$nn.dists) ^ control_weight$q
        w[[i]] = w[[i]] / apply(w[[i]], 1, sum)
      } else if(control_weight$method == 'Modified'){
        if(is.null(control_weight$r)){
          stop('please specify a valid radius for Modified Shepard Weight by control_weight$r = xx')
        }
        w[[i]] = (max(0, control_weight$r - o_nn2[[i]]$nn.dists) / (control_weight$r * o_nn2[[i]]$nn.dists)) ^ control_weight$q
        w[[i]] = w[[i]] / apply(w[[i]], 1, sum)
      } else {
        stop('weight method only suppors "Shepard" and "Modified" currently.')
      }

    }


    for(j in name_cov){
      #firstly, locate the index of element from where we got this covariate
      i = index_name_cov[which(j == name_cov)]
      values = loc_cov[[i]][, j]
      
      if(is.numeric(values)){
        #generate a matrix with the same dimension as o_nn2[[i]]$nn.idx, and
        #mat_values[i, j] = values[o_nn2[[i]]$nn.idx[i, j]]
        mat_values = matrix(values[o_nn2[[i]]$nn.idx], nrow = n.masks[s])
        output[[s]][[j]] = apply(mat_values * w[[i]], 1, sum)
      } else {
        values = as.character(values)
        #generate a matrix with the same dimension as o_nn2[[i]]$nn.idx, and
        #mat_values[i, j] = values[o_nn2[[i]]$nn.idx[i, j]]
        mat_values = matrix(values[o_nn2[[i]]$nn.idx], nrow = n.masks[s])
        v = character(n.masks[s])
        
        if(control_char$weighted_nn){
          #for each row, select the element from the 'nearest neighbours' in terms of weight
          for(n in 1:n.masks[s]){
            tem_v = mat_values[n,]
            tem_w = w[[i]][n,]
            tem_agg = aggregate(tem_w, list(v = tem_v), sum)
            v[n] = tem_agg$v[which.max(tem_agg$x)]
          }
        } else {
          #if not use weighted closest, we could simply select the KNN method for classification

          #since the KNN is based on nn2's result,the k here must be smaller or equal to the number of
          #"closet points" generated by nn2
          stopifnot(control_char$k <= ncol(o_nn2[[i]]$nn.idx))

          #drop the (k + 1)'th and after columns
          mat_values = mat_values[, 1:control_char$k, drop = FALSE]


          if(control_char$k != 1){
            #for each row, select the value with the largest frequency from its k's cloest neighbors 
            v = apply(mat_values, 1, function(x) names(which.max(table(x))))
          } else {
            #if k == 1, then it means we only take the nearst neighbour, and mat_values is a n*1 matrix, super easy
            v = as.vector(mat_values)


          }
        }

        output[[s]][[j]] = v
      }
      
      
      #end of co-variate j
    }
    
    
    #end of session s
  }
  return(do.call('rbind', output))
  
}


par_extend_create = function(loc_cov = NULL, mask = NULL, control_convert_loc2mask = list(),
                             session_cov = NULL, trap_cov = NULL){
  
  
  if(any(!is.null(loc_cov), !is.null(session_cov), !is.null(trap_cov))){
    par.extend = list()
    
    #if location related covariates provided, convert it to mask-level data frame
    if(!is.null(loc_cov)){
      control_convert_loc2mask$loc_cov = loc_cov
      control_convert_loc2mask$mask = mask
      
      mask_cov = do.call('location_cov_to_mask', control_convert_loc2mask)
    } else {
      mask_cov = NULL
    }

    par.extend$data = list()
    par.extend$data$session = session_cov
    par.extend$data$trap = trap_cov
    par.extend$data$mask = mask_cov
  } else (
    par.extend = NULL
  )
  
  return(par.extend)
}


demo_traps = function(){
  o = data.frame(x = rep(c(0,5) ,each = 3), y = rep(c(0, 5, 10), 2))
  return(o)
}

demo_traps2 = function(){
  o = data.frame(x = rep(20, 3), y = c(20, 25, 30))
  return(o)
}

demo_session_cov = function(){
  o = data.frame(session = 1:3, weather = c('sunny', 'rainy', 'sunny'))
  return(o)
}

demo_trap_cov = function(){
  o = data.frame(trap = 1:6, brand = rep(c('sony', 'panasonic'), each = 3))
  return(o)
}

demo_trap_cov2 = function(){
  o = data.frame(trap = 1:3, brand = c('sony', 'panasonic', 'panasonic'))
  return(o)
}

demo_loc_cov = function(){
  o = data.frame(x = rep(c(-20, 2.5, 25), each = 3), y = rep(c(-20, 5, 30), 3),
                 noise = c(6, 10, 11, 7, 12, 10, 11, 9, 8),
                 forest_volumn = c(rep('high', 2), rep('median', 3), rep('low', 4)))
  return(o)
}


sim_args_generator = function(sim_name){
  output = list()
  #generate some common settings
  traps = demo_traps()
  control_create_mask = list(buffer = 30)
  par_extend_model = NULL
  survey.length = NULL
  ss.opts = NULL
  cue.rates = NULL
  n.sessions = NULL
  session_cov = demo_session_cov()
  trap_cov = demo_trap_cov()
  loc_cov = demo_loc_cov()
  
  
  if(sim_name == 'dist_hn'){
    param = list(g0 = 1, sigma = 5.33, alpha = 5.62, D = 2300)
    detfn = 'hn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
  } else if(sim_name == 'bearing_dist_hn'){
    param = list(g0 = 1, sigma = 5.39, kappa = 51.12, alpha = 5.02, D = 2265)
    detfn = 'hn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
  } else if(sim_name == 'bearing_hn'){
    param = list(g0 = 1, sigma = 5.53, kappa = 53.88, D = 2172)
    detfn = 'hn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
  } else if(sim_name == 'ihd'){
    param = list(g0 = 1, sigma = 5.38, D = c(8.3, -0.2, -0.1))
    detfn = 'hn'
    par_extend_model = list(D = ~forest_volumn)
    session_cov = NULL
    trap_cov = NULL    
  } else if(sim_name == 'ihd_ext'){
    param = list(g0 = 1, sigma = c(1.95, -0.25), D = c(8.3, -0.1))
    detfn = 'hn'
    par_extend_model = list(D = ~noise, sigma = ~brand)
    session_cov = NULL
  } else if(sim_name == 'mul_ses'){
    param = list(g0 = 0.7135, sigma = 3.3, kappa = 14.8, alpha = 3.77, D = 2533)
    traps = list(traps, demo_traps2())
    control_create_mask = list(buffer = 15)
    detfn = 'hn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
    
  } else if(sim_name == 'mul_ses_ext'){
    param = list(g0 = c(1.172, -0.8), sigma = c(1.1383, 0.1407), kappa = 14.8, alpha = 3.77, D = 2533)
    traps = list(traps, demo_traps2())
    control_create_mask = list(buffer = 15)
    detfn = 'hn'
    session_cov = session_cov[1:2,]
    trap_cov = rbind(trap_cov, demo_trap_cov2())
    trap_cov$session = c(rep(1, 6), rep(2, 3))
    loc_cov = NULL
    par_extend_model = list(g0 = ~weather, sigma = ~brand)
    
  } else if(sim_name == 'simple_hhn'){
    param = list(sigma = 3.66, lambda0 = 4.29, D = 2657)
    detfn = 'hhn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
    
  } else if(sim_name == 'hhn_cue'){
    param = list(sigma = 3.66, lambda0 = 4.29, D = 500)
    detfn = 'hhn'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
    cue.rates = c(5, 4, 11, 2, 3)
    survey.length = 1
  } else if(sim_name == 'simple_hr'){
    param = list(g0 = 0.88, sigma = 7, z = 7.5, D = 2665)
    detfn = 'hr'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
  } else if(sim_name == 'ss'){
    param = list(b0.ss = 89.07, b1.ss = 4.13, sigma.ss = 9.52, D = 2657)
    detfn = 'ss'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
    ss.opts = list(cutoff = 60)
  } else if(sim_name == 'ss_toa'){
    param = list(b0.ss = 90, b1.ss = 4, sigma.ss = 9.34, sigma.toa = 0.002, D = 2430)
    detfn = 'ss'
    session_cov = NULL
    trap_cov = NULL  
    loc_cov = NULL
    ss.opts = list(cutoff = 60)
  } else if(sim_name == 'ind_bearing_dist'){
    param = list(g0 = 0.68, sigma = 9.8, kappa = 10, alpha = c(0.56, 0.16), D = c(4.15, 0.62), mu = 8.8)
    detfn = 'hn'
    control_create_mask = list(buffer = 15)
    loc_cov = NULL
    par_extend_model = list(D = ~weather, alpha = ~brand)
    n.sessions = 3
    survey.length = c(1, 2, 1)
  } else if(sim_name == 'ind_toa_hhn'){
    param = list(sigma = 2.13, lambda0 = 8.33, sigma.toa = 0.0011, D = c(6.81, -0.12), mu = 8.2)
    detfn = 'hhn'
    control_create_mask = list(buffer = 15)
    session_cov = NULL
    trap_cov = NULL
    par_extend_model = list(D = ~noise)
    n.sessions = 2
    
  } else if(sim_name == 'ind_ss'){
    param = list(b0.ss = c(4.3, 0.1), b1.ss = 0.8872, sigma.ss = 10.1, D = c(4.3, -0.1), mu = 8.3)
    detfn = 'ss'
    control_create_mask = list(buffer = 15)
    session_cov = NULL
    par_extend_model = list(D=~noise, b0.ss=~brand)
    n.sessions = 3
    ss.opts = list(cutoff = 60)
  } else if(sim_name == 'ind_ss_log'){
    ##################################################################################################
    #seems that the estimation is biased, check it further
    param = list(b0.ss = 5.35, b1.ss = 0.08, sigma.ss = 9.788, D = 67, mu = 8.3)
    detfn = 'ss'
    control_create_mask = list(buffer = 15)
    session_cov = NULL
    trap_cov = NULL
    loc_cov = NULL
    n.sessions = 3
    ss.opts = list(cutoff = 75, ss.link = 'log')
    ####################################################################################################
  } else if(sim_name == 'ind_ss_sp'){
    ##################################################################################################
    #seems that the estimation is biased, check it further
    param = list(b0.ss = 100.25, b1.ss = 0.082, sigma.ss = 9.788, D = 52, mu = 9.1)
    detfn = 'ss'
    control_create_mask = list(buffer = 15)
    session_cov = NULL
    trap_cov = NULL
    loc_cov = NULL
    n.sessions = 3
    ss.opts = list(cutoff = 75, ss.link = 'spherical')
    ##################################################################################################
  }
  
  
  #include all arguments that may differ from default
  output$detfn = detfn
  output$param = param
  output$par_extend_model = par_extend_model
  output$traps = traps
  output$control_create_mask = control_create_mask
  output$session_cov = session_cov
  output$trap_cov = trap_cov
  output$loc_cov = loc_cov
  output$survey.length = survey.length
  output$ss.opts = ss.opts
  output$cue.rates = cue.rates
  output$n.sessions = n.sessions
  
  return(output)
  
}


fit_args_generator_from_sim = function(sim_name, fit_args){

  output = fit_args

  if(sim_name == 'dist_hn'){
    output$fix = list(g0 = 1)
  } else if(sim_name == 'bearing_dist_hn'){
    output$fix = list(g0 = 1)
  } else if(sim_name == 'bearing_hn'){
    output$fix = list(g0 = 1)
  } else if(sim_name == 'ihd'){
    output$fix = list(g0 = 1)
  } else if(sim_name == 'ihd_ext'){
    output$fix = list(g0 = 1)
  } else if(sim_name == 'mul_ses'){
    output$sv = list(kappa = 20)
  } else if(sim_name == 'mul_ses_ext'){
    output$sv = list(kappa = 20)
  } else if(sim_name == 'simple_hhn'){
    #nothing need to be done
  } else if(sim_name == 'hhn_cue'){
    #nothing need to be done
  } else if(sim_name == 'simple_hr'){
    output$sv = list(z = 5)
  } else if(sim_name == 'ss'){
    output$sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10)
  } else if(sim_name == 'ss_toa'){
    output$sv = list(b0.ss = 90, b1.ss = 4, sigma.ss = 10)
  } else if(sim_name == 'ind_bearing_dist'){
    #nothing need to be done
  } else if(sim_name == 'ind_toa_hhn'){
    #nothing need to be done
  } else if(sim_name == 'ind_ss'){
    #nothing need to be done
  } else if(sim_name == 'ind_ss_log'){
    #nothing need to be done
  } else if(sim_name == 'ind_ss_sp'){
    #nothing need to be done
  }
  
  return(output)
}


param_transform = function(param, df_link, back = FALSE){
  par_names = names(param)
  which_trans = par_names[which(sapply(param, function(x) length(x) == 1))]
  
  if(back){
    for(i in which_trans){
      link = df_link[which(df_link$par == i), 'link']
      param[[i]] = unlink.fun(link, param[[i]])
    }
  } else {
    for(i in which_trans){
      link = df_link[which(df_link$par == i), 'link']
      param[[i]] = link.fun(link, param[[i]])
    }
  }
  
  
  return(param)
}


default_df_link = function(){
  pars = fulllist.par.generator()
  output = data.frame(par = pars, link = character(length(pars)))
  for(i in 1:length(pars)){
    output[i,'link'] = default.link(pars[i])
  }
  return(output)
}


diag_block_combine = function(lst){
  n_lst = sapply(lst, nrow)
  n = sum(n_lst)
  out = matrix(0, nrow = n, ncol = n)
  
  start_point = 1
  
  for(i in 1:length(lst)){
    end_point = start_point + n_lst[i] - 1
    indices = seq(from = start_point, to = end_point)
    out[indices, indices] = lst[[i]]
    start_point = start_point + n_lst[i]
  }
  
  return(out)
  
}


scale_convert = function(from, to){
  max_to = max(to)
  min_to = min(to)
  
  if(length(from) == 1){
    output = 0.5 * (max_to + min_to)
  } else {
    min_from = min(from)
    dist_og = max(from) - min_from
    dist_to = max_to - min_to
    if(dist_og == 0){
      output = 0.5 * (max_to + min_to)
    } else {
      output = (from - min_from) * dist_to / dist_og + min_to
    }
  }


  return(output)
}



circleFun <- function(centre = data.frame(x = 0, y = 0), r = 0.5, npoints = 50){
  n = nrow(centre)
  if(length(r) != n) r = rep(r, length.out = n)
  output = vector('list', n)
  for(i in 1:n){
    tt <- seq(0, 2*pi, length.out = npoints)
    xx <- centre[i, 1] + r[i] * cos(tt)
    yy <- centre[i, 2] + r[i] * sin(tt)
    output[[i]] = data.frame(x = xx, y = yy, cir_index = i)
  }
  
  output = do.call('rbind', output)
  
  return(output)
}


homo_digit = function(x){
  x = as.character(x)
  n_digit = nchar(x)
  max_digit = max(n_digit)
  
  n_zero = max_digit - n_digit
  
  y = character(length(x))
  for(i in 1:length(x)){
    y[i] = paste(rep('0', n_zero[i]), collapse = "")
  }
  output = paste(y, x, sep = "")
  return(output)
}




#predict density with information of locations coordinates, only used in the plot currently

predict_D_for_plot = function(fit, session_select = 1, new_data = NULL, D_cov = NULL, xlim = NULL, ylim = NULL,
                              x_pixels = 50, y_pixels = 50, se_fit = FALSE, log_scale = FALSE, set_zero = NULL, 
                              control_convert_loc2mask = NULL,
                              control_boot = list(correct_bias = FALSE, from_boot = TRUE), ...){
  
  
  if(!is.null(set_zero) & se_fit){
    warning('During the calculation of standard error, std_zero will be ignored.')
  }
  
  original_mask = FALSE
  
  if(is.null(new_data)){
    
    #if new_data is not provided, but xlim and ylim provided, we use xlim and ylim to build mask instead of 
    #extracting mask from fit
    if(any(!is.null(xlim), !is.null(ylim))){
      if(!all(!is.null(xlim), !is.null(ylim))) stop('please provide both xlim and ylim.')
      x = seq(from = xlim[1], to = xlim[2], length.out = x_pixels)
      y = seq(from = ylim[1], to = ylim[2], length.out = y_pixels)
      mask = data.frame(x = rep(x, each = y_pixels), y = rep(y, x_pixels))
    } else {
      mask = get_mask(fit)[[session_select]]
      original_mask = TRUE
    }
    
  } else {
    #if new_data is provided, then xlim and ylim will just do what they are supposed to do,
    #to trim the plot instead of building "mask"
    stopifnot(any(is(new_data, 'data.frame'), is(new_data, 'matrix')))
    stopifnot(all(c('x', 'y') %in% colnames(new_data)))
    mask = new_data[, c('x', 'y')]
  }
  
  
  if('D' %in% get_par_extend_name(fit)){
    
    #build the old_covariates
    ##we interpolate it again no matter there is new mask grid or not because in theory, user
    ##could use the same mask grid but different control_convert_loc2mask
    
    old_covariates = as.data.frame(mask)
    
    old_loc_cov = get_loc_cov(fit)
    
    #when we cannot find loc_cov from the model fitting object, it is possible the D is
    #not mask-level extended, it is also possible that the model fitting object comes
    #from simulation_study
    mask_level_dat_extract = FALSE
    if(is.null(old_loc_cov)){
      old_loc_cov = get_par_extend_data(fit)$mask
      if(!is.null(old_loc_cov)){
        mask_level_dat_extract = TRUE
        if('session' %in% colnames(old_loc_cov)){
          old_loc_cov = subset(old_loc_cov, session == session_select)
        }
        
        #when locations/mask in prediction is assigned by xlim/ylim or new_data, we use
        #the original mask_level data as loc_cov, otherwise, we do not need to do
        #conversion in the first place
        if(!original_mask){
          
          if('session' %in% colnames(old_loc_cov)){
            old_loc_cov = old_loc_cov[,which(colnames(old_loc_cov)!='session'),drop = FALSE]
          }
          if('mask' %in% colnames(old_loc_cov)){
            old_loc_cov = old_loc_cov[,which(colnames(old_loc_cov)!='mask'),drop = FALSE]
          }
          
          tem_mask = get_mask(fit)
          if(is(tem_mask, 'list')) tem_mask = tem_mask[[session_select]]
          stopifnot(nrow(old_loc_cov) == nrow(tem_mask))
          old_loc_cov = cbind(tem_mask, old_loc_cov)
          
        }
      }
      
    }
    
    
    if(!is.null(old_loc_cov)){
      #like mentioned right above, if we are using original mask data, and mask-level covariates data
      #we do not need to do conversion at all
      if(mask_level_dat_extract & original_mask){
        cov_mask = old_loc_cov
      } else {
        if(is.null(control_convert_loc2mask)){
          control_convert_loc2mask = vector('list', 2)
          names(control_convert_loc2mask) = c('mask', 'loc_cov')
        }
        control_convert_loc2mask$mask = list(mask)
        control_convert_loc2mask$loc_cov = old_loc_cov
        
        cov_mask = do.call('location_cov_to_mask', control_convert_loc2mask)
      }
      
      if(any(colnames(cov_mask) %in% c('session', 'mask'))){
        old_covariates = cbind(old_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      } else {
        old_covariates = cbind(old_covariates, cov_mask)
      }
      
    }
    
    
    
    session_data_in_model = get_par_extend_data(fit)$session
    ##for session related covariates, it is simple, just take them from the input argument of fit
    
    if(!is.null(session_data_in_model)){
      tem = session_data_in_model[which(session_data_in_model$session == session_select), , drop = FALSE]
      tem = tem[, -which(colnames(tem) == 'session'), drop = FALSE]
      for(i in colnames(tem)) old_covariates[[i]] = tem[1,i]
    }
    
    
    if(!is.null(D_cov) | !is.null(new_data)){
      #firstly deal with all scenarios that there may be any new covariate provided 
      if(!is.null(D_cov)){
        stopifnot(is(D_cov, 'list'))
        stopifnot(any(c('session', 'location') %in% names(D_cov)))
      }
      
      #considering the possibility that user may only want to change part of covariates and remains other as the same
      #as the model, for example, 3 covariates related to D, 'weather' (session related), 'noise' and 'forest_type'(loc related)
      #and user only want to change weather, or to change noise, or to change noise and forest_type, and etc.
      
      #so we create 2 data frame, one for all new covariates, and one for old covariates. We must make sure these two data frame
      #contains the same mask points, and then replace the covariates in the "old" data frame with the same column in the "new",
      #if the same covariate appears in the "new" data frame.
      if(!is.null(new_data)){
        #in new_data, user could include any location related covariates directly, and it also contains x and y,
        #so we could directly use it instead of mask
        new_covariates = as.data.frame(new_data)
      } else {
        new_covariates = as.data.frame(mask)
      }
      
      
      
      #build the new_covariates based on all information we could have
      if(!is.null(D_cov$location)){
        if(is.null(control_convert_loc2mask)){
          control_convert_loc2mask = vector('list', 2)
          names(control_convert_loc2mask) = c('mask', 'loc_cov')
        }
        control_convert_loc2mask$mask = list(mask)
        control_convert_loc2mask$loc_cov = D_cov$location
        
        
        cov_mask = do.call('location_cov_to_mask', control_convert_loc2mask)
        new_covariates = cbind(new_covariates, cov_mask[, -which(colnames(cov_mask) %in% c('session', 'mask')), drop = FALSE])
      }
      
      
      #since we only plot one session, the number of row for D_cov$session should be only 1
      if(!is.null(D_cov$session)){
        stopifnot(nrow(D_cov$session) == 1)
        for(i in colnames(D_cov$session)) new_covariates[[i]] = D_cov$session[1,i]
        
      }
      
      #update the old_covariates by the new_covariates
      
      for(i in colnames(old_covariates)){
        if(all(i != 'x', i != 'y', i %in% colnames(new_covariates))){
          old_covariates[[i]] = new_covariates[[i]]
        }
      }
    }
    
    par_info = get_data_param(fit)
    par_info = subset(par_info, par == 'D')
    gam.model = get_gam(fit, 'D')
    if(is(fit, 'ascr_boot')){
      values_link = as.vector(coef(fit, types = 'linked', pars = 'D', correct_bias = control_boot$correct_bias))
    } else {
      values_link = as.vector(coef(fit, types = 'linked', pars = 'D'))
    }
    
    values_link_og = values_link
    
    if(!is.null(set_zero)) values_link[set_zero] = 0
    
    
    tem = get_extended_par_value(gam.model, par_info$n_col_full, par_info$n_col_mask,
                                 values_link, old_covariates, DX_output = TRUE)
    D.mask = unlink.fun(link = par_info$link, value = tem$output)
    
    if(se_fit){
      DX = tem$DX
      
      if(is(fit, 'ascr_boot')){
        if(log_scale){
          type = 'linked'
        } else {
          type = 'fitted'
        }
        
        D.se = as.vector(stdEr(fit, types = type, new_covariates = old_covariates, pars = 'D',
                               from_boot = control_boot$from_boot, show_fixed_par = FALSE))
        
      } else {
        vcov_matrix = vcov(fit, types = 'linked', pars = 'D', show_fixed_par = FALSE)
        D.se = sqrt(delta_for_pred(DX, values_link_og, vcov_matrix, log_scale))
      }
      

    }
    
  } else {
    #when D is not extended, it is literally a constant, just take the 'fitted' estimated D by coef() and stdEr()
    D.mask <- rep(coef(fit, types = "fitted", pars = 'D'), nrow(mask))
    if(se_fit){
      if(log_scale){
        type = "linked"
      } else {
        type = "fitted"
      }

      D.se = as.vector(rep(stdEr(fit, types = type, pars = 'D', show_fixed_par = FALSE), nrow(mask)))
    }
    
  }
  
  output = as.data.frame(mask)
  output$est = D.mask
  if(se_fit) output$std = D.se
  
  
  return(output)
}



res_split = function(res, n.sessions){
  n_pars = ncol(res) - 1 - n.sessions
  ## extracting and removing maximum gradient component.
  maxgrads <- res[, n_pars + 1]
  res_esa = res[, (n_pars + 2):ncol(res), drop = FALSE]
  res <- res[, 1:n_pars, drop = FALSE]
  
  
  return(list(res = res, maxgrads = maxgrads, res_esa = res_esa))
}


var_from_res = function(res, fixed = NULL){
  par_names = colnames(res)
  
  if(!grepl('^esa', par_names[1])){
    to_be_rm = which(ori_name(par_names) %in% fixed)
    if(length(to_be_rm) != 0){
      res = res[,-to_be_rm, drop = FALSE]
      par_names = par_names[-to_be_rm]
    }
  }

  n_pars = ncol(res)
  
  se <- apply(res, 2, sd, na.rm = TRUE)
  names(se) <- par_names
  
  corr <- diag(n_pars)
  dimnames(corr) <- list(par_names, par_names)
  vcov <- diag(n_pars)
  diag(vcov) = se^2
  dimnames(vcov) <- list(par_names, par_names)
  if(n_pars > 1){
    for (i in 1:(n_pars - 1)){
      for (j in (i + 1):n_pars){
        corr[i, j] <- corr[j, i] <- cor(res[, i], res[, j], use = "na.or.complete")
        vcov[i, j] <- vcov[j, i] <- corr[i, j]*se[i]*se[j]
      }
    }
  }
  
  return(vcov)
  
}


mce_from_res = function(res, coefs, M, seed_mce){
  n_pars = ncol(res)
  par_names = colnames(res)
  bias <- apply(res, 2, mean, na.rm = TRUE) - as.vector(coefs)
  ## Bootstrap to calculate MCE for bias and standard errors.
  bias.mce <- numeric(n_pars)
  se.mce <- numeric(n_pars)
  names(bias.mce) <- par_names
  names(se.mce) <- par_names
  if (M > 0){
    set.seed(seed_mce)
    converged <- which(!is.na(res[, 1]))
    n.converged <- length(converged)
    mce.boot <- matrix(sample(converged, size = n.converged*M,
                              replace = TRUE), nrow = M,
                       ncol = n.converged)
    #browser()
    for (i in par_names){
      par.boot <- matrix(res[mce.boot, i], nrow = M, ncol = n.converged)
      bias.mce[i] <- sd(apply(par.boot, 1, mean) - coefs[i] - bias[i])
      se.mce[i] <- sd(apply(par.boot, 1, sd))
    }
  } else {
    bias.mce <- NA
    se.mce <- NA
  }
  
  return(list(bias = bias, bias.mce = bias.mce, se.mce = se.mce))
}


#used for solving the issue of determination of "types" and "pars" in kinds of methods
types_pars_sol = function(types, pars, new_covariates){

  #if new_covariates is provided, we need "fitted"
  #if(!is.null(new_covariates) & (!"fitted" %in% types)) types = c(types, 'fitted')
  #if types is still nothing, set it as 'linked'
  if(is.null(types)) types = 'linked'
  
  #confirm types is within the valid set of options
  if(any(!types %in% c('all', 'fitted', 'linked', 'derived'))){
    stop("Argument 'types' must be a subset of {'fitted', 'linked', 'derived', 'all'}.")
  } 
  
  if('esa' %in% pars){
    pars = pars[-which(pars == 'esa')]
    #if pars = 'esa' only, then regard it as NULL, force types to 'derived' only no matter
    #what it is assigned.
    if(length(pars) == 0){
      pars = NULL
      types = 'derived'
    } else {
      
      #if user assigns any other "pars" and only assigns types = "derived", add "linked" to types
      if(length(types) == 1){
        if(types == 'derived') types = c(types, 'linked')
      }
      if(!'derived' %in% types) types = c(types, 'derived')
    }
  }
  
  #after 'esa' modification, if we still have assigned 'pars', but types is assigned as
  #"derived" only, we need to add the default setting "linked" to it
  #together with the 'esa' modification above, the basic logic is that "pars" has priority
  #comparing to "types"
  if(length(pars) > 0 & all(types == 'derived')) types = c(types, 'linked')

  if ("all" %in% types){
    types <- c("fitted", "derived", "linked")
  }

  return(list(types = types, pars = pars))
  
}

################################################################################################################################


boot_res_to_CI = function(res, level){
  n_par = ncol(res)
  
  p_upper = 0.5 + 0.5 * level
  p_lower = 0.5 - 0.5 * level
  
  output = matrix(0, nrow = n_par, ncol = 2)
  colnames(output) = paste(c(p_lower, p_upper) * 100, "%", sep = " ")
  rownames(output) = colnames(res)
  
  for(i in 1:n_par){
    output[i,] = quantile(res[,i], c(p_lower, p_upper))
  }
  
  return(output)
}

#transform the linked boot result to fitted type
res_transform = function(res, new_covariates, pars, object, back_trans = TRUE){
  name_extend = get_par_extend_name(object)
  df_param = get_data_param(object)
  col_name_og = ori_name(colnames(res))
  output = vector('list', length(pars))
  names(output) = pars
  
  extend_par_covariates = get_par_extend_covariate(object)
  
  if(is.null(name_extend) & !is.null(new_covariates)){
    warning('No parameter is extended, argument "new_covariates" will be ignored.')
    new_covariates = NULL
  }
  
  name_extned_covariate_provided = ext_par_in_new_df(name_extend, new_covariates, extend_par_covariates)
  
  
  for(j in pars){
    values_link = res[, which(col_name_og == j), drop = FALSE]
    par_info = subset(df_param, par == j)
    link = par_info$link
    
    if(j %in% name_extend){

      if(!is.null(new_covariates) & (j %in% name_extned_covariate_provided)){
        gam = get_gam(object, j)
        values_fitted = get_extended_par_value(gam, par_info$n_col_full,
                                               par_info$n_col_mask, values_link, new_covariates,
                                               matrix_par_value = TRUE)
        colnames(values_fitted) = rep(j, ncol(values_fitted))
        if(!back_trans){
          colnames(values_fitted) = paste(colnames(values_fitted), "link", sep = "_")
        }
      } else {
          #if there is no new covariates assigned to this extended parameter (maybe no new_covariates at all,
          #or at least one of the covariates for this parameter is not provided), 
          #then regards all of the relevant beta as identity linked parameters
          values_fitted = values_link
          if(back_trans){
            colnames(values_fitted) = gsub("_link$", "", colnames(values_fitted))
          }
          link = 'identity'
      }
      
      
    } else {
      values_fitted = values_link
      if(back_trans){
        colnames(values_fitted) = j
      } else {
        colnames(values_fitted) = paste(j, "link", sep = "_")
      }
      
    }
    
    if(back_trans){
      output[[j]] = unlink.fun(link = link, value = values_fitted)
    } else {
      output[[j]] = values_fitted
    }
    
  }
  
  output = do.call('cbind', output)
  return(output)
}

#a method to do bias correction for parametric bootstrap
#see https://math.mit.edu/~dav/05.dir/class25-slides-all.pdf
res_mod_for_CI = function(res, est, correct_bias){
  if(correct_bias){
    res = -1 * res
    res = sweep(res, 2, 2 * est, FUN = "+")
  }
  return(res)
}

