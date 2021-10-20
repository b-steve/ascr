###Interpolation using idw
###
###
cov.bind = function(cov.var, cov.list) {
  names(cov.list) = 1:length(cov.list)
  col = c("X", "Y", cov.var)
  
  df = data.frame()
  for (i in 1:length(cov.list)) {
    cov.df = cov.list[[i]]
    if (cov.var %in% names(cov.df)) {
      df = rbind(df, cov.df[col])
      df = as.data.frame(df)
    }
  }
  df
}

plot.prediction = function(mask,
                           cov.var=NULL,
                           point.var=NULL,
                           nmax = 1,
                           maxdist = 10000,
                           cov.list=NULL,
                           point.df=NULL) {
  mask = as.data.frame(mask)
  df = mask[, c("x", "y")]
  coordinates(df) = ~ x + y
  gridded(df) = TRUE
  
  plot.list = list()
  output.df = mask[, c("x", "y")]
  if (!is.null(cov.var)) {
    formula = paste(cov.var, "~ 1")
    formula.num = 1
    for (i in cov.var) {
      df1 = cov.bind(i, cov.list)
      var.df = df1
      coordinates(var.df) = ~ X + Y
      
      if (is.numeric(var.df[[i]])) {
        # Inverse Distance Weighting
        idw = idw(
          formula = as.formula(formula[formula.num]),
          locations = var.df,
          newdata = df,
          nmax = as.numeric(nmax),
          maxdist = as.numeric(maxdist)
        )
        idw.output = as.data.frame(idw)
        names(idw.output)[1:3] = c("x", "y", "Prediction")
        output.df[[i]] = idw.output[,3]
        
        plot.list[[i]] = ggplot(data = idw.output, aes(x = x, y = y)) +
          geom_tile(data = idw.output, aes(fill = Prediction)) +
          scale_fill_gradient(low = "#FEEBE2", high = "#7A0177") +
          coord_equal() +
          ggtitle(i) 
        
        formula.num = formula.num + 1
        
      }
      else if (is.character(var.df[[i]])){
        nearest = nn2(df1[, c("X", "Y")],
                      mask[, c("x", "y")],
                      k = 1,
                      searchtype = "radius",
                      radius = maxdist)
        nearest$nn.idx[which(nearest$nn.idx == 0)] = NA
        output.df[[i]] = df1[,3][nearest$nn.idx]
        
        plot.df = cbind(mask[, c("x", "y")], df1[,3][nearest$nn.idx])
        names(plot.df) = c("x", "y", "Prediction")
        
        plot.list[[i]] = ggplot(data = plot.df, aes(x = x, y = y)) +
          geom_tile(aes(fill = Prediction)) +
          coord_equal() +
          ggtitle(i)
        
        formula.num = formula.num + 1
        
      }
    }
  }
  
  if (!is.null(point.var)){
    var.num = 1
    for (i in point.var){
      point = point.df[point.df$feature==i,c("X","Y")]
      nearest = nn2(point, 
                    mask[, c("x", "y")], 
                    k = 1)
      nearest$nn.dists[which(nearest$nn.idx == 0)] = NA
      output.df[[i]]=nearest$nn.dists
      
      plot.df = cbind(mask[, c("x", "y")],nearest$nn.dists)
      names(plot.df) = c("x", "y", "Distance_Prediction")
      
      
      plot.list[[i]] = ggplot(data = plot.df, aes(x = x, y = y)) +
        geom_tile(aes(fill = Distance_Prediction)) +
        geom_point(data = point, aes(x = X, y = Y), colour = "white") +
        coord_equal() +
        ggtitle(i)+
        coord_flip(ylim = c(min(mask[["y"]]),max(mask[["y"]])), 
                   xlim = c(min(mask[["x"]]), max(mask[["x"]])))
      var.num = var.num + 1
    }
  }
  
  names(output.df)[-1:-2] = c(cov.var,point.var)
  names(plot.list) = c(cov.var,point.var)
  multi.return = list(output.df, plot.list)
  
  return(multi.return)
}