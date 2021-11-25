#' Create an ascr.data object
#' 
#' Sorting the data used for functions in R package \emph{ascr}.
#' 
#' Create a mask object to use with the function \code{\link{fit.ascr}} according 
#' to the \code{traps} and other argument provided. 
#' 
#' Also, create a list of data frame of the predicted covariate value at each
#' mask point to use with the function \code{\link{fit.ascr}} when fitting on 
#' inhomogeneous density surfaces. The function \code{\link{show.prediction}} is
#' used.
#' 
#' @param capt A list containing capture histories and supplementary data for 
#' each detected individual. It is most easily created using \code{\link{create.capt}}.
#' 
#' @param traps A matrix with two columns. Each row provides Cartesian 
#' coordinates for the location of a detector. Alternatively, this can be a 
#' list of such matrices if detections from multiple detector arrays 
#' (or 'sessions') are being used to fit a single model.
#' 
#' @param mask A matrix with two columns, or a list of such matrices for a 
#' multi-session model. Each row provides Cartesian coordinates for the location
#' of a mask point. It is most easily created using \code{\link{create.mask}}.
#' If missing, mask will be created using the provided \code{traps}, 
#' \code{capt} and \code{buffer}.
#' 
#' @param buffer The minimum distance between trap locations and the edge of the 
#' generated mask. This is used with \code{\link{create.mask}} to create the mask
#' point. Only applied when \code{mask} is missing.
#' 
#' @param cov.var A vector of character string of the covariate variables used 
#' to fit the model using \code{\link{fit.ascr}}. 
#' 
#' @param point.var A vector of the character string of the distance covariate
#' used to fit the model using \code{\link{fit.ascr}}.
#' 
#' @param cov.list A list of all the covariate data frames. Each row 
#' provides Cartesian coordinates named "X" and "Y" and its 
#' covariate values.
#' 
#' @param point.df A data frames of the distance covariates. Each row 
#' provides Cartesian coordinates named "X" and "Y" and the 
#' corresponding feature named "feature".
#' 
#' @param nmax the number of nearest observations that should be used to predict
#' the covariate value. This is used with \code{\link{show.prediction}}.
#' Only applied when \code{cov.var} is numeric variable.
#' 
#' @param maxdist only observations within a distance of maxdist from the 
#' prediction location are used to predict the covariate value. This is used 
#' with \code{\link{show.prediction}}.only applied when \code{cov.var} is numeric
#' variable.
#' 
#' @param ... Further arguments to be passed to create the mask point using 
#' \code{\link{create.mask}}.
#' 
#' @return An object of class \code{ascr.data}. The object is a list with 
#' elements \code{capt}, \code{traps}, \code{mask}, \code{cov.prediction.df}, etc. 
#' 
#' 
#' @export

alldata <- function(capt=NULL, traps=NULL, mask=NULL, buffer=NULL,cov.var=NULL,point.var=NULL,cov.list=NULL,point.df=NULL,nmax=10,maxdist=1000,...){
  datalist=list(capt=capt,
                traps=traps,
                mask=mask,
                cov.list=cov.list,
                point.df=point.df,
                cov.var=cov.var,
                point.var=point.var)
  if (is.null(mask)){datalist$mask = create.mask(traps = datalist$traps, buffer = buffer, ...)}
  
  if (!is.null(cov.var)&!is.null(cov.list) | !is.null(point.var)&!is.null(point.df)){
    if (class(datalist$mask) == "list"){
      prediction = lapply(datalist$mask,function(x){show.prediction(mask=x, datalist = datalist,nmax = nmax, maxdist = maxdist)})
      datalist$cov.prediction.df = lapply(prediction, function(x){x[[1]]})
      datalist$cov.prediction.plot = lapply(prediction, function(x){x[[2]]})
    } else {
      prediction = show.prediction(datalist$mask, datalist = datalist, nmax = nmax, maxdist = maxdist)
      datalist$cov.prediction.df = prediction[[1]]
      datalist$cov.prediction.plot = prediction[[2]]
    }
  }
  
  class(datalist) <- "ascrdata"
  return(datalist)
}
