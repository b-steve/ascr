#' Plot an ascrdata object
#'
#' \description{
#' Generate several types of plots from an ascrdata object including: 
#'   
#'   1. trap locations and capture histories using \link{show.capt}.
#'   
#'   2. mask and trap layout using \link{show.survey}.
#'         
#'   3. the predicted covariate values at each mask point using 
#'      \link{show.prediction}. See 'Details' for further information.
#' }
#' 
#' When \code{covariate} is chosen for \code{type}, \code{datalist}
#' should at least contain either \code{cov.var} or \code{point.var} and the
#' corresponding \code{cov.list} or \code{point.df}. See \link{alldata} for details.
#' 
#' @param datalist A list containing capture histories, traps locations, etc.
#'     It is most easily created using \link{alldata}. 
#' @param type A character string specifying the type of plot to be generated.
#'     One of "capt", "survey" or "covariate".
#' @param session The session(s) for which the mask point and trap
#'     locations are to be plotted. If \code{NULL}, the default, then plots will 
#'     be shown for all sessions.
#' @param cov A vector of character string specifying the covariates to be plotted. 
#'     Only used when \code{covariate} is chosen for \code{type}.
#'     By default, "all" will plot all covariates stored in \code{datalist}.
#'
#' @export

plot.ascrdata = function(datalist,type = "covariate", session = NULL, cov = "all"){
  if (type == "capt") return(show.capt(capt = datalist$capt, traps = datalist$traps, mask = datalist$mask, session = session))
  
  if (type == "survey"){
    if (is.null(session)){session = "all"}
    return(show.survey(traps = datalist$traps, mask = datalist$mask, session = session))}
  
  if (type == "covariate") {
    ask.save <- par("ask")
    par(ask = TRUE)
    ## Making sure par is restored on function exit.
    on.exit(par(ask = ask.save))
    
    ##plotting all covariates
    if (cov == "all"){
      n.cov = length(datalist$cov.var)+length(datalist$point.var)
      ##specifying session
      if (!is.null(session)){
        for (i in 1:n.cov){
          print(datalist$cov.prediction.plot[[session]][[i]]+labs(subtitle = paste("session", session)))
        }
      } else{ ##plotting all sessions
        n.session = length(datalist$traps)
        for (i in 1:n.session){
          for (j in 1:n.cov){
            print(datalist$cov.prediction.plot[[i]][[j]]+labs(subtitle = paste("session", i)))
          }
        }
      }
    } else{ ##specifying covariates
      if (!is.null(session)){
        for (i in 1:length(cov)){
          print(datalist$cov.prediction.plot[[session]][[cov[i]]]+labs(subtitle = paste("session", session)))
        }
      } else{ ##plotting all sessions
        n.session = length(datalist$traps)
        for (i in 1:n.session){
          for (j in 1:length(cov)){
            print(datalist$cov.prediction.plot[[i]][[cov[j]]]+labs(subtitle = paste("session", i)))
          }
        }
      }
    }
  }
}

