#' Plotting mask and trap layout
#'
#' Plots the mask points and trap locations used in a model fitted
#' with the function \link{admbsecr}.
#'
#' @param ... Further arguments to be passed to \link{plot}.
#' @inheritParams locations
#'
#' @examples
#' show.survey(example$fits$simple.hn)
#'
#' @export
show.survey <- function(fit, ...){
    plot(fit$args$mask, pch = ".", cex = 3, asp = 1, ...)
    points(fit$args$traps, pch = 16, col = "red")
}

#' Plotting the detection probability surface
#'
#' Plots the detection probability surface, based on trap locations
#' and the estimated detection function from a model fitted using
#' \link{admbsecr}.
#'
#' @inheritParams locations
#'
#' @examples
#' show.detsurf(example$fits$simple.hn)
#'
#' @export
show.detsurf <- function(fit){
    p.det <- p.dot(fit)
    mask <- get.mask(fit)
    unique.x <- sort(unique(mask[, 1]))
    unique.y <- sort(unique(mask[, 2]))
    z <- matrix(NA, nrow = length(unique.x), ncol = length(unique.y))
    n.mask <- nrow(mask)
    for (i in 1:n.mask){
        x <- mask[i, 1]
        y <- mask[i, 2]
        index.x <- which(x == unique.x)
        index.y <- which(y == unique.y)
        z[index.x, index.y] <- p.det[i]
    }
    print(wireframe(z, zlim = c(0, 1),
                    zlab = list("Detection probability", rot = 95),
                    xlab = "x", ylab = "y", shade = TRUE))
    invisible(TRUE)
}
