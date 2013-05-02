\name{sim.capthist.ss}
\alias{sim.capthist.ss}
\title{Simulated signal strength capture history matrix}
\usage{
  sim.capthist.ss(traps, popn, detectpars, log.link,
    re = FALSE)
}
\arguments{
  \item{traps}{trap locations.}

  \item{popn}{simulated population.}

  \item{detectpars}{detection function parameters.}

  \item{log.link}{logical, if \code{TRUE} a log link
  function is used for expected signal strengths.}

  \item{re}{logical, if \code{TRUE} each sound is generated
  with a different source strength.}
}
\description{
  Simulating a signal strength capture history matrix.
}

