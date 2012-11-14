\name{sim.capthist.ss}
\alias{sim.capthist.ss}
\title{Simulated signal strength capture history matrix}
\usage{
  sim.capthist.ss(traps, popn, detectpars)
}
\arguments{
  \item{traps}{trap locations.}

  \item{popn}{simulated population.}

  \item{detectpars}{detection function parameters.}
}
\description{
  Simulating a signal strength capture history matrix.
  Signal strength detection function uses a log link
  function, and thus is different to \code{sim.capthist}
  where \code{detectfn = 10}.
}

