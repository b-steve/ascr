\name{sim.capthist.dist}
\alias{sim.capthist.dist}
\title{Simulated distance error model capture history matrix}
\usage{
  sim.capthist.dist(traps, popn, detectpars)
}
\arguments{
  \item{traps}{trap locations.}

  \item{popn}{simulated population.}

  \item{detectpars}{detection function parameters.}
}
\description{
  Simulating a capture history matrix for a distance error
  model. This function specifically deals with the case
  where there are two traps in the same location, each of
  which has different detection function parameters.
}

