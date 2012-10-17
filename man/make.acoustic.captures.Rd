\name{make.acoustic.captures}
\alias{make.acoustic.captures}
\title{Assigning ID numbers to sounds.}
\usage{
  make.acoustic.captures(mics, clicks, dt)
}
\arguments{
  \item{mics}{a matrix containing the coordinates of trap
  locations.}

  \item{clicks}{a data frame containing (at least): (i)
  \code{$tim$}, the precise time of arrival of the received
  sound, and (ii) \code{$trap} the trap at which the sound
  was recorded.}

  \item{dt}{a \code{K} by \code{K} matrix (where \code{K}
  is the number of traps) containing the time taken for
  sound to travel between each pair of traps.}
}
\value{
  A data frame. Specifically, the \code{clicks} dataframe,
  now with a new variable, \code{ID}.
}
\description{
  Identifies recaptures and assigns ID numbers to sounds
  recorded for an SECR model.
}
\details{
  Detected sounds are assumed to come from the same animal
  if times of arrival at different microphones are closer
  together than the time it would take for sound to travel
  between these microphones.
}
\author{
  David Borchers, Ben Stevenson
}

