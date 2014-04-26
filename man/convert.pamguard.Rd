\name{convert.pamguard}
\alias{convert.pamguard}
\title{Create a capture history object from a PAMGuard output file}
\usage{
convert.pamguard(dets, mics, sound.speed = 330)
}
\arguments{
  \item{dets}{Detection output dataframe from PAMGuard.}

  \item{mics}{A matrix containing the coordinates of
  microphone locations.}

  \item{sound.speed}{The speed of sound in metres per
  second.}
}
\description{
Converts a PAMGuard output file to a capture history object
suitable for use with the \link{admbsecr} function. This
uses \link{make.acoustic.captures} to allocate call
identities to detections.
}

