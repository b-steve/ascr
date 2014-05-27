# admbsecr

Traditional capture-recapture approaches to estimating animal abundance or density ignore an obvious spatial component of capture probability; organisms close to traps are more likely to be captured than those that are far away. Explicitly accounting for an individual's location provides additional information from which to infer animal density. Spatially explicit capture-recapture (SECR) methods have been developed for this purpose. An advantage of these over traditional capture-recapture methodology is that they allow for animal density estimation using passive detectors (e.g., cameras or microphones) over a single sampling occasion.

In the simplest case, distances between traps provide the spatial information required to implement SECR methods. In some situations, passive detectors provide supplementary information which can be used to better estimate the exact location of an individual. This could be the precise time of arrival and/or received strength of an acoustic signal, the estimated angle and/or distance between an animal and the trap, or even the exact location of the animal itself. Currently available software implementations of SECR methods are unable make use of such information.

[AD Model Builder] (http://admb-project.org/) (ADMB) is a statisical software package most widely used for nonlinear modelling, and appears to be well suited to the implementation of maximum likelihood SECR methods. Although growing in popularity since becoming freely available, open-source software in 2008, ADMB is used by a minority of statisticians and ecologists, who, in general, are far more comfortable with the popular programming language and software environment [R](http://www.r-project.org).

The aim of admbsecr is to bridge both of these gaps. Using the R function `admbsecr()`, a user is able to fit SECR models that incorporate additional spatial information. This calls ADMB to fit the model and return the results to the R session.

## Installation

To install:

* Package dependencies will not install automatically as admbsecr is not on CRAN. To install these, run:
```
install.packages(c("CircStats", "lattice", "matrixStats", "plyr", "Rcpp", "R2admb", "secr"))
```

* For the [stable version on R-Forge](https://r-forge.r-project.org/projects/admbsecr/):
```
install.packages("admbsecr", repos = "http://R-Forge.R-project.org")
```
This requires the newest version of R.

If you receive the warning `package 'admbsecr' is not available (for R version 3.x.x)`, you may have attempted to install during R-Forge's daily build process. Check [here](https://r-forge.r-project.org/R/?group_id=1506) to see the build status of admbsecr and wait for a "Current" status.

* For the development version (i.e., this repository; stability not guaranteed):
```
library(devtools)
install_github("admbsecr", "b-steve")
```
For this option Windows users will need a compatible version of [Rtools](http://cran.r-project.org/bin/windows/Rtools/) installed.

* Note that this R package does not require an ADMB installation.

## Troubleshooting

If you are experiencing difficulties with getting things running, please get in touch. We are attempting to iron out any problems with the installation procedure. Below are some known issues and some possible fixes.

* It is best to use the newest version of R. In particular, installing from R-forge (as above) will not work unless R is up to date.

* There is a known issue with some versions of OSX, which causes a segmentation fault when the `admbsecr()` function is called.

* If a permission error is thrown when the `admbsecr()` function is called, it probably means you do not have permission to call the ADMB executable used to fit the model. In R, run `system.file(package = "admbsecr")` to find the directory in which the admbsecr package is installed, then navigate to the appropriate directory (according to your OS) in `ADMB/bin`. Allow execute permission to the executable `secr` therein; e.g., using `chmod u+x secr` at a shell command line on Unix or OSX.

## Contact

Ben Stevenson, bcs5@st-andrews.ac.uk

## Acknowledgements

I am grateful to David Borchers and Hans Skaug for their continued assistance with this project.

Additional thanks to the National Geographic Society/Waitt Grants Program (Grant #W184-11).
