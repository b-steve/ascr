# admbsecr

Traditional capture-recapture approaches to estimating animal abundance or density ignore an obvious spatial component of capture probability; organisms close to traps are more likely to be captured than those that are far away. Explicitly accounting for an individual's location provides additional information from which to infer animal density. Spatially explicit capture-recapture (SECR) methods have been developed for this purpose. An advantage of these over traditional capture-recapture methodology is that they allow for animal density estimation using passive detectors (e.g., cameras or microphones) over a single sampling occasion.

In the simplest case, distances between traps provide the spatial information required to implement SECR methods. In some situations, passive detectors provide supplementary information which can be used to better estimate the exact location of an individual. This could be the precise time of arrival and/or received strength of an acoustic signal, the estimated angle and/or distance between an animal and the trap, or even the exact location of the animal itself. Currently available software implementations of SECR methods are unable make use of such information.

AD Model Builder (ADMB) is a statisical software package most widely used for nonlinear modelling, and appears to be well suited to the implementation of maximum likelihood SECR methods. Although growing in popularity since becoming freely available, open-source software in 2008, ADMB is used by a minority of statisticians and ecologists, who, in general, are far more comfortable with the popular programming language and software environment R.

The aim of admbsecr is to bridge both of these gaps. Using the R function `admbsecr()`, a user is able to fit SECR models that incorporate additional spatial information. This calls ADMB (through use of the package R2admb) to fit the model and return the results to the R session.

## Installation

To install:

1. Download and install [AD Model Builder](http://admb-project.org/).

2. Install the admbsecr R library.

* For the [development version on R-Forge](https://r-forge.r-project.org/projects/admbsecr/):
```
install.packages("admbsecr", repos = "http://R-Forge.R-project.org")
```

* For the very latest version (i.e., this repository):
```
library(devtools)
install_github("admbsecr", "b-steve")
```

NOTE: Due to what I assume is a bug with R-Forge, the `install.packages()` command above seems to be installing a very old version of admbsecr. A workaround is to click [here](http://r-forge.r-project.org/R/?group_id=1506), download the appropriate file (determined by your operating system), and unzip it to your R package library directory.

## Troubleshooting

* Ensure the command `admb` is recognised at the terminal command line. Mac/Unix users may have to add the following to ~/.bashrc, ~/.bash_login, or a similar script that runs on startup.
```bash
#!/bin/bash
export ADMB_HOME=/path/to/admb/
PATH=$PATH:$ADMB_HOME/bin
export PATH
```

* Windows users may need to add the environment variable `R_SHELL` with a value that points to `cmd.exe`, which is found in the System32 folder.

* An error while calling `install.packages()` (as above) is probably due to an outdated version of R. Update and try again.

## Acknowledgements

I am grateful to David Borchers and Hans Skaug for their continued assistance with this project.

Additional thanks to the National Geographic Society/Waitt Grants Program (Grant #W184-11).
