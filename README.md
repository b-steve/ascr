# admbsecr

Traditional capture-recapture approaches to estimating animal abundance or density ignore an obvious spatial component of capture probability; organisms close to traps are more likely to be captured than those that are far away. Explicitly accounting for an individual's location provides additional information from which to infer animal density. Spatially explicit capture-recapture (SECR) methods have been developed for this purpose. An advantage of these over traditional capture-recapture methodology is that they allow for animal density estimation using passive detectors (e.g., cameras or microphones) over a single sampling occasion.

In the simplest case, distances between traps provide the spatial information required to implement SECR methods. In some situations, passive detectors provide supplementary information which can be used to better estimate the exact location of an individual. This could be the precise time of arrival and/or received strength of an acoustic signal, the estimated angle and/or distance between an animal and the trap, or even the exact location of the animal itself. Currently available software implementations of SECR methods are unable make use of such information.

[AD Model Builder] (http://admb-project.org/) (ADMB) is a statisical software package most widely used for nonlinear modelling, and appears to be well suited to the implementation of maximum likelihood SECR methods. Although growing in popularity since becoming freely available, open-source software in 2008, ADMB is used by a minority of statisticians and ecologists, who, in general, are far more comfortable with the popular programming language and software environment [R](http://www.r-project.org).

The aim of admbsecr is to bridge both of these gaps. Using the R function `admbsecr()`, a user is able to fit SECR models that incorporate additional spatial information. This calls ADMB to fit the model and return the results to the R session.

## Installation

* Run the script found [here](https://raw.githubusercontent.com/b-steve/admbsecr/master/inst/scripts/install.r).

* Optimisation is carried out using an ADMB executable. To ensure that the executable is running correctly, run
```
test.admbsecr(quick = TRUE)
```

* If you wish, for more rigorous testing, run
```
test.admbsecr(quick = FALSE)
```
This can take a little while.

* This package has not been widely tested on all platforms. If there are any errors with the installation or in the above tests please get in touch using the contact details below.

* Note that this R package does not require an ADMB installation.

## Contact

Ben Stevenson, bcs5[at]st-andrews.ac.uk

## Acknowledgements

I am grateful to David Borchers and Hans Skaug for their continued assistance with this project.

Additional thanks to the National Geographic Society/Waitt Grants Program (Grant #W184-11).
