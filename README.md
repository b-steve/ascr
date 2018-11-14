# ascr

Acoustic surveys provide a cheap, efficient means of monitoring wildlife populations. Animal density is one demographic parameter that is often of particular interest---however, it is difficult to estimate from such data for various reasons. For example,

1. Not every individual within the region of interest is detected, and so those that are not detected must somehow be accounted for.

2. The source locations of the detected cues are not observed, and so the extent of the survey area that is monitored by the acoustic detectors is not known.

Spatial capture-recapture (SCR) models appear well-suited to overcoming these challenges. SCR estimates detectability across the survey area, and so it is possible to estimate the density of cues that have not been detected. Furthermore, SCR accounts for the cues' unknown source locations by modelling these as latent variables. However, data from acoustic surveys are a little different to "typical" SCR data (for example, that use cameras or physical traps to detect individuals):

1. It is often cues (rather than individuals) that are the unit of detection.

2. Cues cannot necessarily be attributed to individuals, so it is not always known how many unique animals have been detected.

3. Acoustic detectors often collect additional data such as the received signal strength and time of arrival of detected cues, and estimated distances or bearings from the detector to detected cues' source locations. Such data are informative about the locations of detected animals, and can thus improve inference.

SCR methodology has been developed to account for the these data features (see References section, below). The ascr ("Acoustic SCR") package provides a software implementation of these to facilitate their use. Parameter estimation is carried out by an executable compiled by the powerful optimisation software package AD Model Builder (ADMB)---but users are not required to have any familiarity with this, or even have it installed. However, as the ascr package contains an ADMB executable it is not permitted on the Comprehensive R Archive Network (CRAN). See below for installation instructions.

The ascr package was formerly named admbsecr (https://github.com/b-steve/admbsecr). Any code written for use with admbsecr should also run cleanly with ascr. Any functions that were named differently in admbsecr have been aliased with the new corresponding functions in ascr.

## Installation

* Run the script found [here](https://raw.githubusercontent.com/b-steve/ascr/master/inst/scripts/install.r).

* Load the package by running
```
library(ascr)
```

* Optimisation is carried out using an ADMB executable. To ensure that the executable is running correctly, run
```
test.ascr(quick = TRUE)
```

* If you wish, for more rigorous testing, run
```
test.ascr(quick = FALSE)
```
This can take a little while.

* This package has not been widely tested on all platforms. If there are any errors with the installation or in the above tests please get in touch using the contact details below.

## References

Almost all methods presented in the references below have been implemented in the ascr package:

* Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas, L., and Marques, T. A. (2015). A unifying model for capture-recapture and distance sampling surveys of wildlife populations. *Journal of the American Statistical Association*, *110*, 195--204.

* Dawson, D. K., and Efford, M. G. (2009). Bird population density estimated from acoustic signals. *Journal of Applied Ecology*, *46*, 1201--1209.

* Efford, M. G., Dawson, D. K., and Borchers, D. L. (2009). Population density estimated from locations of individuals on a passive detector array. *Ecology*, *90*, 2676--2682.

* Kidney, D., Rawson, B. M., Borchers, D. L., Stevenson, B. C., Marques, T. A., and Thomas, L. (2016). An efficient acoustic density estimation method with human detectors applied to gibbons in Cambodia. *PLoS ONE*, *11*, e0155066.

* Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift, R. J., Gillespie, D. M., and Measey, G. J. (2015). A general framework for animal density estimation from acoustic detections across a fixed microphone array. *Methods in Ecology and Evolution*, *6*, 38--48.

The ascr package also implements methods that are not currently published, but can be found in my PhD thesis:

* Stevenson, B. C. (2016). *Methods in spatially explicit capture-recapture* (PhD thesis).

## Contact

Ben Stevenson, ben.stevenson[at]auckland.ac.nz

## Acknowledgements

Thanks to the EPSRC and the NERC (Grant #EP/1000917/1), and the National Geographic Society/Waitt Grants Program (Grant #W184-11).
