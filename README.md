# ascr

 Spatial capture-recapture (SCR) models estimate animal density from records of when and where individuals were detected. The first SCR methods for acoustic data were developed by Efford, Dawson, and Borchers (2009).

 The ascr package fits SCR models to estimate animal density from acoustic surveys. In particular, these models can incorporate additional data like received signal strengths, times of arrival, estimated distances, and estimated bearings, which are informative about animals' locations. Other more general software packages exist for SCR (e.g., the comprehensive secr R package), but they cannot handle all data types listed above.

 Models that can be fitted in `ascr` include those described by Efford, Dawson, and Borchers (2009), Borchers et al (2015), and Stevenson et al (2015). The package also implements some unpublished methods that are described in Stevenson (2016), a PhD thesis.

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

If there are any errors with the installation or in the above tests
please get in touch using the contact details below.

NOTE: OS X users will create the ADMB executable themselves:

1. Download and install [ADMB](http://www.admb-project.org/downloads/).
2. In R, run `system.file(package = "ascr")` to find the directory where `ascr` is installed.
3. Navigate to this directory in a terminal, then to the subirectory `ADMB/src`.
4. Run `admb secr`.
5. This will create an executable named `secr`. Move it to the subdirectory `ADMB/bin/mac`.
6. Test that the executable works using the `test.ascr()` code above.

## References

* Borchers, D. L., Stevenson, B. C., Kidney, D., Thomas, L., and Marques, T. A. (2015). A unifying model for capture-recapture and distance sampling surveys of wildlife populations. *Journal of the American Statistical Association*, *110*: 195--204.

* Efford, M. G., Dawson, D. K., and Borchers, D. L. (2009). Population density estimated from locations of individuals on a passive detector array. *Ecology*, *90*: 2676--2682.

* Stevenson, B. C., Borchers, D. L., Altwegg, R., Swift, R. J., Gillespie, D. M., and Measey, G. J. (2015). A general framework for animal density estimation from acoustic detections across a fixed microphone array. *Methods in Ecology and Evolution*, *6*: 38--48.

* Stevenson, B. C. (2016). *Methods in spatially explicit capture-recapture* (PhD thesis).

## Contact

Ben Stevenson, ben.stevenson[at]auckland.ac.nz

