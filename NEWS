acre 0.1.0
--------------

* Initial release of acre.

* Implemented model could fit the data as follows:

    * acoustic capture history stand along or with the following additional information:
    
        * "animal_ID" -- individual identification.

        * "toa" -- time-of-arrival.

        * "bearing" -- angle/direction.

        * "ss" -- signal strength.

        * "dist" -- distance.
    
    * extra covariates may affect the animal density or the probability of a call to be detected via the capture hisotry information:
    
        * "session" or time related covariates.
        
        * "trap" or detector related covariates.
        
        * "location" related covariates.
        

* options of the the model:

    * detection function -- the model could support detection functions of 'hn', 'hhn', 'hr', 'th', 'lth' and 'ss'.
    
    * specified start value of each coefficient.
    
    * specified boundary of each coefficient.
    
    * specified fixed value of each coefficient.
    
    * specified cut off threshold, link function for signal strength data.
    
    * specified link function and formula for the relationship between the capture history related coefficients and the extra covariates.
    
    * specified survey length for each session in a multi-session study.
    
    * specified sound speed.
    
    
* simulation:
    
    * simulate capture history from an output object of a fitted model.
    
    * simulate capture history from specified values of the coefficients.
    
* plot:
    
    * survey area with detectors locations.
    
    * capture history, stand alone or with additional information of "animal_ID", "ss", "bearing" or "dist".
    
    * location related covariates.
    
    * animal density surface from an output object of a fitted model. (through a different function instead of "plot()")
    
    

