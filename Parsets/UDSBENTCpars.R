

mxval <- 17 ## maximum value for the outcome
mnval <- 0 ## minimum value for the outcome

incage <- FALSE  ## whether curve increases with respect to NACCAGE
inceduc <- as.logical(1-incage)  ## whether curve increases with respect to EDUC 
## default to opposite of direction with NACCAGE

chopT <- TRUE ## Whether or not any data chopping is required
chopvals <- function(x){
  x <- x[!(x$NACCAGE == 40 & x[ ,varset[indx]] < 12), ]
  x <- x[x[ ,varset[indx]]>5, ]
  x
} ## the stuff behind the ! (not) sign indicates 
## what will be chopped for fitting of the SD curve

weightSD <- FALSE ## Whether to weight SD fit in scam2 according to N in the window
linweightSD <- FALSE ## Whether to have linear weights for NACCAGE in scam1

linEd <- TRUE ## Whether to restrict education effect to linear in SCAM model

incSD <- TRUE ## Whether to have SD as increasing with age

gamfitSD <- FALSE ## If TRUE, then don't restrict to monotonic and perform regular GAM fit
