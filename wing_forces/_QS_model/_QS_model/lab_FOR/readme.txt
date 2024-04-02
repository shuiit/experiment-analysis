
wing_FOR (frame of reference)
-----------------------------

aerodynamicForceSimpleQS.m
calculates simple aerodynamic force using coeffs from whitney and wood (i.e. not Dickinson). 
uses wing speed u in the WING frame of reference


lab_FOR
-------
the function quasiSteadyForce.m is the master function
it calculates QS force using one of three QS models
I used the "simple" one (modelType=2), which is the same as in the function in wing_FOR.
