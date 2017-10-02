# rmetasim  population genetic simulation engine

This software started as a c++ library (Strand, 2002) and has been modified to be a c++ backend and R frontend.  The idea is to provide a flexible environment for designing and running forward-time genetic simulations.  Populations can have complex demography and there can be arbitrary dispersal among them.  Microsatellites, Infinite alleles and sequences can be simulated.

# important recent changes (March 2016):

Several functions (landscape.amova.locus, landscape.locus,
	landscape.locus.states, landscape.mismatchdist,
	landscape.states) that access genetic information from
	landscapes had a strange order for parameters.  That has been
	resolved to be the same as the rest of rmetasim, where the
	first argument is the landscape.  If you have code that uses
	these and is based on rmetasim v<3.0, you will need to
	rearrange the parameters for these function calls, unless
	parameters were explicitly named in the calls.


Strand, A.  Metasim 1.0: an individual-based environment for
    simulating population genetics of complex population
    dynamics. Mol. Ecol. Notes, 2:373-376, 2002
