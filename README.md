# ggl
This package contains codes to estimate and model galaxy-galaxy lensing 
statistics between source shape and lens catalogues.

The codes presented are:

* delsigest.f -- Fortran code to estimate Delta Sigma(R) from shape and 
lens catalogues (no external libraries required).

* lensmod.py -- Python code to model a set of lensing and clustering 
statistics (needs Python libraries: numpy, scipy, astropy)

The other accompanying files are:

* sources_buzzard_pix128_zp0pt80_1pt00.dat.gz -- shape sub-sample from 
the Buzzard mocks.

* lenses_buzzard_lrg_pix128_zs0pt40_0pt60.dat.gz -- lens LRG sub-sample 
from the Buzzard mocks.

* pz_sources_buzzard_zp0pt80_1pt00.dat -- N(z) of shape sub-sample.

* pz_lenses_buzzard_lrg_zs0pt40_0pt60.dat -- N(z) of lens sub-sample.

* pkcambhalofit_zrange_buzzard.dat -- CAMB halofit model power spectrum 
generated over a set of redshifts, with the same fiducial cosmology as 
the Buzzard mocks.

* delsigest.dat -- output of delsigest.f.

* delsigmod.dat -- Delta Sigma(R) output of lensmod.py with b=1.
