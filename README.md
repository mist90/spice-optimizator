# spice-optimizator
Optimizator of Spice-parameters for NMOS Level 3 Model 

## Instruction:
* Set up ini configuration file.
* Run command:
python3 main.py path-to-ini-file
* If optimization will end successfully then program print out new values of spice-parameters:
VT0, NFS, KP, THETA, KAPPA
and plot graphs on graphical window.

## Configuration ini-file.
* Set up your experimental dataset of points on VA-characteristics - see examples in configs directory.
* Set up start values of spice parameters - see examples in configs directory.
* Full list of spice-parameters: VT0, NSUB, NFS, Rd, Rs, KP, U0, Weff, Leff, TOX, VMAX, THETA,
XJ, PHI, DELTA, GAMMA, ETA, KAPPA.
* For parameters to be optimized you may set limit min and max values:
VT0min, VT0max, NFSmin, NFSmax, KPmin, KPmax, THETAmin, THETAmax, KAPPAmin, KAPPAmax.

Dependency: Python3, numpy, scipy, matplotlib
