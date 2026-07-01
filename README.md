# Spatiotemporal-control-of-charge-1-topological-defects-in-polar-active-matter
source code for the paper of the same name in PRR, DOI: https://doi.org/10.1103/xc83-kz86

## which OpenFPM version
used OpenFPM 4.1.0 (installation from Docker container on Linux) 
OpenFPM on GitHub: https://github.com/mosaic-group/openfpm
OpenFPM website: https://openfpm.mpi-cbg.de

## file structure
use:
- non_homogeneous_activity.cpp for all simulation without control, set parameters as in example file run_psib_m0_3_to_m0_45.pl (in paper: Fig.1c)
- non_homogeneous_activity_activity_control_rd.cpp for all simulation without control, set parameters as in example file run_rd_8_act_m100_to_m40.pl (in paper: Fig. 3a,b)
- non_homogeneous_activity_hanangata.cpp for all simulation without control, set parameters as in example file run_hanangata.pl (in paper: Fig.4d)
- Makefile to compile (change first line depending on where you saved the example.mk file from OpenFPM, see documentation there)

