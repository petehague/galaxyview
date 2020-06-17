# README #

The code is intended to extract velocity fields, rotation curves, and other information from APOSTLE galaxies

### Usage ###

Edit the first line of run.sh so that PARENTDIR points to the location you want the output to be written to

If there are any issues with compiling or running the code, replace g++ with icpc in the first line of the makefile

### Using with RainfallMCMC ###

To compile the MCMC modules, put the Rainfall folder in the root directory of this repository 