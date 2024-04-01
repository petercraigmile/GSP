# GSP (Geostatistical stationary processes)

R and R/C++ code for simulating, fitting, and predicting geostatistical
stationary processes.  At the moment, the code assumes that each
process is isotropic (the covariance depends on the distance between
spatial locations and not the angle).


Further remarks:

- There are two versions: (i) an R/C++ version that uses the Rcpp
  [http://www.rcpp.org] library to call C++ functions from R; (ii) a R
  only version.  The R/C++ version is quicker.  The 'Rcpp' folder
  contains the R and C++ functions for the R/C++ version.  The
  'R_only' folder contains the R functions that differ from the R/C++.

- Releases contains complete versions of the code, organized by date.
The releases are created using the makefile.

- The parameterization of the covariance function 'GSP.Matern' does
  not match with 'GSP.exp', 'GSP.Diggle', and 'GSP.Gaussian'

- This code could be optimized further.

- Most R and C++ functions need further documentation.



