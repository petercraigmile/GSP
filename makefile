
RFILES		=	Rcpp/GSP_covariance_functions.R \
			Rcpp/GSP_site_functions.R \
			Rcpp/GSP_simulation.R \
			Rcpp/GSP_likelihood_methods.R \
			Rcpp/GSP_prediction.R

CPPFILES	=	Rcpp/GSP.cpp

RNOCPP		=	R_only/R_only_GSP_covariance_functions.R \
			R_only/R_only_GSP_site_functions.R \
			Rcpp/GSP_simulation.R \
			Rcpp/GSP_likelihood_methods.R \
			Rcpp/GSP_prediction.R


build:		$(RFILES)
		mkdir releases/`date "+%Y_%m_%d"`
		cat $(RFILES) >  releases/`date "+%Y_%m_%d"`/GSP_`date "+%Y_%m_%d"`.R
		cat $(CPPFILES) > releases/`date "+%Y_%m_%d"`/GSP_`date "+%Y_%m_%d"`.cpp
		cat $(RNOCPP) > releases/`date "+%Y_%m_%d"`/GSP_R_only_`date "+%Y_%m_%d"`.R

