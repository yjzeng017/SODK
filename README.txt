# 0) SUPPLEMENTS:
The supplements consist of code, data and additional appendices (.pdf files). 
For a description of the code and data content and how to use it see below.

# 1) SETUP: all the numerical studies are based on R and Python environments, thus readers should install both R and python platforms
	a) For R: Standard libraries should be available on CRAN. The following versions of the packages and R were used:
		- R version 4.4.1
		- reticulate 1.38.0 (interface to 'Python' modules, classes and functions.)
		- georob 0.3-19 (for GPRV method)
		- hetGP 1.1.6 (for GPH method)
	
	b) For Python: Standard package should be installed by using pip. The following versions of the packages and Python were used:
		- Python version 3.9.18
		- gpy 1.13.2 (for GPST method)
		- gpflow 2.9.2 (for SODK methods)
	In addtion, these packages depends on other packakes, such as numpy, tensorflow and so on. 
	These dependencies will be automatically installed when you install ''gpy'' or ''gpflow''.
	
	c) more R packages for numerical studies:
		- ggplot2, plot3D, plotly, geomtextpath, patchwork, viridis, (for visulization)
		- lhs (for Latin Hypercube Samples)
		- doSNOW, foreach (for parallel computing)
		- GeoModels, rrcov  (for bivariate Matern model and EH swapping, respectively)
		- gstat (for Jura dataset)
		- dplyr	
	
# 2) SETUP code:
The root folder "SODK" contains all codes and data.

######################################################################################
NOTE: 	a) When run the codes, the user should specifiy the root path to "SODK" folder for necessary
	b) The user should also specify the path of Python environment.
######################################################################################

# 3) FILE INFORMATION code:

SODK method(./SODK):
	- SODK.R: contain main function for SODK method
	- penalty-plot.R: produce the plots for shrinkage effects of normal-gamma prior (Figure 1).

Method(./SODK/benchmark/method):
	- comprison-methods.R: used for GPRV, GPH and GPST methods in main text, and additional GPST-MCMC, SODK-MCMC for runtime analysis in Supplementary
	- GPST-MCMC.py and ODK-MCMC.py: python scripts for GPST-MCMC, SODK-MCMC methods. 
				When load the script "runtime-methods.R", these Python scripts will be also loaded.

Simulation(./SODK/benchmark/experiment):
	- simulation 1: sim1.R (Figure 2)
	- simulation 2: 
		a) sim2-helping-functions.R: various helping functions for simulation 2, evaluation, contamination, ...
		b) sim2-outlier-detection.R: used for outlier detection performance of various methods (Figure 3)
		c) sim2-prediction.R: used for spatial prediction performance of various methods (Table 1)
Real example(./SODK/benchmark/experiment):
	- real-jura.R: used for investigating Jura data set (Figure 4(a) and Figure 4(c))
	- real-airfoil.R: used for investigating Airfoil simulation example (Figure 4(b) and Table 2)
	- Airfoil.txt: Airfoil simulation data


Supplementary(./SODK/benchmark/supp):
	- convergence.R: produce convergence plots for ODK method in various simulation scenarios in Supplementary A.1
	- parameter-sensitivity-analysis.R: used for sensitivity analysis of the regularization paramter omega in ODK method in Supplementary A.2
	- threshold-analysis.R: used for thresholding analysis in  Supplementary A.3
	- runtime-analysis.R, runtime-analysis-MCMC.R: used for comprison of runtime in Supplementary A.4
	- more-numerical-results (./SODK/benchmark/supp/more-numerical-results):
			-- spatial.R: an additional spatial example
			-- Solid-End-Milling.R: an additional computer experiment example: Solid end milling
			-- transcriptomic-patterns.R: identify the spatial patterns of transcriptomic examples
			-- transcriptomic-analysis.R: analysis the spatial data from ' transcriptomic-patterns.R'.

Results: the results are saved as .Rdata and figure.
	- ./SODK/result: contains .Rdata files.
	- ./SODK/figure: contains figures
