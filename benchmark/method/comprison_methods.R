#############################################################################
# The considered methods for comprison in the manuscript, including               
# - Gaussian process based on robust variogram via robust REML (GPRV)         
# - Gaussian process with heterogeneous noise (GPH)                            
# - Gaussian process with student's likelihood (Laplace approximation)(GPST)  
#############################################################################

#############################################################################
# The considered methods in runtime analysis (in Supplementary):   
# - Gaussian process based on robust variogram via robust REML (GPRV)         
# - Gaussian process with heterogeneous noise (GPH)        
# - Gaussian process with student's likelihood:                                
#   -- Laplace approximation (GPST-LA)                                       
#   -- Expectation Propagation approximation (GPST-EP)                        
#   -- MCMC (GPST-MCMC)
# 
# - Fully Bayesian SODK with MCMC sampling for estimation (SODK-MCMC)
#############################################################################

library(georob) # GPRV
library(hetGP) # for GPH
library(reticulate) # for interface to Python 

# NOTE: User should specify root path
root_path = "D:/SODK"
# NOTE: User should specify the path of Python environment
python_path = "D:/Anaconda3/python.exe" 
use_python(python_path)

# check the Python environment and the Python package 'GPy' for GPST
if(!py_module_available('GPy')){
  stop("Error message: the GPy module should be available. \n")
}

# check the Python environment and the Python package 'gpflow' for GPST-MCMC, SODK-MCMC
if(!py_module_available('gpflow')){
  stop("Error message: the GPy module should be available. \n")
}

# User should specify the path of the Python scripts: 
#   GPST-MCMC.py and SODK-MCMC.py

# load the Python scripts
source_python(paste(root_path, "/benchmark/method/GPST-MCMC.py", sep = '')) 
source_python(paste(root_path, "/benchmark/method/SODK-MCMC.py", sep = ''))



############################################################################
# GPRV in Kunsch et al. (2011)
############################################################################
# main function
GPRV = function(X, y, cov='Matern5_2', tuning.psi=2, nugget_estim=TRUE){
  # X: design matrix
  # y: response
  # cov: covariance function:
  #      - gauss: Gaussian kernel
  #      - exp: exponential kernel
  #      - Matern3_2: matern class with smooth parameter nu = 3/2
  #      - Matern5_2: matern class with smooth parameter nu = 5/2
  # nugget_estim: a logical scale defining whether to estimate the nugget (noise)
  # return: list of GPRV model and outlier
  
  n = nrow(as.matrix(X))
  d = ncol(as.matrix(X))
  if(cov == 'gauss') c
  if(cov == 'exp') variogram.model = 'RMexp'
  if(cov == 'Matern3_2')  { variogram.model = 'RMmatern'; nu= 3/2}
  if(cov == 'Matern5_2') { variogram.model = 'RMmatern'; nu = 5/2 }
  
  if(nugget_estim) nugget = 0.2
  if(!nugget_estim) nugget=1e-4
  
  colnames(X) = NULL
  data = data.frame(x=X, y=y)
  if(d>1) x.names = paste('x.', 1:d, sep = '')
  if(d==1) x.names = 'x'
  if (cov %in% c('Matern3_2', 'Matern5_2'))  param = c(variance=1, nugget=nugget, scale=1, nu=nu)
  else   param = c(variance=1, nugget=nugget, scale=1)

  locations_str = sprintf("~%s", paste(x.names, collapse = " + "))
  model = georob(y~1, data=data, locations = as.formula(locations_str),
         variogram.model=variogram.model, 
         param=param,
         fit.param = c(variance=TRUE, nugget=nugget_estim, scale=TRUE),
         tuning.psi = tuning.psi,
  )
  # variance = model$variogram.object[[1]]$param['variance'], 
  # lengthscale = model$variogram.object[[1]]$param['scale'],
  
  # outlier identification by using boxplot for residuals
  box = boxplot(model$residuals, plot = FALSE)
  out = which(model$residuals < box$stats[1] | model$residuals > box$stats[5])
  names(out) = NULL
  return(list(model=model, X=X, y=y, out = out))
}

# prediction function
GPRV.predict = function(object, Xnew){
  # object: result of GPRV
  # Xnew: new inputs
  # return: prediction at Xnew
  y=object$y
  colnames(Xnew) = NULL
  y_pred = predict(object$model, newdata = data.frame(x=Xnew))
  return(y_pred$pred)
}

######################################################################
# GPH in Gramacy and Ludkovski (2016)
# hetGP in Binois et al. (2018) 
#####################################################################
# main function
GPH = function(X, y, cov='Matern3_2', iso=FALSE){
  # X: design matrix
  # y: response
  # cov: covariance funtion:
  #      - gauss: Gaussian covariance
  #      - Matern3_2: Matern class with smooth parameter nu=3/2
  #      - Matern5_2: Matern class with smooth parameter nu=5/2
  # iso: a logical scalar: 
  #      TRUE for isotropic covariance
  #      FALSE for anisotropic covariance.
  # return list of model and out
  
  if(cov == 'gauss') covtype = 'Gaussian'
  if(cov == 'Matern3_2') covtype = 'Matern3_2'
  if(cov == 'Matern5_2') covtype = 'Matern5_2'
  
  X = as.matrix(X)
  y = as.vector(y)
  if(iso) model = mleHetGP(X, y, covtype = covtype, init = list(theta=1), settings = list(linkThetas='none'))
  if(!iso) model = mleHetGP(X, y, covtype = covtype, init = list(theta=rep(1, dim(X)[2])), settings = list(linkThetas='none'))
  
  # outlier identification by using boxplot for residuals
  y_pred = predict(object = model, x = X)$mean
  rds = y - y_pred
  box = boxplot(rds, plot = FALSE)
  out = which(rds > box$stats[5] | rds < box$stats[1])
  names(out) = NULL
  return(list(model=model, out = out))
}

# prediction function
GPH.predict = function(object, Xnew){
  # object: an object of class hetGP; e.g., as returned by mleHetGP
  # Xnew: a new design
  # return: prediction at Xnew
  model = object$model
  Xnew = as.matrix(Xnew)
  y_pred = predict(object = model, x=Xnew)$mean
  return(y_pred)
}

###########################################
# GPST in Vanhatalo et al. (2009)
###########################################
# main function
GPST = function(X, y, cov='Matern5_2', iso=TRUE, nu=4, scale = 1, approx='Laplace'){
  # X: design matrix
  # y: response
  # cov: covariance funtion:
  #      - exp: exponential covariance
  #      - sq_exp: square exponential covariance
  #      - gauss: Gaussian covariance
  #      - Matern3_2: Matern class with smooth parameter nu=3/2
  #      - Matern5_2: Matern class with smooth parameter nu=5/2
  # iso: a logical scalar: 
  #      TRUE for isotropic covariance
  #      FALSE for anisotropic covariance.
  #
  # nu, scale: the parameters for student t distribution
  # approx: approximation method:
  #         - Laplace: Laplace approximation
  #         - EP: Expectation Propagation approximation
  # return: list of model and out
  
  X = as.matrix(X)
  y = as.matrix(y)
  n = nrow(X)
  d = ncol(X)
  
  GPy = import('GPy')
  
  if(iso) {lengthscale=1; ARD = FALSE}
  else { lengthscale=as.vector(rep(1, d)); ARD = TRUE }
  
  if(cov == 'exp') kernel = GPy$kern$Exponential(input_dim = d, lengthscale=1)
  if(cov == 'sq_exp') kernel = GPy$kern$RBF(input_dim = d, lengthscale=lengthscale, ARD = ARD)
  if(cov == 'gauss') kernel = GPy$kern$RBF(input_dim = d, lengthscale=1, ARD = FALSE)
  if(cov == 'Matern3_2') kernel = GPy$kern$Matern32(input_dim = d, lengthscale=lengthscale, ARD = ARD)
  if(cov == 'Matern5_2') kernel = GPy$kern$Matern52(input_dim = d, lengthscale=lengthscale, ARD = ARD)
  
  t_distribution = GPy$likelihoods$StudentT(deg_free=nu, sigma2=scale)
  if(approx == 'Laplace') approx_inf = GPy$inference$latent_function_inference$Laplace()
  if(approx == 'EP') approx_inf = GPy$inference$latent_function_inference$EP()
  model = GPy$core$GP(X=X, Y=y, kernel=kernel, likelihood=t_distribution, inference_method=approx_inf)
  model$optimize()
  
  # variance = model$kern$parameters[0]
  # lengthscale = model$kern$parameters[1]
  
  # outlier identification by using boxplot for residuals
  y_pred = model$predict(X)[[1]][, 1]
  residuals = y - y_pred
  box = boxplot(residuals, plot = FALSE)
  out = which(residuals < box$stats[1] | residuals > box$stats[5])
  names(out) = NULL
  return(list(model=model, out = out))
}

# prediction function
GPST.predict = function(object, Xnew){
  # object: result of GPST
  # Xnew: new input
  # return: prediction at Xnew
  model = object$model
  Xnew = as.matrix(Xnew)
  y_pred = model$predict(Xnew)[[1]][, 1]
  return(y_pred)
}


###############################################################################
# MCMC-based GPST and SODK for comparison of runtime analysis in Supplementary.
# GPST-MCMC and SODK-MCMC are carried out by using the Python package 'gpflow'.
###############################################################################

###########################################################################
# GPST-MCMC: GP with student t likelihood - MCMC sampling for estimation
##########################################################################
# main function
GPSTMCMC = function(X, y, cov='Matern5_2', iso=TRUE, nu = 4, scale=.1, num_sample = 5000, num_burnin = 2000){
  # X: design matrix
  # y: response
  # cov: covariance funtion:
  #      - exp: exponential covariance
  #      - sq_exp: square exponential covariance
  #      - gauss: Gaussian covariance
  #      - Matern3_2: Matern class with smooth parameter nu=3/2
  #      - Matern5_2: Matern class with smooth parameter nu=5/2
  # iso: a logical scalar: 
  #      TRUE for isotropic covariance
  #      FALSE for anisotropic covariance.
  #
  # nu, scale: the parameters for student t distribution
  # num_sample: the number of samples at each sampling iteration
  # num_burnin: the number of burn-in samples at each sampling iteration
  # return: list of model and out
  
  X = as.matrix(X)
  y = as.vector(y)
  colnames(X) = NULL
  n = nrow(X)
  d = ncol(X)
  
  gpflow = import('gpflow')
  
  if(iso) lengthscales = 1
  if(!iso) lengthscales = as.vector(rep(1, d))
  if(cov == 'exp') kernel = gpflow$kernels$Exponential(lengthscales=lengthscales)
  if(cov == 'sq_exp') kernel = gpflow$kernels$SquaredExponential(lengthscales=lengthscales)
  if(cov == 'gauss') kernel = gpflow$kernels$RBF(lengthscales=lengthscales)
  if(cov == 'Matern3_2') kernel = gpflow$kernels$Matern32(lengthscales=lengthscales)
  if(cov == 'Matern5_2') kernel = gpflow$kernels$Matern52(lengthscales=lengthscales)
  
  model = GPSTMCMC_py(X=X, y=y, kernel = kernel, nu=nu, scale=scale, num_sample = num_sample, num_burnin = num_burnin)
  parameters = lapply(gpflow$utilities$parameter_dict(model), function(x){x$numpy()})
  scale = parameters$.likelihood.variance
  box = robustbase::adjbox(scale, plot = FALSE)
  out = which(scale > box$fence[2])
  return(list(model=model, out = out))
}

# prediction function
GPSTMCMC.predict = function(object, Xnew){
  # object: result of GPSTMCMC
  # Xnew: new inputs
  # return: prediction at Xnew
  model = object$model
  Xnew = as.matrix(Xnew)
  y_pred = model$predict_f_samples(Xnew)$numpy()[, 1]
  return(y_pred)
}

###########################################################################
# SODK-MCMC: SODK - MCMC sampling for estimation
##########################################################################
# main function
SODKMCMC = function(X, y, cov = 'Matern5_2', iso=TRUE, omega = 30, phi = 1, jitter=1e-4, num_sample=5000, num_burnin=2000){
  
  if(!(reticulate::py_available())){
    stop("Error message: the python module should be available. \n")
  }
  if(!(reticulate::py_module_available('gpflow'))){
    stop("Error message: the gpflow module should be available. \n")
  }
  
  X = as.matrix(X)
  y = as.vector(y)
  colnames(X) = NULL
  n = nrow(X)
  d = ncol(X)
  
  gpflow = import('gpflow')
  
  if(iso) lengthscales = 1
  if(!iso) lengthscales = as.vector(rep(1, d))
  
  if(cov == 'powexp') {
    kernel = gpflow$kernels$SquaredExponential(lengthscales=lengthscales)
  }
  if(cov == 'gauss') kernel = gpflow$kernels$RBF(lengthscales=lengthscales)
  if(cov == 'exp') kernel = gpflow$kernels$Exponential(lengthscales=lengthscales)
  if(cov == 'Matern3_2') kernel = gpflow$kernels$Matern32(lengthscales=lengthscales)
  if(cov == 'Matern5_2') kernel = gpflow$kernels$Matern52(lengthscales=lengthscales)
  
  model = ODKMCMC_py(X=X, y=as.matrix(y), kernel = kernel, omega=omega, jitter=jitter, num_sample = num_sample, num_burnin = num_burnin)
  parameters = lapply(gpflow$utilities$parameter_dict(model), function(x){x$numpy()})
  scale = parameters$.likelihood.variance
  box = robustbase::adjbox(scale, plot = FALSE)
  out = which(scale > box$fence[2])
  return(list(model=model, out = out))
}

# prediction function
SODKMCMC.predict = function(object, Xnew){
  # object: result of ODKMCMC
  # Xnew: new inputs
  # return: prediction at Xnew
  
  model = object$model
  Xnew = as.matrix(Xnew)
  y_pred = model$predict_f_samples(Xnew)$numpy()[, 1]
  return(y_pred)
}
