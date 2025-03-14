# Supplementary: additional simulations: spatial data

# user must specify the correct path

# user must specify the correct root path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
#
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))

library('reticulate')
library(MASS)
library(doSNOW)  # progress bar
library(foreach) # parallel computing


# measures for outlier detection
M_S_JD = function(out, true_out, n){
  normal_point = (1:n)[-true_out]
  M = length(setdiff(true_out, out))/length(true_out)
  S = length(intersect(out, normal_point))/length(normal_point)
  
  if(identical(intersect(out, true_out), true_out )) JD = 1
  if(!identical(intersect(out, true_out), true_out )) JD = 0
  
  return(c(M,S, JD))
}


# Parameter setting
l = c(5, 10, 15) # domain [0,l]^2, n = 4*l^2
out.types = c('extreme', 'skewed', 'uniform')
beta = c(0.05, 0.10)
parameters = expand.grid(1:100, out.types, beta, l)
colnames(parameters) = c("seeds", 'out.type', 'beta', 'l')


spatial_simulation = function(i, par, cov='gauss'){
  # i = 3
  # par = parameters
  
  seed = par[i, 1]
  out.type = par[i, 2]
  beta = par[i, 3]
  l = par[i, 4]

  cat('zmore simulations: spatial example. The', i, '/', nrow(par), '-th replication with', 'n =', 4*l^2, 'contamination = ', beta, 'outlier type =', as.character(out.type), '\n')
  set.seed(seed)
  
  n = as.integer(4*l^2)
  X = cbind(runif(2*n , 0, l), runif(2*n , 0, l))
  
  kernel =  kernlab::rbfdot(sigma = 0.5)
  K = kernlab::kernelMatrix(kernel, X)
  K = K@.Data
  y = MASS::mvrnorm(n=1, mu = rep(0, 2*n), Sigma = K)
  
  test.loc = sample(1:dim(X)[1], size = n)
  X_train = X[-test.loc, ]
  y_train = y[-test.loc]
  X_test = X[test.loc, ]
  y_test = y[test.loc]
  out = sort(sample(1:n, size = as.integer(n*beta)))
  if(out.type == 'extreme') bias = rnorm(length(out), 0, 5)
  if(out.type == 'skewed') bias = rnorm(length(out), 3, 1)
  if(out.type == 'uniform') bias = runif(length(out), 2, 3)
  y_train[out] = y_train[out] + bias
  
  ###################################
  # performance of prediction
  ###################################
  # CK with outliers
  temp = try({
    CK.m = SODK(X=X_train, y=y_train, cov = cov, phi = 0)
    CK.pred = SODK.predict(object = CK.m, Xnew = as.matrix(X_test))
    CK.RMSE = sqrt(mean((y_test - CK.pred)**2))
    cat("CK model has been successfully executed. \n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    CK.RMSE = NA
  }
  
  # GPH (noises are smoothed and modeled by another GP)
  temp = try({
    GPH.m = GPH(X=X_train, y = y_train, cov = cov)
    GPH.MSJD = M_S_JD(GPH.m$out, out, n)
    GPH.pred = GPH.predict(GPH.m, Xnew = X_test)
    GPH.RMSE = sqrt(mean((y_test - GPH.pred)**2))
    cat("GPH model has been successfully executed. \n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPH.RMSE = NA
    GPH.MSJD = rep(NA, 3)
  }
  
  # GPRV (robust variogram estimation)
  temp = try({
    GPRV.m = GPRV(X=X_train, y=y_train, cov = cov)
    GPRV.pred = GPRV.predict(object = GPRV.m, Xnew = X_test)
    GPRV.RMSE = sqrt(mean((y_test - GPRV.pred)**2))
    GPRV.MSJD = M_S_JD(GPRV.m$out, out, n)
    cat("GPRV model has been successfully executed. \n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPRV.RMSE = NA
    GPRV.MSJD = rep(NA, 3)
  }
  
  # GPST (Laplace approximation)
  temp = try({
    GPST.m = GPST(X=X_train, y=y_train, cov = cov)
    GPST.MSJD = M_S_JD(GPST.m$out, out, n)
    GPST.pred = GPST.predict(object = GPST.m, Xnew = X_test)
    GPST.RMSE = sqrt(mean((y_test - GPST.pred)**2))
    cat("GPST model has been successfully executed. \n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPST.RMSE = NA
    GPST.MSJD = rep(NA, 3)
  }
  
  # proposed SODK model (H-likelihood procedure)
  temp = try({
    SODK.m = SODK(X=X_train, y=y_train, cov = cov, phi = 1, phi_estim = TRUE, omega = 30, u_threshold = 1e-2)
    SODK.MSJD = M_S_JD(SODK.m$out, out, n)
    SODK.pred = SODK.predict(object = SODK.m, Xnew = as.matrix(X_test))
    SODK.RMSE = sqrt(mean((y_test - SODK.pred)**2))
    cat("SODK model has been successfully executed. \n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    SODK.RMSE = NA
    SODK.MSJD = rep(NA, 3)
  }
  
  result = c(CK.RMSE, GPH.RMSE, GPRV.RMSE, GPST.RMSE, SODK.RMSE,
             GPH.MSJD, GPRV.MSJD, GPST.MSJD, SODK.MSJD)
  
  return(result)
}


##########################################################################################
# set up progress bar
npar = dim(parameters)[1]
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(10)
registerDoSNOW(cl)


spatial.result = foreach(i=1:nrow(parameters), .combine = 'rbind', .options.snow = opts,
                        .packages = c('reticulate', 'MASS', 'kernlab', 'hetGP', 'georob')
) %dopar% spatial_simulation(
  i=i,
  par=parameters,
)

stopCluster(cl)


spatial.result = data.frame(spatial.result, parameters$out.type, parameters$beta, parameters$l)
colnames(spatial.result) = c('CK.RMSE', 'GPH.RMSE', 'GPRV.RMSE', 'GPST.RMSE', 'SODK.RMSE', 
                        'GPH.M', 'GPH.S', 'GPH.JD',
                        'GPRV.M', 'GPRV.S', 'GPRV.JD',
                        'GPST.M', 'GPST.S', 'GPST.JD',
                        'SODK.M', 'SODK.S', 'SODK.JD',
                        'out.type', 'beta', 'l'
                        )

save(spatial.result, file = paste(root_path, '/result/supp/supp-spatial-', Sys.Date(), '.Rdata', sep = ''))

# # show the average RMSE
library(dplyr)
RMSE.mean = na.omit(spatial.result) %>% group_by(l, beta, out.type) %>%
  summarise(across(everything(), mean)) %>% group_by(l, beta, out.type) %>%
  summarise(across(everything(), function(x){round(x, 2)}))

RMSE.sd = na.omit(spatial.result) %>% group_by(l, beta, out.type) %>%
  summarise(across(everything(), sd)) %>% group_by(l, beta, out.type) %>%
  summarise(across(everything(), function(x){round(x, 2)})) 
