###########################################
# simulation 2- outlier detection 
###########################################

# user must specify the correct root path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"

# load the R scripts for simulation
source(paste(root_path,'/benchmark/experiment/sim2-helping-functions.R', sep = ''))
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))

# load the packages for parallel computation
library(doSNOW)  # progress bar
library(foreach) # parallel computing

########### start to simulation #####################
# Parameter setting
swapping = c('EH', 'Random')
design = c('Grid', 'LHS')
beta = seq(5, 20, 5)/100 # contamination level
rho = 0.7
nu = 2.5
parameters = expand.grid(1:100, beta, swapping, design)
colnames(parameters) = c("seeds", "beta", "swapping", 'design')


########### start to simulation #####################
prediction = function(i, par, cov='Matern5_2', nu=2.5, rho=0.7, k=5, test.size=100){
  # i: i-th simulation
  # par: parameters setting:
  #      - par[i, 1]: seed
  #      - par[i, 2]: beta
  #      - par[i, 3]: swapping: 'EH'- switch observations that are most different according to first (global) eigenvector score
  #                             'Random' - switches observations completely random
  #      - par[i, 4]: design: 'Grid' for grids design, 'LHS' for LHS design
  # cov: matern class for all considered models
  # nu: smoothness parameter for bivariate Matern model
  # k: if one point is swapped, then remove its k nearest neighbour points from the following swapping process
  # test.size: size of testing sample
  
  
  seed = par[i, 1]
  beta = par[i, 2]
  swapping = par[i, 3]
  design = par[i, 4]
  cat('Simulation 2: The', i, '/', nrow(par), '-th replication with', 'nu=', nu, 'rho=', rho, 'beta=', beta, 'swapping=', as.character(swapping), '\n')
  
  if(design == 'Grid') location = expand.grid(seq(0, 10, 0.5), seq(0, 10, 0.5))
  if(design == 'LHS') location = lhs::randomLHS(n=441, 2)*10
  matern_param = list(mean_1 = 0, mean_2 = 0, smooth_1 = nu, smooth_2 = nu, smooth_12 = nu,
                      scale_1 = 1, scale_2 = 1, scale_12 = 1,
                      sill_1 = 1, sill_2 = 1, nugget_1 = 0, nugget_2 = 0, pcol = rho)
  
  
  
  data = contaminated_data_randomfields(seed = seed, coords = location, param = matern_param, beta=beta, type = swapping, k = k)
  X = as.matrix(data[, c('x', 'y')])
  y = as.vector(data[, 'variable1'])
  n = nrow(X)
  
  set.seed(seed)
  out = which(data$out == 1)
  test.loc = sample((1:n)[-out], size = test.size)
  X_train = X[-test.loc, ]
  y_train = y[-test.loc]
  X_test = X[test.loc, ]
  y_test = y[test.loc]
  ###################################
  # performance of prediction
  ###################################
  # CK
  temp = try({
    a1 = as.numeric(Sys.time())
    CK.m = SODK(X=X_train, y=y_train, cov = cov, phi = 0)
    a2 = as.numeric(Sys.time())
  }, silent = FALSE)
  
  if('try-error' %in% class(temp)) {
    CK.RMSE = NA
  } else{
    CK.time = eval_runtimes(a1, a2, "secs")
    CK.pred = SODK.predict(object = CK.m, Xnew = as.matrix(X_test))
    CK.RMSE = sqrt(mean((y_test - CK.pred)**2))
    cat("CK model has been successfully executed in ", CK.time, "secs.", "\n")
  }
  
  # GPH
  temp = try({
    b1 = as.numeric(Sys.time())
    GPH.m = GPH(X=X_train, y = y_train, cov = cov)
    b2 = as.numeric(Sys.time())
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPH.RMSE = NA
  } else{
    GPH.time = eval_runtimes(b1, b2, "secs")
    GPH.pred = GPH.predict(GPH.m, Xnew = X_test)
    GPH.RMSE = sqrt(mean((y_test - GPH.pred)**2))
    cat("GPH model has been successfully executed in ", GPH.time, "secs.", "\n")
  }
  
  # GPRV
  temp = try({
    c1 = as.numeric(Sys.time())
    GPRV.m = GPRV(X=X_train, y=y_train, cov = cov)
    c2 = as.numeric(Sys.time())
    GPRV.pred = GPRV.predict(object = GPRV.m, Xnew = X_test)
    GPRV.time = eval_runtimes(c1, c2, "secs")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPRV.RMSE = NA
  } else {
    GPRV.RMSE = sqrt(mean((y_test - GPRV.pred)**2))
    cat("GPRV model has been successfully executed in ", GPRV.time, "secs.", "\n")
  }
  
  # GPST
  temp = try({
    d1 = as.numeric(Sys.time())
    GPST.m = GPST(X=X_train, y=y_train, cov = cov)
    d2 = as.numeric(Sys.time())
    GPST.time = eval_runtimes(d1, d2, "secs")
    GPST.pred = GPST.predict(object = GPST.m, Xnew = X_test)
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPST.RMSE = NA
  } else{
    GPST.RMSE = sqrt(mean((y_test - GPST.pred)**2))
    cat("GPST model has been successfully executed in ", GPST.time, "secs.", "\n")
  }
  
  # SODK
  temp = try({
    f1 = as.numeric(Sys.time())
    SODK.m = SODK(X=X_train, y=y_train, cov = cov, phi = 1, phi_estim = TRUE, omega = 30, u_threshold = 1e-2)
    f2 = as.numeric(Sys.time())
  }, silent = FALSE)
  
  if('try-error' %in% class(temp)) {
    SODK.RMSE = NA
  } else{
    SODK.time = eval_runtimes(f1, f2, "secs")
    SODK.pred = SODK.predict(object = SODK.m, Xnew = as.matrix(X_test))
    SODK.RMSE = sqrt(mean((y_test - SODK.pred)**2))
    cat("SODK model has been successfully executed in ", SODK.time, "secs.", "\n")
  }
  
  gc() # release the cache
  
  result = c(CK.RMSE, GPH.RMSE, GPRV.RMSE, GPST.RMSE, SODK.RMSE)
  return(result)
}


# set up progress bar
npar = dim(parameters)[1]
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


# adapt number of kernels used
cl <- makeCluster(15)
registerDoSNOW(cl)

prediction.result = foreach(i=1:nrow(parameters), .combine = 'rbind', .options.snow = opts,
                            .packages = c('reticulate','lhs','GeoModels','georob', 'hetGP','rrcov','robustbase')
) %dopar% prediction(
  i=i,
  par=parameters,
  k = 5,
  test.size=100
)

stopCluster(cl)


prediction.result = data.frame(prediction.result, parameters$beta, parameters$swapping, parameters$design)
colnames(prediction.result) = c('CK.RMSE', 'GPH.RMSE', 'GPRV.RMSE', 'GPST.RMSE', 'SODK.RMSE', 'beta', 'swapping', 'design')


save(prediction.result, file = paste(root_path, './result/sim2-RMSE-', Sys.Date(), '.Rdata', sep = ''))

# 


# # show the average RMSE
library(dplyr)
RMSE.mean = na.omit(prediction.result) %>% group_by(beta, swapping, design) %>%
  summarise(across(everything(), mean)) %>% group_by(swapping, beta, design) %>%
  summarise(across(everything(), function(x){round(x, 2)}))


RMSE.sd = na.omit(prediction.result) %>% group_by(beta, swapping, design) %>%
  summarise(across(everything(), sd)) %>% group_by(swapping, beta, design) %>%
  summarise(across(everything(), function(x){round(x, 2)}))
