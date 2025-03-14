###########################################
############# simulation 2 ################
###########################################

# user must specify the correct root path load the R scripts
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))
source(paste(root_path,'/benchmark/experiment/sim2-helping-functions.R', sep = ''))

# load the required packages
library(lhs)
library(doSNOW)  # progress bar
library(foreach) # parallel computing
library(ggplot2) # for visualization


########### start to simulation #####################
# Parameter setting
swapping = c('Random','EH')
beta = c(0.05, 0.10, 0.15) # contamination level
rho = 0.7
nu = 2.5
sample_size = c(300, 400, 500, 600) #

parameters = expand.grid(1, beta, swapping, sample_size)
colnames(parameters) = c("seeds", "beta", "swapping", 'sample_size')
npar = dim(parameters)[1]
#
MCMCruntime.analysis = function(i, par, cov='Matern5_2', k=5){
  # i: i-th simulation
  # par: parameters setting:
  #      - par[i, 1]: seed
  #      - par[i, 2]: beta
  #      - par[i, 3]: swapping: 'EH'- switch observations that are most different according to first (global) eigenvector score
  #                             'Random' - switches observations completely random
  #      - par[i, 4]: sample size
  # cov: matern class for all considered models
  # nu: smoothness parameter for bivariate Matern model
  # k: if one point is swapped, then remove its k nearest neighbour points from the following swapping process
  
  
  seed = par[i, 1]
  beta = par[i, 2]
  swapping = par[i, 3]
  n = par[i, 4]
  nu=2.5
  rho=0.7
  
  cat('Runtime analysis: The ', i, '/', nrow(par), '-th replication: \n', sep = '')
  cat('design=LHS,',' beta=', beta, ', swapping=', as.character(swapping), ', n=', n, '\n', sep = '')
      
  location = randomLHS(n, 2)*10
  matern_param = list(mean_1 = 0, mean_2 = 0, smooth_1 = nu, smooth_2 = nu, smooth_12 = nu,
                      scale_1 = 1, scale_2 = 1, scale_12 = 1,
                      sill_1 = 1, sill_2 = 1, nugget_1 = 0, nugget_2 = 0, pcol = rho)
  
  
  
  data = contaminated_data_randomfields(seed = seed, coords = location, param = matern_param, beta=beta, type = swapping, k = k)
  X = as.matrix(data[, c('x', 'y')])
  y = as.vector(data[, 'variable1'])

  ################################### start to fit ###################################

  # GPH
  temp = try({
    b1 = as.numeric(Sys.time())
    GPH.m = GPH(X=X, y = y,cov = cov)
    b2 = as.numeric(Sys.time())
    GPH.time = eval_runtimes(b1, b2, "secs")
    GPH.FNR = classification(data = data, outliers = GPH.m$out)[1]
    GPH.FPR = classification(data = data, outliers = GPH.m$out)[2]
    GPH.F1 = F_measure(fnr = GPH.FNR, fpr = GPH.FPR, beta=beta, n = n)
    cat("GPH model has been executed and took ", GPH.time, "secs.", "\n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPH.FNR = GPH.FPR = GPH.F1 = GPH.time = NA
  } 
  # GPRV
  temp = try({
    c1 = as.numeric(Sys.time())
    GPRV.m = GPRV(X=X, y=y, cov = cov)
    c2 = as.numeric(Sys.time())
    GPRV.time = eval_runtimes(c1, c2, "secs")
    GPRV.FNR = classification(data = data, outliers = GPRV.m$out)[1]
    GPRV.FPR = classification(data = data, outliers = GPRV.m$out)[2]
    GPRV.F1 = F_measure(fnr = GPRV.FNR, fpr = GPRV.FPR, beta=beta, n = n)
    cat("GPRV model has been executed and took ", GPRV.time, "secs.", "\n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPRV.FNR = GPRV.FPR = GPRV.F1 = GPRV.time = NA
  }
  
  # GPST
  temp = try({
    d1 = as.numeric(Sys.time())
    GPST.m = GPST(X=X, y=y, cov = cov, nu=4,approx = 'Laplace')
    d2 = as.numeric(Sys.time())
    GPST.time = eval_runtimes(d1, d2, "secs")
    GPST.FNR = classification(data = data, outliers = GPST.m$out)[1]
    GPST.FPR = classification(data = data, outliers = GPST.m$out)[2]
    GPST.F1 = F_measure(fnr = GPST.FNR, fpr = GPST.FPR, beta=beta, n = n)
    cat("GPST model has been executed and took ", GPST.time, "secs.", "\n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPST.FNR =GPST.FPR = GPST.F1 = GPST.time = NA
  } 
  
  # SODK
  temp = try({
    f1 = as.numeric(Sys.time())
    SODK.m = SODK(X=X, y=y, cov = cov, phi = .1, phi_estim = TRUE, omega = 30, u_threshold = 1e-2)
    f2 = as.numeric(Sys.time())
    SODK.time = eval_runtimes(f1, f2, "secs")
    SODK.FNR = classification(data = data, outliers = SODK.m$out)[1]
    SODK.FPR = classification(data = data, outliers = SODK.m$out)[2]
    SODK.F1 = F_measure(fnr = SODK.FNR, fpr = SODK.FPR, beta=beta, n = n)
    cat("SODK model has been executed and took ", SODK.time, "secs.", "\n")
  }, silent = FALSE)
  
  if('try-error' %in% class(temp)) {
    SODK.FNR = SODK.FPR = SODK.F1 = SODK.time = NA
  }
  
  ########### MCMC methods #################
  # GPST-MCMC 
  temp = try({
    d1 = as.numeric(Sys.time())
    GPSTMCMC.m = GPSTMCMC(X=X, y=y, cov = cov, nu=4, num_sample = 5000, num_burnin = 2000)
    d2 = as.numeric(Sys.time())
    GPSTMCMC.time = eval_runtimes(d1, d2, "secs")
    GPSTMCMC.FNR = classification(data = data, outliers = GPSTMCMC.m$out)[1]
    GPSTMCMC.FPR = classification(data = data, outliers = GPSTMCMC.m$out)[2]
    GPSTMCMC.F1 = F_measure(fnr = GPSTMCMC.FNR, fpr = GPSTMCMC.FPR, beta=beta, n = n)
    cat("GPST-MCMC model has been executed and took ", GPSTMCMC.time, "secs.", "\n")
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPSTMCMC.FNR = GPSTMCMC.FPR = GPSTMCMC.F1 = GPSTMCMC.time = NA
  }
  
  # SODK-MCMC
  temp = try({
    f1 = as.numeric(Sys.time())
    SODKMCMC.m = SODKMCMC(X=X, y=y, cov = cov, omega = 30, num_sample = 5000, num_burnin = 2000)
    f2 = as.numeric(Sys.time())
    SODKMCMC.time = eval_runtimes(f1, f2, "secs")
    SODKMCMC.FNR = classification(data = data, outliers = SODKMCMC.m$out)[1]
    SODKMCMC.FPR = classification(data = data, outliers = SODKMCMC.m$out)[2]
    SODKMCMC.F1 = F_measure(fnr = SODKMCMC.FNR, fpr = SODKMCMC.FPR, beta=beta, n = n)
    cat("ODKMCMC model has been executed and took ", ODKMCMC.time, "secs.", "\n")
  }, silent = FALSE)
  
  if('try-error' %in% class(temp)) {
    SODKMCMC.FNR = SODKMCMC.FPR = SODKMCMC.F1 = SODKMCMC.time = NA
  }
  
  gc() # release the cache
  
  result = c(GPH.FNR, GPH.FPR, GPH.F1, GPH.time, 
             GPRV.FNR, GPRV.FPR, GPRV.F1, GPRV.time,
             GPST.FNR, GPST.FPR, GPST.F1, GPST.time,
             SODK.FNR, SODK.FPR, SODK.F1, SODK.time,
             GPSTMCMC.FNR, GPSTMCMC.FPR, GPSTMCMC.F1, GPSTMCMC.time,
             SODKMCMC.FNR, SODKMCMC.FPR, SODKMCMC.F1, SODKMCMC.time)
  return(result)
}


MCMCruntime.result = foreach(i=1:nrow(parameters), .combine = 'rbind', 
                         .packages = c('reticulate','lhs','GeoModels','georob', 'hetGP','rrcov')
) %do% MCMCruntime.analysis(
  i=i,
  par=parameters,
  k = 5
)

stopCluster(cl)


MCMCruntime.result = data.frame(MCMCruntime.result, parameters$beta, parameters$swapping, parameters$sample_size)
colnames(MCMCruntime.result) = c('GPH.FNR', 'GPH.FPR', 'GPH.F1', 'GPH.time',
                                'GPRV.FNR', 'GPRV.FPR', 'GPRV.F1', 'GPRV.time', 
                                'GPST.FNR', 'GPST.FPR', 'GPST.F1', 'GPST.time', 
                                'SODK.FNR', 'SODK.FPR', 'SODK.F1', 'SODK.time',
                                'GPSTMCMC.FNP', 'GPSTMCMC.FPR', 'GPSTMCMC.F1','GPSTMCMC.time', 
                                'SODKMCMC.FNR', 'SODKMCMC.FPR', 'SODKMCMC.F1', 'SODKMCMC.time',
                                'beta', 'swapping', 'sample_size')

# save(MCMCruntime_result, file = paste(root_path, '/result/supp/MCMCruntime-',Sys.Date(), '.Rdata', sep = ''))


############### plots of FNR, RPR, F1 score #######################
MCMCruntime.plot <- data.frame(
  Method = rep(c('GPH',"GPRV", 'GPST', 'SODK', 'GPST-MCMC', 'SODK-MCMC'), each = 24),
  beta = rep(paste('beta = ', MCMCruntime.result$beta,sep = ''), 6),
  swapping = rep(MCMCruntime_result$swapping, 6),
  n = rep(MCMCruntime_result$sample_size, 6),
  Runtime = c(MCMCruntime.result$GPH.time, MCMCruntime.result$GPRV.time, MCMCruntime.result$GPST.time, 
              MCMCruntime.result$SODK.time, MCMCruntime.result$GPSTMCMC.time, MCMCruntime.result$SODKMCMC.time)
)


Runtime = ggplot(MCMCruntime.plot, aes(x = n, y = Runtime, color = Method)) + 
  geom_line(aes(color = Method), linewidth=1) + 
  geom_point(aes(shape = Method),size=3) + 
  labs(x=expression(paste('Sample size (x', , 10^2,')', sep = '')) , y = "Runtime[s]") +
  facet_grid(~swapping~beta) + 
  scale_x_continuous(breaks = c(300, 400, 500, 600),
                     labels = c('3', '4', '5', '6')) +
  scale_color_manual(values = c(2, 3, 4, 6, 5, 7)) +
  scale_shape_manual(values = c(17, 15, 3, 4, 7, 16)) + 
  theme_bw() + 
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))

