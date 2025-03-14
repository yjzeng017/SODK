setwd("D:/SODK")
source('./benchmark/experiment/sim2/sim2-helping-functions.R')
library(doSNOW)  # progress bar
library(foreach) # parallel computing

# Parameter setting
swapping = c('EH', 'Random')
design = c('Grid', 'LHS')
beta = seq(1, 20, 1)/100 # contamination level
rho = 0.7
nu = 2.5
parameters = expand.grid(1:100, beta, swapping, design)
colnames(parameters) = c("seeds", "beta", "swapping", 'design')

# 
pb <- txtProgressBar(max = dim(parameters)[1], style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl = makeCluster(15)
registerDoSNOW(cl)

foreach(i=1:dim(parameters)[1], .options.snow = opts, .final = function(x) '\n Done. \n', .packages = c('lhs','GeoModels','rrcov')) %dopar% {
  seed = parameters[i, 1]
  beta_i = parameters[i, 2]
  swapping_i = parameters[i, 3]
  design_i = parameters[i, 4]
  
  if(design_i == 'Grid') location = expand.grid(seq(0, 10, 0.5), seq(0, 10, 0.5))
  if(design_i == 'LHS') location = lhs::randomLHS(n=441, 2)*10
  
  # specify the parameters for bivariate Matern model
  matern_param = list(mean_1 = 0, mean_2 = 0, smooth_1 = nu, smooth_2 = nu, smooth_12 = nu,
                      scale_1 = 1, scale_2 = 1, scale_12 = 1,
                      sill_1 = 1, sill_2 = 1, nugget_1 = 0, nugget_2 = 0, pcol = rho)
  
  # generate the simulated data
  data = contaminated_data_randomfields(seed = seed, coords = location, param = matern_param, beta=beta_i, type = swapping_i, k = 5)
  data['beta'] = rep(beta_i, dim(data)[1])
  data['swapping'] = rep(swapping_i, dim(data)[1])
  data['design'] = rep(design_i, dim(data)[1])
  
  save(data, file = paste('./benchmark/experiment/sim2/sim2-data-OD/data_', i,'.Rdata', sep = ''))
  return(NULL)
}
stopCluster(cl)
