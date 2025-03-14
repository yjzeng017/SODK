# Supplementary: Thresholding analysis


# user must specify the correct path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
# load R scripts
source(paste(root_path, '/benchmark/experiment/sim2-helping-functions.R', sep = ''))
source(paste(root_path, '/SODK.R', sep = ''))
#
library(doSNOW)  # progress bar
library(foreach) # parallel computing


# Parameter setting
swapping = c('EH', 'Random')
design = c('Grid')
beta = seq(5, 15, 5)/100 # contamination level
parameters = expand.grid(1:100, beta, swapping, design)
colnames(parameters) = c("seeds", "beta", "swapping", 'design')

########### start to simulation #####################
convergence.analysis = function(i, par, cov='Matern5_2', k=5, max_iter=30){
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
  # max_iter: the maximal number of iterations
  
  seed = par[i, 1]
  beta = par[i, 2]
  swapping = par[i, 3]
  design = par[i, 4]
  
  nu=2.5
  rho=0.7
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
  
  # the outlier detection based SODK for various thresholds
  temp = try({
    SODK.m = SODK(X=X, y=y, cov = cov, phi = 1, phi_estim = TRUE, omega = 30, u_threshold = 0, change_between_iter = 0, max_iter = max_iter)
    theta_gap = SODK.m$theta_gap
    u_gap = SODK.m$u_gap
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    theta_gap = rep(NA, max_iter)
    u_gap = rep(NA, max_iter)
  }

  return(c(u_gap, theta_gap))
}


# message(npar)

##########################################################################################
# set up progress bar
npar = dim(parameters)[1]
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(15)
registerDoSNOW(cl)

convergence.results = foreach(i=1:nrow(parameters), .combine = 'rbind', .options.snow = opts,
                              .packages = c('reticulate','lhs', 'GeoModels', 'rrcov')
) %dopar% convergence.analysis(
  i=i,
  par=parameters,
  k = 5,
)
stopCluster(cl)

convergence.results = data.frame(convergence.results, parameters$beta, parameters$swapping, parameters$design)
colnames(convergence.results) = c(paste('u_gap', 1:30, sep = '_'), paste('theta_gap', 1:30, sep = '_'), 'beta', 'swapping', 'design')
# save(convergence.results, file=paste(root_path, '/result/supp/convergence-', Sys.Date(), '.Rdata', sep = ''))


############### plots of FNR, RPR, F1 score #######################
library(ggplot2)
library(dplyr)

convergence.plot <- data.frame(
  swapping = rep(convergence.results$swapping, 30),
  beta = rep(paste('beta=', convergence.results$beta, sep = ''), 30),
  u_gap = unlist(convergence.results[,1:30]),
  theta_gap = unlist(convergence.results[, 31:60]),
  iter = rep(1:30, each = 600)
)

u.plot = ggplot(convergence.plot, aes(x = iter, y = u_gap, group = iter)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) + 
  labs(x=paste('The number of iteration'), y=expression(paste('||', u[n]^t, '-', u[n]^'t-1', '||')))


theta.plot = ggplot(convergence.plot, aes(x = iter, y = theta_gap, group = iter)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x=paste('The number of iteration'), y=expression(paste('||', theta^t, '-', theta^'t-1', '||')))


gridExtra::grid.arrange(u_gap, theta_gap, ncol = 1)
