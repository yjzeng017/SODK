# Supplementary: Thresholding analysis


# user must specify the correct path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
# load R scripts
source(paste(root_path, '/benchmark/experiment/sim2-helping-functions.R', sep = ''))
source(paste(root_path, '/SODK.R', sep = ''))
# 
library(doSNOW)  # progress bar
library(foreach) # parallel computing


########### start to simulation #####################
# Parameter setting
swapping = c('EH', 'Random')
design = c('Grid')
# beta = seq(2, 10, 2)/100 # contamination level
beta = seq(5, 15, 5)/100 # contamination level
rho = 0.7
nu = 2.5
u_threshold = c(1, 0.1, 0.01, 0.001, 0.0001)
parameters = expand.grid(1:100, beta, swapping, design, u_threshold)
colnames(parameters) = c("seeds", "beta", "swapping", 'design', 'u_threshold')

# 
threshold.analysis = function(i, par, cov='Matern5_2', k=5){
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

  seed = par[i, 1]
  beta = par[i, 2]
  swapping = par[i, 3]
  design = par[i, 4]
  u_threshold = par[i, 5]
  
  rho = 0.7
  nu = 2.5
  
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
  
  temp = try({
    f1 = as.numeric(Sys.time())
    SODK.m = SODK(X=X, y=y, cov = cov, phi = 1, phi_estim = TRUE, omega = 30, u_threshold = u_threshold)
    f2 = as.numeric(Sys.time())
  }, silent = FALSE)
    
  if('try-error' %in% class(temp)) {
    SODK.FNR = SODK.FPR = SODK.F1 = NA
  } else{
    SODK.time = eval_runtimes(f1, f2, "secs")
    out = SODK.m$out
    SODK.FNR = classification(data = data, outliers = out)[1]
    SODK.FPR = classification(data = data, outliers = out)[2]
    SODK.F1 = F_measure(fnr = SODK.FNR, fpr = SODK.FPR, beta=beta, n = n)
  }
  
  return(c(SODK.FNR, SODK.FPR, SODK.F1))
}


##########################################################################################
# set up progress bar
npar = dim(parameters)[1]
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(15)
registerDoSNOW(cl)

threshold.results = foreach(i=1:nrow(parameters), .combine = 'rbind', .options.snow = opts,
                  .packages = c('reticulate','lhs', 'GeoModels', 'rrcov')
) %dopar% threshold.analysis(
  i=i,
  par=parameters,
  k = 5,
)
stopCluster(cl)


threshold.results = data.frame(threshold.results, parameters$beta, parameters$swapping, parameters$u_threshold)
colnames(threshold.results) = c('FNR', 'FPR', 'F1', 'beta', 'swapping', 'u_threshold')
# save(threshold.results, file = paste(root_path, '/result/supp/threshold-', Sys.Date(),'.Rdata', sep = ''))


############### plots of FNR, RPR, F1 score #######################
library(ggplot2)
library(dplyr)

# threshold.results.mean = na.omit(threshold.results) %>% 
#   group_by(beta, swapping, u_threshold) %>%
#   summarise(across(everything(), mean))


threshold.plot <- data.frame(
  swapping = threshold.results$swapping,
  threshold = factor(rep(c('1', '1e-1', '1e-2', '1e-3', '1e-4'), each = 600)),
  beta = paste('beta=',threshold.results$beta, sep = ''),
  FNR = threshold.results$FNR,
  FPR = threshold.results$FPR,
  F1 = threshold.results$F1
)

FPR = ggplot(threshold.plot, aes(x = threshold, y = FPR, group = threshold)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) + 
  labs(x='')

FNR = ggplot(threshold.plot, aes(x = threshold, y = FNR, group = threshold)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) + 
  labs(x='')


F1 = ggplot(threshold.plot, aes(x = threshold, y = F1, group = threshold)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x='Threshold', y='F1 scores')


gridExtra::grid.arrange(FPR, FNR, F1, ncol = 1)
