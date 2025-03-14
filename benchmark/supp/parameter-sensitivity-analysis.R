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
beta = seq(5, 15, 5)/100 # contamination level
omega = c(5, 10, 15, 20, 25, 30, 35, 40)
parameters = expand.grid(1:3, beta, swapping, design, omega)
colnames(parameters) = c("seeds", "beta", "swapping", 'design', 'omega')

sensitivity.analysis = function(i, par, cov='Matern5_2', k=5, test.size=100){
  # i: i-th simulation
  # par: parameters setting:
  #      - par[i, 1]: seed
  #      - par[i, 2]: beta
  #      - par[i, 3]: swapping: 'EH'- switch observations that are most different according to first (global) eigenvector score
  #                             'Random' - switches observations completely random
  #      - par[i, 4]: design: 'Grid' for grids design, 'LHS' for LHS design
  #      - par[i, 5]: omega: regularization parameter
  
  # cov: matern class for all considered models
  # nu: smoothness parameter for bivariate Matern model
  # k: if one point is swapped, then remove its k nearest neighbour points from the following swapping process
  
  
  seed = par[i, 1]
  beta = par[i, 2]
  swapping = par[i, 3]
  design = par[i, 4]
  omega = par[i, 5]
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
    f1 = as.numeric(Sys.time())
    SODK.m = SODK(X=X, y=y, cov = cov, phi = 1, phi_estim = TRUE, omega = omega, u_threshold = 1e-2)
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
  
  # set.seed(seed)
  # out = which(data$out == 1)
  # test.loc = sample((1:n)[-out], size = test.size)
  # X_train = X[-test.loc, ]
  # y_train = y[-test.loc]
  # X_test = X[test.loc, ]
  # y_test = y[test.loc]
  # 
  # # CK with outliers
  # temp = try({
  #   SODK.m = SODK(X=X_train, y=y_train, cov = cov, phi = 1, phi_estim = TRUE, omega = omega, u_threshold = 1e-2)
  # }, silent = FALSE)
  # 
  # if('try-error' %in% class(temp)) {
  #   SODK.RMSE = NA
  # } else{
  #   SODK.pred = SODK.predict(object = SODK.m, Xnew = as.matrix(X_test))
  #   SODK.RMSE = sqrt(mean((y_test - SODK.pred)**2))
  # }
  
  return(c(SODK.FNR, SODK.FPR, SODK.F1))
}


##########################################################################################
# set up progress bar
npar = dim(parameters)[1]
pb <- txtProgressBar(max = npar, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(18)
registerDoSNOW(cl)

sensitivity.results = foreach(i=1:nrow(parameters), .combine = 'rbind', .options.snow = opts,
                            .packages = c('reticulate','lhs', 'GeoModels', 'rrcov')
) %dopar% sensitivity.analysis(
  i=i,
  par=parameters,
  k = 5,
)
stopCluster(cl)


sensitivity.results = data.frame(sensitivity.results, parameters$beta, parameters$swapping, parameters$design, parameters$omega)
colnames(sensitivity.results) = c('FNR', 'FPR', 'F1', 'beta', 'swapping', 'design', 'omega')
# save(sensitivity.results, file = paste(root_path, '/result/supp/sensitivity-', Sys.Date(), '.Rdata', sep = ''))


############### plots of FNR, RPR, F1 score #######################
library(ggplot2)
library(dplyr)

sensitivity.plot <- data.frame(
  swapping = sensitivity.results$swapping,
  omega = factor(sensitivity.results$omega),
  beta = paste('beta=',sensitivity.results$beta, sep = ''),
  FNR = sensitivity.results$FNR,
  FPR = sensitivity.results$FPR,
  F1 = sensitivity.results$F1
)

FPR = ggplot(sensitivity.plot, aes(x = omega, y = FPR, group = omega)) + 
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

FNR = ggplot(sensitivity.plot, aes(x = omega, y = FNR, group = omega)) + 
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

F1 = ggplot(sensitivity.plot, aes(x = omega, y = F1, group = omega)) + 
  geom_boxplot() + 
  facet_grid(~swapping~beta) + 
  theme_bw() +
  theme(legend.position = 'none', 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x=expression(paste('Regularization parameter ',omega)), y='F1 scores')

gridExtra::grid.arrange(FPR, FNR, F1, ncol = 1)
