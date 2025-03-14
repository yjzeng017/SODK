###########################################
# simulation 2- outlier detection 
###########################################

setwd("D:/SODK") # user must specify the correct root path
data_path = './benchmark/experiment/sim2/sim2-data-OD/'
source('./benchmark/method/comprison_methods.R')
source('./benchmark/experiment/sim2/sim2-helping-functions.R')
source('./SODK.R')
library(doSNOW)  # progress bar
library(foreach) # parallel computing

# Parameter setting
swapping = c('EH', 'Random')
design = c('Grid', 'LHS')
beta = seq(1, 15, 1)/100 # contamination level
rho = 0.7
nu = 2.5
parameters = expand.grid(1:100, beta, swapping, design)
colnames(parameters) = c("seeds", "beta", "swapping", 'design')

# main function
outlier_detection = function(i, par, cov='Matern5_2', nu=2.5, rho=0.7, k=5){
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
  
  cat('Simulation 2: The', i, '/', nrow(par), '-th replication with', 
      'nu=', nu, 'rho=', rho, 'beta=', beta, 'swapping=', as.character(swapping), 'design=',as.character(design), '\n')
  
  if(design == 'Grid') location = expand.grid(seq(0, 10, 0.5), seq(0, 10, 0.5))
  if(design == 'LHS') location = lhs::randomLHS(n=441, 2)*10
  
  # specify the parameters for bivariate Matern model
  matern_param = list(mean_1 = 0, mean_2 = 0, smooth_1 = nu, smooth_2 = nu, smooth_12 = nu,
                      scale_1 = 1, scale_2 = 1, scale_12 = 1,
                      sill_1 = 1, sill_2 = 1, nugget_1 = 0, nugget_2 = 0, pcol = rho)
  
  # generate the simulated data
  data = contaminated_data_randomfields(seed = seed, coords = location, param = matern_param, beta=beta, type = swapping, k = k)
  X = as.matrix(data[, c('x', 'y')])
  y = as.vector(data[, 'variable1'])
  n = nrow(X)
  
  ###################################
  # performance of outlier detection
  ###################################
  # GPH
  temp = try({
    b1 = as.numeric(Sys.time())
    GPH.m = GPH(X=X, y = y,cov = cov)
    b2 = as.numeric(Sys.time())
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPH.FNR = GPH.FPR = GPH.F1 = NA
  } else{
    GPH.time = eval_runtimes(b1, b2, "secs")
    GPH.FNR = classification(data = data, outliers = GPH.m$out)[1]
    GPH.FPR = classification(data = data, outliers = GPH.m$out)[2]
    GPH.F1 = F_measure(fnr = GPH.FNR, fpr = GPH.FPR, beta=beta, n = n)
    cat("GPH model has been executed and took ", GPH.time, "secs.", "\n")
  }
  
  # GPRV
  temp = try({
    c1 = as.numeric(Sys.time())
    GPRV.m = GPRV(X=X, y=y, cov = cov)
    c2 = as.numeric(Sys.time())
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPRV.FNR = GPRV.FPR = GPRV.F1 = NA
  } else {
    GPRV.time = eval_runtimes(c1, c2, "secs")
    GPRV.FNR = classification(data = data, outliers = GPRV.m$out)[1]
    GPRV.FPR = classification(data = data, outliers = GPRV.m$out)[2]
    GPRV.F1 = F_measure(fnr = GPRV.FNR, fpr = GPRV.FPR, beta=beta, n = n)
    cat("GPRV model has been executed and took ", GPRV.time, "secs.", "\n")
  }
  
  # GPST
  temp = try({
    d1 = as.numeric(Sys.time())
    GPST.m = GPST(X=X, y=y, cov = cov)
    d2 = as.numeric(Sys.time())
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    GPST.FNR =GPST.FPR = GPST.F1 = NA
  } else{
    GPST.time = eval_runtimes(d1, d2, "secs")
    GPST.FNR = classification(data = data, outliers = GPST.m$out)[1]
    GPST.FPR = classification(data = data, outliers = GPST.m$out)[2]
    GPST.F1 = F_measure(fnr = GPST.FNR, fpr = GPST.FPR, beta=beta, n = n)
    cat("GPST model has been executed and took ", GPST.time, "secs.", "\n")
  }
  
  # SODK
  temp = try({
    f1 = as.numeric(Sys.time())
    SODK.m = SODK(X=X, y=y, cov = cov, phi = 1, phi_estim = TRUE, omega = 30, u_threshold = 1e-2)
    f2 = as.numeric(Sys.time())
  }, silent = FALSE)
  if('try-error' %in% class(temp)) {
    SODK.FNR = SODK.FPR = SODK.F1 = NA
  } else{
    SODK.time = eval_runtimes(f1, f2, "secs")
    SODK.FNR = classification(data = data, outliers = SODK.m$out)[1]
    SODK.FPR = classification(data = data, outliers = SODK.m$out)[2]
    SODK.F1 = F_measure(fnr = SODK.FNR, fpr = SODK.FPR, beta=beta, n = n)
    cat("SODK model has been executed and took ", SODK.time, "secs.", "\n")
  }
  
  gc() # release the cache
  
  result = c(GPH.FNR, GPH.FPR, GPH.F1,
             GPRV.FNR, GPRV.FPR, GPRV.F1, 
             GPST.FNR, GPST.FPR, GPST.F1, 
             SODK.FNR, SODK.FPR, SODK.F1)
  return(result)
}


# set up progress bar
pb <- txtProgressBar(max = 300, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
cl = makeCluster(10)
registerDoSNOW(cl)
result = foreach(i=group_j, .combine = 'rbind', .options.snow = opts,
                   .packages = c('lhs','GeoModels','georob', 'hetGP','rrcov', 'reticulate','Matrix')
) %dopar% outlier_detection(i=i, data_path=data_path)

stopCluster(cl)

OD.result = data.frame(result, parameters$beta, parameters$swapping, parameters$design)
colnames(OD.result) = c(
    'GPRV.FNR', 'GPRV.FPR', 'GPRV.F1',
    'GPH.FNR', 'GPH.FPR', 'GPH.F1',
    'GPST.FNR', 'GPST.FPR', 'GPST.F1',
    'SODK.FNR', 'SODK.FPR', 'SODK.F1', 'beta', 'swapping','design')

##################### plot the average FPR, FNR, F1 #####################
library(dplyr)
library(ggplot2)
library(patchwork)

OD.result.mean = na.omit(OD.result) %>% 
  group_by(beta, swapping, design) %>%
  summarise(across(everything(), mean))


OD.plot <- data.frame(
  swapping = rep(OD.result.mean$swapping, 4),
  design = rep(OD.result.mean$design, 4),
  Method = rep(c("GPRV", "GPH", 'GPST', 'SODK'), each = dim(OD.result.mean)[1]),
  beta = rep(OD.result.mean$beta, 4),
  FNR = c(OD.result.mean$GPRV.FNR, OD.result.mean$GPH.FNR, OD.result.mean$GPST.FNR, OD.result.mean$SODK.FNR),
  FPR = c(OD.result.mean$GPRV.FPR, OD.result.mean$GPH.FPR, OD.result.mean$GPST.FPR, OD.result.mean$SODK.FPR),
  F1 = c(OD.result.mean$GPRV.F1, OD.result.mean$GPH.F1, OD.result.mean$GPST.F1, OD.result.mean$SODK.F1)
)

FPR = ggplot(OD.plot, aes(x = beta, y = FPR, color = Method)) +
  geom_line(aes(color = Method), linewidth=1) +
  geom_point(aes(shape = Method), size=2) +
  labs(x='',y = "FPR") +
  facet_grid(~design~swapping) +
  scale_color_manual(values = c(2, 3, 4, 5)) + 
  scale_x_continuous(
    breaks = seq(0, 0.20, by = 0.05),
    labels = c("0", "5%", "10%", "15%", "20%")
  )  + 
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))


FNR = ggplot(OD.plot, aes(x = beta, y = FNR, color = Method)) +
  geom_line(aes(color = Method), linewidth=1) +
  geom_point(aes(shape = Method), size=2) +
  labs(x = '', y = "FNR") +
  scale_x_continuous(
    breaks = seq(0, 0.20, by = 0.05),
    labels = c("0", "5%", "10%", "15%", "20%")
  ) + 
  facet_grid(~design~swapping) +
  scale_color_manual(values = c(2, 3, 4, 5)) + 
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12))


combined <- FPR + FNR & theme(legend.position = "none") & labs(x=expression(paste('Contamination level ', beta, sep = '')))
combined + plot_layout(guides = "collect")

ggsave('./figures/sim2-FPR-FNR.eps', width = 10, height = 4.0, device = cairo_ps)

F1 = ggplot(OD.plot, aes(x = beta, y = F1, color = Method)) +
  geom_line(aes(color = Method), linewidth=1) +
  geom_point(aes(shape = Method), size=2) +
  labs(x='', y = "F1 score") +
  facet_grid(~swapping+design) +
  scale_x_continuous(
    breaks = seq(0, 0.20, by = 0.05),
    labels = c("0", "5%", "10%", "15%", "20%")
  ) + 
  scale_color_manual(values = c(2, 3, 4, 5)) + 
  labs(x=expression(paste('Contamination level ', beta, sep = ''))) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        legend.position = 'bottom')

ggsave('./figures/sim2-F1.eps', width = 10, height = 4.0, device = cairo_ps)

