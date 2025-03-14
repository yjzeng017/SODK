##############################################
# Real example: Airfoil simulation
##############################################

# user must specify the correct root path load the R scripts
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))

# load packages
# library(patchwork)
library(plotly)
library(plot3D)
library(dplyr)
library(doSNOW)  # progress bar
library(foreach) # parallel computing


###################### start to analysis #######################
air = read.table(paste(root_path,'/benchmark/experiment/Airfoil.txt', sep = ''), header = TRUE)
X = air[, 1:4]
y_clean = air$gliding_ratio
y = y_clean
X_scale = X
X_scale[, 1] = X_scale[, 1]/0.095
X_scale[, 2] = (X_scale[, 2] - 0.30)/(0.50 - 0.30)
X_scale[, 3] = (X_scale[, 3] - 0.06)/(0.30 - 0.06)
X_scale[, 4] = X_scale[, 4]/10

# outlier identification  
cov = 'Matern3_2'
GPH.m = GPH(X=X_scale, y=y, cov = cov, iso = FALSE)
GPRV.m = GPRV(X=X, y=y, cov=cov, nugget_estim = FALSE)
GPST.m = GPST(X, y, cov=cov,iso = FALSE)
SODK.m = SODK(X=X_scale, y=y, cov = cov, iso = FALSE, phi = .01, phi_estim = TRUE, u_threshold = 1e-2, max_iter = 50)
# 
GPH.out = GPH.m$out
GPRV.out = GPRV.m$out
GPST.out = GPST.m$out
SODK.out = SODK.m$out
out.sum = sort(unique(c(GPH.out, GPRV.out, GPST.out, SODK.out)))


# spatila plot of outliers
cate = rep('Normal', nrow(X))
cate[out.sum] = 'Outlier'
shape = rep(0, nrow(X))
shape[out.sum] = 1


df.out = data.frame(
  x1 = X[, 1],
  x2 = X[, 2],
  x3 = X[, 3],
  x4 = X[, 4],
  category = factor(cate),
  shape = factor(shape)
)

outlier.plot = plot_ly(df.out, x = ~x1, y = ~x3, z = ~x4,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 7),
        color = ~category, 
        # colors = c('#1f77b4', '#ff7f0e'),
        symbol = ~category,
        symbols = c('circle-open', 'diamond')
) %>%
  layout(legend = list(x=0.85, y = 0.85), font = list(size=18),
         scene=list(xaxis=list(title = 'M'),
                    yaxis=list(title = 'T'),
                    zaxis=list(title = 'alpha')))

# NOTE: the .PNG image must be saved in the Viewer of RStudio


##################### prediction ######################
test.size = seq(20, 60, 20)
par = expand.grid(1:100, test.size)
pb <- txtProgressBar(max = nrow(par), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(15)
registerDoSNOW(cl)


air.prediction = foreach(i=1:nrow(par), .combine = 'rbind', .options.snow = opts,
                     .packages = c('reticulate', 'georob', 'hetGP')) %dopar% {
                       cat('Airfoil simulation: the ', i, '/', nrow(par), '-th repetition.', '\n')
                       set.seed(par[i, 1])
                       test.size = par[i, 2]
                       n = nrow(X)
                       test.loc = sample((1:n)[-out.sum], size = test.size)
                       X_test = X_scale[test.loc, ]
                       y_test = y[test.loc]
                      
                       # keep outliers 
                       test.loc = sample((1:n)[-out.sum], size = test.size)
                       X_train = X_scale[-test.loc, ]
                       y_train = y[-test.loc]
                       X_test = X_scale[test.loc, ]
                       y_test = y[test.loc]
                       
                       # CK 
                       temp = try({
                         CK = SODK(X=X_train, y=y_train, cov = cov, iso = FALSE, phi = 0)
                         y.pred.CK = SODK.predict(CK, Xnew = X_test)
                         CK.RMSE = sqrt(mean((y_test - y.pred.CK)**2))
                         cat("CK is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CK.RMSE= NA
                       }
                       
                       # GPH
                       temp = try({
                         GPH.m = GPH(X=X_train, y=y_train, cov=cov, iso = FALSE)
                         y.pred.GPH = GPH.predict(GPH.m, Xnew = X_test)
                         GPH.RMSE = sqrt(mean((y_test - y.pred.GPH)^2))
                         cat("GPH is done.\n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         GPH.RMSE= NA
                       }
                       
                       # GPRV
                       temp = try({
                         GPRV.m = GPRV(X=X_train, y=y_train, cov=cov, nugget_estim = FALSE)
                         y.pred.GPRV = GPRV.predict(GPRV.m, Xnew = X_test)
                         GPRV.RMSE = sqrt(mean((y_test - y.pred.GPRV)**2))
                         cat("GPRV id done.\n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         GPRV.RMSE= NA
                       }
                       
                       # GPST
                       temp = try({
                         GPST.m = GPST(X=X_train, y=y_train, cov = cov, iso = FALSE)
                         y.pred.GPST = GPST.predict(GPST.m, Xnew = X_test)
                         GPST.RMSE = sqrt(mean((y_test - y.pred.GPST)**2))
                         cat("GPST is done.\n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         GPST.RMSE= NA
                       }
                       
                       # SODK
                       temp = try({
                         SODK.m = SODK(X=X_train, y=y_train, cov = cov, phi = .001, iso = FALSE, phi_estim = TRUE, u_threshold = 0, max_iter = 50)
                         y.pred.SODK = SODK.predict(SODK.m, Xnew = X_test)
                         SODK.RMSE = sqrt(mean((y_test - y.pred.SODK)**2))
                         cat("OKD is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         SODK.RMSE= NA
                       } 
                       
                       # remove outliers
                       ## CKD-GPH
                       temp = try({
                         CKD.GPH = SODK(X=X_scale[-c(test.loc, GPH.out), ], y=y[-c(test.loc, GPH.out)], cov = cov, iso = FALSE, phi = 0)
                         y.pred.CKD.GPH = SODK.predict(CKD.GPH, Xnew = X_test)
                         CKD.GPH.RMSE = sqrt(mean((y_test - y.pred.CKD.GPH)^2))
                         cat("CKD-GPH is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CKD.GPH.RMSE= NA
                       }
                       
                       ## CKD-GPRV
                       temp = try({
                         CKD.GPRV = SODK(X=X_scale[-c(test.loc, GPRV.out), ], y=y[-c(test.loc, GPRV.out)], cov = cov, iso = FALSE, phi = 0)
                         y.pred.CKD.GPRV = SODK.predict(CKD.GPRV, Xnew = X_test)
                         CKD.GPRV.RMSE = sqrt(mean((y_test - y.pred.CKD.GPRV)^2))
                         cat("CKD-GPRV is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CKD.GPRV.RMSE= NA
                       }
                       
                       ## CKD-GPST
                       temp = try({
                         CKD.GPST = SODK(X=X_scale[-c(test.loc, GPST.out), ], y=y[-c(test.loc, GPST.out)], cov = cov, iso = FALSE, phi = 0)
                         y.pred.CKD.GPST = SODK.predict(CKD.GPST, Xnew = X_test)
                         CKD.GPST.RMSE = sqrt(mean((y_test - y.pred.CKD.GPST)^2))
                         cat("CKD-GPST is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CKD.GPST.RMSE= NA
                       }
                       
                       ## CKD-SODK
                       temp = try({
                         CKD.SODK = SODK(X=X_scale[-c(test.loc, SODK.out), ], y=y[-c(test.loc, SODK.out)], cov = cov, iso = FALSE, phi = 0)
                         y.pred.CKD.SODK = SODK.predict(CKD.SODK, Xnew = X_test)
                         CKD.SODK.RMSE = sqrt(mean((y_test - y.pred.CKD.SODK)^2))
                         cat("CKD-SODK is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CKD.SODK.RMSE= NA
                       }
                       
                       # CKD-ALL
                       temp = try({
                         CKD.ALL = SODK(X=X_scale[-c(test.loc, out.sum), ], y=y[-c(test.loc, out.sum)], cov =cov, iso = FALSE, phi = 0)
                         y.pred.CKD.ALL = SODK.predict(CKD.ALL, Xnew = X_test)
                         CKD.ALL.RMSE = sqrt(mean((y_test - y.pred.CKD.ALL)^2))
                         cat("CKD-ALL is done. \n")
                       }, silent = FALSE)
                       if('try-error' %in% class(temp)) {
                         CKD.ALL.RMSE= NA
                       }
                       
                       RMSE.D = c(CKD.GPH.RMSE, CKD.GPRV.RMSE, CKD.GPST.RMSE, CKD.SODK.RMSE, CKD.ALL.RMSE)
                       RMSE = c(CK.RMSE, GPH.RMSE, GPRV.RMSE, GPST.RMSE, SODK.RMSE)
                       
                       return(c(RMSE, RMSE.D))
                     }

stopCluster(cl)

air.prediction = data.frame(air.prediction, test.size=par[, 2])
colnames(air.prediction) = c('CK', 'GPH', 'GPRV', 'GPST', 'SODK', 'CKD-GPH', 'CKD-GPRV', 'CKD-GPST', 'CKD-SODK', 'CKD-ALL', 'test.size')

save(air.prediction, file = paste(root_path, '/result/Airfoil-RMSE-', Sys.Date(), '.Rdata', sep = ''))

meanRMSE = na.omit(air.prediction) %>% group_by(test.size) %>%
   summarise(across(everything(), mean)) %>%  group_by(test.size) %>% summarise(across(everything(), function(x){round(x, 2)}))
                                                                                                                                                            
sdRMSE = na.omit(air.prediction) %>% group_by(test.size) %>%
  summarise(across(everything(), sd)) %>%  group_by(test.size) %>% summarise(across(everything(), function(x){round(x, 2)}))

