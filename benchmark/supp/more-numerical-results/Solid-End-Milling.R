# more-simulations: solid-milling


# Solid End Milling simulation

# user must specify the correct root path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
#
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))

#
library(plotly)
library(plot3D)
library(dplyr)
library(doSNOW)
library(foreach)

solid = read.table(paste(root_path, '/benchmark/supp/more-numerical-results/Solid-End-Milling.txt', sep = ''), header = TRUE)
X = solid[, 1:6]
X_scale = X
y_clean = solid$y
y = y_clean


################################################################
# outlier detection 
cov = 'Matern5_2'

# GPH
GPH.m = GPH(X=X_scale, y=y, cov = cov, iso = FALSE)
GPH.out = GPH.m$out
## GPRV
GPRV.m = GPRV(X=X, y=y, cov=cov, nugget_estim = FALSE)
GPRV.out = GPRV.m$out
## GPST
GPST.m = GPST(X, y, cov=cov,iso = FALSE)
GPST.out = GPST.m$out
## ODK
ODK.m = ODK(X=X_scale, y=y, cov = cov, iso = FALSE, phi = .1, phi_estim = TRUE, u_threshold = 1e-2, max_iter = 50)
ODK.out = ODK.m$out
#
out.sum = sort(unique(c(GPH.out, GPRV.out, GPST.out, ODK.out)))
# out.sum = ODK.out

# length(out.sum)

############# spatila plot for visualizing the outliers ############
npar = expand.grid(1:100)
pb <- txtProgressBar(max = nrow(npar), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# adapt number of kernels used
cl <- makeCluster(10)
registerDoSNOW(cl)


solid_prediction = foreach(i=1:nrow(npar), .combine = 'rbind', .options.snow = opts,
                         .packages = c('reticulate', 'georob', 'hetGP')) %dopar% {
                           
                           cat('Airfoil simulation: the ', i, '/', nrow(npar), '-th repetition.', '\n')
                           set.seed(npar[i, 1])

                           test.size = 12
                           n = nrow(X)
                           # cov = 'Matern3_2'
                           
                           test.loc = sample((1:n)[-out.sum], size = test.size)
                           X_train = X_scale[-test.loc, ]
                           y_train = y[-test.loc]
                           X_test = X_scale[test.loc, ]
                           y_test = y[test.loc]
                           
                           # with outliers 
                           temp = try({
                             CK = ODK(X=X_train, y=y_train, cov = cov, iso = FALSE, phi = 0)
                             y_pred_CK = ODK.predict(CK$model, Xnew = X_test)
                             CK.RMSE = sqrt(mean((y_test - y_pred_CK)**2))
                             cat("CK model has been executed. \n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             CK.RMSE= NA
                           }
                           
                           temp = try({
                             GPH.m = GPH(X=X_train, y=y_train, cov=cov, iso = FALSE)
                             y_pred_GPH = GPH.predict(GPH.m$model, Xnew = X_test)
                             GPH.RMSE = sqrt(mean((y_test - y_pred_GPH)**2))
                             cat("GPH model has been executed.\n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             GPH.RMSE= NA
                           }
                           
                           temp = try({
                             GPRV.m = GPRV(X=X_train, y=y_train, cov=cov, nugget_estim = FALSE)
                             y_pred_GPRV = GPRV.predict(GPRV.m$model, Xnew = X_test)
                             GPRV.RMSE = sqrt(mean((y_test - y_pred_GPRV)**2))
                             cat("GPRV model has been executed.\n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             GPRV.RMSE= NA
                           }
                           
                           temp = try({
                             GPST.m = GPST(X=X_train, y=y_train, cov = cov, iso = FALSE)
                             y_pred_GPST = GPST.predict(GPST.m$model, Xnew = X_test)
                             GPST.RMSE = sqrt(mean((y_test - y_pred_GPST)**2))
                             cat("GPST model has been executed. \n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             GPST.RMSE= NA
                           }
                           temp = try({
                             ODK.m = ODK(X=X_train, y=y_train, cov = cov, phi = .1, iso = FALSE, phi_estim = TRUE, u_threshold = 0, max_iter = 50)
                             y_pred_ODK = ODK.predict(ODK.m$model, Xnew = X_test)
                             ODK.RMSE = sqrt(mean((y_test - y_pred_ODK)^2))
                             cat("OKD model has been executed. \n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             ODK.RMSE= NA
                           } 
                           
                           RMSE = c(CK.RMSE, GPH.RMSE, GPRV.RMSE, GPST.RMSE, ODK.RMSE)
                           return(RMSE)
                         }

stopCluster(cl)

colMeans(na.omit(solid_prediction))

apply(na.omit(solid_prediction), MARGIN = 2, FUN = sd)


save(solid_prediction, file = './benchmark/supp/more-numerical-results/solid-RMSE.Rdata')

