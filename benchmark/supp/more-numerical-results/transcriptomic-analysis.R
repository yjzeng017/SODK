# transcriptomic examples


# user must specify the correct root path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
#
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))
# 
library(dplyr)
# 
library(foreach)
library(doSNOW)


################# analysis ############3

# load data
load(file = './benchmark/supp/more-numerical-results/MOB_patterns.Rdata')
load(file = './benchmark/supp/more-numerical-results/BC_patterns.Rdata')


# adapt number of kernels used
cov = 'gauss'
patterns = c('Pattern1', 'Pattern2')
examples = c('MOB', 'BC')
npar = expand.grid(1:100, patterns, examples)
pb <- txtProgressBar(max = nrow(npar), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
# adapt number of kernels used

cl = makeCluster(10)
registerDoSNOW(cl)


transcriptomic_prediction = foreach(i=1:nrow(npar), .combine = 'rbind', .options.snow = opts,
                         .packages = c('reticulate', 'georob', 'hetGP')) %dopar% {
                           
                           cat('transcriptomic examples: the ', i, '/', nrow(npar), '-th repetition.', '\n')
                           set.seed(npar[i, 1])
                           pattern = npar[i, 2]
                           example = npar[i, 3]
                           
                           
                           if(example == 'MOB') {
                             X = MOB_patterns$coordinate
                             if(pattern == 'Pattern1') y = MOB_patterns$patterns[, 1]
                             if(pattern == 'Pattern2') y = MOB_patterns$patterns[, 2]
                           }
                             
                           if(example == 'BC') {
                             X = BC_patterns$coordinate
                             if(pattern == 'Pattern1') y = BC_patterns$patterns[, 1]
                             if(pattern == 'Pattern2') y = BC_patterns$patterns[, 2]
                           }
                           
                           
                           n = nrow(X)
                           test.size = as.integer(n * 0.2)
                           test.loc = sample(1:n, size = test.size)
                           X_train = X[-test.loc, ]
                           y_train = y[-test.loc]
                           
                           out = sample(1:length(y_train), 5)
                           y_train[out] = y_train[out] + 2*sd(y)
                             
                           X_test = X[test.loc, ]
                           y_test = y[test.loc]
                           
                           # with outliers 
                           temp = try({
                             CK = ODK(X=X_train, y=y_train, cov = cov, iso = FALSE, phi = 0)
                             y_pred_CK = ODK.predict(CK$model, Xnew = X_test)
                             CK.RMSE = sqrt(mean((y_test - y_pred_CK)^2))
                             cat("CK model has been executed. \n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             CK.RMSE= NA
                           }
                           
                           temp = try({
                             GPH.m = GPH(X=X_train, y=y_train, cov=cov, iso = FALSE)
                             y_pred_GPH = GPH.predict(GPH.m$model, Xnew = X_test)
                             GPH.RMSE = sqrt(mean((y_test - y_pred_GPH)^2))
                             cat("GPH model has been executed.\n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             GPH.RMSE= NA
                           }
                           
                           temp = try({
                             GPRV.m = GPRV(X=X_train, y=y_train, cov=cov)
                             y_pred_GPRV = GPRV.predict(GPRV.m$model, Xnew = X_test)
                             GPRV.RMSE = sqrt(mean((y_test - y_pred_GPRV)^2))
                             cat("GPRV model has been executed.\n")
                           }, silent = FALSE)
                           if('try-error' %in% class(temp)) {
                             GPRV.RMSE= NA
                           }
                           
                           temp = try({
                             GPST.m = GPST(X=X_train, y=y_train, cov = cov, iso = FALSE)
                             y_pred_GPST = GPST.predict(GPST.m$model, Xnew = X_test)
                             GPST.RMSE = sqrt(mean((y_test - y_pred_GPST)^2))
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

transcriptomic_prediction = data.frame(transcriptomic_prediction, npar$Var2, npar$Var3)
colnames(transcriptomic_prediction) = c('CK', 'GPH', 'GPRV', 'GPST', 'ODK', 'Pattern', 'Example')
save(transcriptomic_prediction, file = paste(root_path, '/result/supp/supp-transcriptomic-RMSE-', Sys.Date(), '.Rdata', sep = ''))


na.omit(transcriptomic_prediction) %>% group_by(Pattern, Example) %>%
  summarise(across(everything(), mean)) %>%  group_by(Pattern, Example) %>% summarise(across(everything(), function(x){round(x, 2)}))

na.omit(transcriptomic_prediction) %>% group_by(Pattern, Example) %>%
  summarise(across(everything(), sd)) %>%  group_by(Pattern, Example) %>% summarise(across(everything(), function(x){round(x, 2)}))
