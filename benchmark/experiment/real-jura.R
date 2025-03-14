###############################################
# Real example: Jura dataset Atteia et al. (1994)
###############################################

# user must specify the correct root path load the R scripts
root_path = "D:/SODK"
source(paste(root_path,'/benchmark/method/comprison_methods.R', sep = ''))
source(paste(root_path,'/SODK.R', sep = ''))

# load packages
library(gstat) # for Jura data
library(ggplot2)
library(viridis)
# library(gridExtra)
library(patchwork)
library(dplyr)
library(doSNOW)  # progress bar
library(foreach) # parallel computing


############### start to analysis ####################
# load jura data
data(jura)
jura.data = rbind(jura.pred, jura.val)
# outlier identification based on all considered methods
GPH.m = GPH(X=cbind(jura.data$Xloc, jura.data$Yloc), y=log(jura.data$Pb), cov = 'Matern3_2')
GPRV.m = GPRV(X=cbind(jura.data$Xloc, jura.data$Yloc), y=log(jura.data$Pb), cov = 'exp')
GPST.m = GPST(X=cbind(jura.data$Xloc, jura.data$Yloc), y=log(jura.data$Pb), cov = 'exp', scale=0.1)
SODK.m = SODK(X=cbind(jura.data$Xloc, jura.data$Yloc), y=log(jura.data$Pb), cov = 'exp', noise=T, phi = 0.1, phi_estim = T, max_iter = 30)

# outliers
GPH.out = GPH.m$out
GPRV.out = GPRV.m$out
GPST.out = GPST.m$out
SODK.out = SODK.m$out
# which(SODK.m$phi * SODK.m$u > 0.01)
out.sum = sort(unique(c(GPH.out, GPRV.out, GPST.out, SODK.out))) # all outliers
length(out.sum)

# Maps based on spatial prediction
CK.m = SODK(X=cbind(jura.data$Xloc, jura.data$Yloc), y=log(jura.data$Pb), cov = 'exp', noise = TRUE, phi = 0)
Pb.pred= SODK.predict(CK.m, Xnew=cbind(jura.grid$Xloc, jura.grid$Yloc))

jura.map = data.frame(
  Xloc = jura.grid$Xloc,
  Yloc = jura.grid$Yloc,
  Pb = Pb.pred
)

jura.map.plot = ggplot(jura.map, aes(x=Xloc, y=Yloc)) +
  geom_raster(aes(fill = Pb)) + 
  scale_fill_viridis(option = 'C', direction = -1) + 
  geom_point(data = jura.data[out.sum,], aes(Xloc, Yloc), size = 3, shape=1, inherit.aes = FALSE) + 
  labs(fill='Log(Pb)') + 
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.border = element_rect(color = 'black', fill = NA),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 20),
  ) +
  labs(x = 'X(km)', y='Y(km)')

ggsave(paste(root_path, '/figures/real-jura-map-', Sys.Date(), '.eps', sep = ''), width = 6, height = 5, device = cairo_ps)


######################### plot the outlier in the boxplot #####################
highlight_points_GPH <- data.frame(
  x = rep(-0.3, length(GPH.out)),
  Value = log(jura.data$Pb)[GPH.out] # Points you want to highlight
)

highlight_points_GPRV <- data.frame(
  x = rep(-0.15, length(GPRV.out)),
  Value = log(jura.data$Pb)[GPRV.out] # Points you want to highlight
)
highlight_points_GPST <- data.frame(
  x = rep(0.15, length(GPST.out)),
  Value = log(jura.data$Pb)[GPST.out] # Points you want to highlight
)

highlight_points_SODK <- data.frame(
  x = rep(0.30, length(SODK.out)),
  Value = log(jura.data$Pb)[SODK.out] # Points you want to highlight
)

out.plot = ggplot(data.frame(x=rep(0, length(log(jura.data$Pb))), Pb=log(jura.data$Pb)), aes(x=x, y=Pb)) +
  geom_boxplot(width=1,outlier.shape = 1,outlier.size = 3) +
  geom_point(data=highlight_points_GPH, aes(x=x,  y=Value),
             color=2, size=3, shape=17) +
  geom_point(data=highlight_points_GPRV, aes(x=x, y=Value),
             color=3, size=3, shape=15) +
  geom_point(data=highlight_points_GPST, aes(x=x,  y=Value),
             color=4, size=3, shape=18) +
  geom_point(data=highlight_points_SODK, aes(x=x, y=Value),
             color=5, size=3, shape=7) +
  scale_x_continuous(limits = c(-1, 1)) +
  labs(x='Transformed data', y='log(Pb)') + 
  theme_bw() +
  theme(axis.title.x = element_text(size=16),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        a)

########################### apply all methods for data spliting ####################
# keep the outliers
test.size = seq(10, 100, by=10)
par = expand.grid(1:100, test.size)

for (i in 1:5) {
  group.i = ((i-1)*200+1):(i*200)
  pb <- txtProgressBar(max = length(group.i), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  # adapt number of kernels used
  cl <- makeCluster(10)
  registerDoSNOW(cl)
  RMSE.keep.out = foreach(i=group.i, .combine = 'rbind', .options.snow = opts,
                            .packages = c('reticulate', 'georob', 'hetGP','Matrix')) %dopar% {
                              
                              set.seed(par[i, 1])
                              size = par[i, 2]
                              cat('Jura: The', i, '/', nrow(par),'-th validation with', 'test sample size ', size, '\n')
                              
                              X = cbind(jura.data$Xloc, jura.data$Yloc)
                              y = log(jura.data$Pb)
                              test.loc = sort(sample((1:nrow(X))[-out.sum], size = size))
                              X.train = X[-test.loc, ]
                              y.train = y[-test.loc]
                              X.test = X[test.loc, ]
                              y.test = y[test.loc]
                              
                              # start to perform all methods
                              ## SK
                              temp = try({
                                SK.m = SODK(X=X.train, y=y.train, cov = 'exp', noise = T, phi = 0)
                                y.pred.SK = SODK.predict(SK.m, Xnew = X.test)
                                SK.rmse = sqrt(mean((y.test - y.pred.SK)^2))
                              }, silent = FALSE
                              )
                              if('try-error' %in% class(temp)) {
                                SK.rmse = NA
                              }
                              
                              ## GPH
                              temp = try({
                                GPH.m = GPH(X=X.train, y=y.train,cov = 'Matern3_2')
                                y.pred.GPH = GPH.predict(object = GPH.m, Xnew = X.test)
                                GPH.rmse = sqrt(mean((y.test - y.pred.GPH)^2))
                              }, silent = FALSE
                              )
                              if('try-error' %in% class(temp)) {
                                GPH.rmse = NA
                              }
                              
                              ## GPRV
                              temp = try({
                                GPRV.m = GPRV(X=X.train, y=y.train, cov = 'exp')
                                y.pred.GPRV = GPRV.predict(object = GPRV.m, Xnew = X.test)
                                GPRV.rmse = sqrt(mean((y.test - y.pred.GPRV)^2))
                              }, silent = FALSE
                              )
                              if('try-error' %in% class(temp)) {
                                GPRV.rmse = NA
                              }
                              
                              ## GPST
                              temp = try({
                                GPST.m = GPST(X=X.train, y=y.train, cov = 'exp', scale = 0.1)
                                y.pred.GPST = GPST.predict(GPST.m, Xnew = X.test)
                                GPST.rmse = sqrt(mean((y.test - y.pred.GPST)^2))
                              }, silent = FALSE
                              )
                              if('try-error' %in% class(temp)) {
                                GPST.rmse = NA
                              }
                              
                              ## SODK
                              temp = try({
                                SODK.m = SODK(X=X.train, y=y.train, cov = 'exp', noise = T, phi = 0.1, phi_estim = T, threshold = F, max_iter = 30)
                                y.pred.SODK = SODK.predict(SODK.m, Xnew = X.test)
                                SODK.rmse = sqrt(mean((y.test - y.pred.SODK)^2))
                              }, silent = FALSE
                              )
                              if('try-error' %in% class(temp)) {
                                SODK.rmse = NA
                              }
                              SODK.m$noise_var
                              SODK.m$out
                              SODK.m$phi
                              
                              gc()
                              return(c(SK.rmse, GPH.rmse, GPRV.rmse, GPST.rmse, SODK.rmse))
                            }
  stopCluster(cl)
  save(RMSE.keep.out, file = paste('D:/SODK/result/jura-RMSE-keep-outlier-', i,'.RData', sep = ''))
}

# summarize the RMSE results
temp = data.frame()
for (i in 1:5) {
  temp.i = get(load(file= paste('D:/SODK/result/jura-RMSE-keep-outlier-', i,'.RData', sep = '')))
  temp = rbind(temp, temp.i)
}

RMSE = data.frame(temp, par[, 2])
colnames(RMSE.keep.out) = c('CK','GPH','GPRV', 'GPST', 'SODK', 'test.size')
meanRMSE.keep.out = na.omit(RMSE.keep.out) %>% group_by(test.size) %>%
  summarise(across(everything(), mean))

RMSE.df = data.frame(
  test.size = rep(meanRMSE.keep.out$test.size, 5),
  Scenario1.RMSE = as.vector(as.matrix(meanRMSE.keep.out[, -which(names(meanRMSE.keep.out) == 'test.size')])),
  Method = rep(c('CK', 'GPH', 'GPRV', 'GPST', 'SODK'), each=length(meanRMSE.keep.out$test.size))
)
RMSE.plot = ggplot(RMSE.df, aes(x=test.size, y=Scenario1.RMSE, color = Method)) +
  geom_line(linewidth=0.8) +
  geom_point(aes(shape = Method), size=3) + 
  labs(x = 'Number of prediction samples', y ='RMSE') + 
  scale_color_manual(values = c(1, 2, 3, 4, 5)) + 
  scale_shape_manual(values = c(19, 17, 15, 18, 7)) + 
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) 

#######################################################
combined <- out.plot + Scenario1.plot  & theme()
combined + plot_layout(guides = "collect")

