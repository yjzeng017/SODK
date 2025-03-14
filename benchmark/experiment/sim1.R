######################## Simulation 2: a six-hump camel example ####################
# NOTE: due to randomness, the plots may be not always the same,
#       but the SODK basically captures the true patterns.


# packages for visulization
library(plot3D)
library(ggplot2)
library(geomtextpath)
library(viridis)
# for LHS design
library(lhs)

#NOTE: user should specify the root path
root_path = "C:/Users/Youjie Zeng/Desktop/SODK"
source(paste(root_path, '/SODK.R', sep = ''))

# six-hump camel function
six.hump = function(x1, x2){
  return(x1**2 * (4 - 2.1*x1**2 + x1**4/3) + x1*x2 + x2**2 * (-4 + 4*x2**2))
}


############## start Simulation 1 ################

# visulize the six-hump camel function
num = 50
M <- mesh(seq(-2, 2, length.out = num), seq(-1, 1, length.out = num))
z = six.hump(M$x, M$y)

persp3D(x=seq(-2, 2, length.out = num), y=seq(-1, 1, length.out = num), z=z,
      xlim = c(-2, 2), ylim = c(-1, 1), zlim=c(-1, 6),
      xlab = '', ylab = '', zlab='',
     box=TRUE, nticks=5, ticktype="detailed",
      border = 'black', axes=TRUE, col = 'white', lwd = 0.1,
      colkey = FALSE, xpd=TRUE, mar = c(0, 2, 0, 2),
     main='Six-hump camel function')
text3D(0, -1.8, 0, expression(x[1]), add = TRUE, cex = 1.5)
text3D(3, -0.5, 0, expression(x[2]), add = TRUE, cex= 1.5)
text3D(-2, -1, 7.5, expression(f(x[1],x[2])), add = TRUE, cex = 1.5)


# contour plot and heatmap of patterns
set.seed(12)

num = 100
X = randomLHS(num, 2)
X[, 1] = X[, 1]*4 - 2
X[, 2] = X[, 2]*2 - 1
X = X[order(X[, 1]), ]
y_clean = six.hump(X[, 1], X[, 2])
epsilon = rep(0, num)
out = seq(10, 100, 20) # location of outliers
epsilon[out] = 1
epsilon = epsilon * sd(y_clean) # shifting 5 values of y_clean with one standard deviation of y_clean
y = y_clean + epsilon # y is contaminated

# grid points for prediction
Xnew =expand.grid(seq(-2, 2, length.out = 50), seq(-1, 1, length.out = 50))
z.true = six.hump(Xnew[, 1], Xnew[, 2])

# SODK prediction at Xnew
kriging.SODK = SODK(X=X, y=y, cov='sq_exp', iso = FALSE, omega=30, phi=1, phi_estim = TRUE)
z.SODK = SODK.predict(kriging.SODK, Xnew = Xnew)
# CK prediction at Xnew
kriging.CK = SODK(X=X, y=y, cov='sq_exp',iso = FALSE, phi=0)
z.CK =  SODK.predict(kriging.CK, Xnew = Xnew)


# plot the patterns from CK prediction, SODK prediction and True function
dummy = data.frame(x=rep(X[out, 1], 3), 
                   y=rep(X[out, 2], 3), 
                   group = rep(c("True", 'SODK',  "CK"), each = 5))

df.pattern = data.frame(
  x = rep(Xnew[, 1], 3),
  y = rep(Xnew[, 2], 3),
  z = c(z.true, z.SODK, z.CK),
  group = factor(rep(c("True", 'SODK',  "CK"), each = 50^2))
)

pattern_plot = ggplot(df.pattern, aes(x=x, y=y, fill=z)) +
  geom_tile() +  
  geom_contour(aes(z = z), color = 'white', linewidth = 0.8) + 
  scale_fill_viridis(option = 'C', direction = -1) + 
  labs(fill='Value', x=expression(x[1]), y=expression(x[2])) + 
  facet_wrap(~ group) +  
  geom_point(data = dummy[dummy$group %in% c('SODK', "CK"), ], 
             aes(x, y), color = 'black', size = 4, ,inherit.aes = FALSE) + 
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste(root_path, "/figures/sim1-patterns-", Sys.Date(), ".eps", sep = ''), width = 10, height = 4.0, device = cairo_ps)


################# robust analysis ##########################

# function to computing kriging weights
kriging.weights = function(K, Knew){
  invK = solve(K)
  n = nrow(K)
  m = nrow(Knew)
  ones.n = matrix(rep(1, n), ncol = 1)
  ones.m = matrix(rep(1, m), ncol = 1)
  tmp1 = Knew %*% invK
  tmp2 =  t(ones.n) %*% invK %*% ones.n
  tmp3 = ones.n %*% t(ones.n) %*% invK / as.numeric(tmp2)
  tmp4 = tmp1 %*% tmp3
  tmp5 = ones.m %*% t(ones.n) %*% invK / as.numeric(tmp2)
  weights = tmp1 - tmp4 + tmp5
  return(weights)
}


gpflow = import('gpflow')
K_CK = gpflow$utilities$add_noise_cov(kriging.CK$model$kernel$K(X=as.matrix(X)), kriging.CK$model$likelihood$variance$numpy()[, 1])
K_new_CK = kriging.CK$model$kernel$K(X=as.matrix(X), X2=as.matrix(Xnew))

K_SODK = gpflow$utilities$add_noise_cov(kriging.SODK$model$kernel$K(X=as.matrix(X)), kriging.SODK$model$likelihood$variance$numpy()[, 1])
K_new_SODK = kriging.SODK$model$kernel$K(X=as.matrix(X), X2=as.matrix(Xnew))


w.CK = kriging.weights(K_CK$numpy()[,], t(K_new_CK$numpy()[,]))
w.SODK = kriging.weights(K_SODK$numpy()[,], t(K_new_SODK$numpy()[,]))
weighted.bias.CK = w.CK %*% epsilon / max(abs(epsilon))
weighted.bias.SODK = w.SODK %*% epsilon / max(abs(epsilon))

# plot the contour of deviation from CK prediction, SODK prediction
df.deviation = data.frame(
  w = c(weighted.bias.CK, weighted.bias.SODK),
  x = rep(Xnew[, 1], 2),
  y = rep(Xnew[, 2], 2),
  group = rep(c("CK", 'SODK'), each = 50^2)
)
dummy = data.frame(x=rep(X[out, 1], 2), 
                   y=rep(X[out, 2], 2), 
                   group = rep(c("CK", 'SODK'), each = 5))

deviation_contour = ggplot(df.deviation, aes(x=x, y=y, z=w)) +
  geom_textcontour(bins=6) + 
  labs(fill='', x=expression(x[1]), y=expression(x[2])) + 
  facet_wrap(~ group) +  
  geom_point(data = dummy, aes(x, y), color = 'black', size = 2,inherit.aes = FALSE) + 
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 16),
        plot.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(paste(root_path, "/figures/sim1-deviation-", Sys.Date(), ".eps", sep = ''), width = 10, height = 5, device = cairo_ps)

