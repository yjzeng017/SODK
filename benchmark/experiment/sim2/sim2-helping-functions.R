###########################################################
# some helping functions in Simulation 2 and they are from
# Puchhammer and Filzmoser (2023). We make some revisions.
###########################################################
library(rrcov) # for EH switching
library(lhs) # for LHS design
library(GeoModels) # for bivariate Matern model

################ DATA CONTAMINATION ###########################

switchRandom = function(n, beta, data_clean, p){
  
  # function switches observations completely random
  # n: number of observations
  # beta: percentage of switched observation
  # data_clean: uncontaminated data
  # p: number of variables of data_clean
  
  # select number of outliers to switch (round to next even number)
  outliers = sample(1:n, size = 2*round(n*beta/2), replace = F)
  text = rep(NA, n)
  
  # switch first outlier observation with last outlier and so forth
  for (i in 1:(length(outliers)/2)){
    tmp = data_clean[outliers[i], 1:p]
    data_clean[outliers[i],1:p] =  data_clean[outliers[length(outliers) - (i-1)],1:p]
    data_clean[outliers[length(outliers) - (i-1)],1:p] = tmp
    text[outliers[i]] = i
    text[outliers[length(outliers) - (i-1)]] = i
  }
  
  # create column indicating if switched or not
  data_switched = cbind(data_clean, "out" = (1:n) %in% outliers)
  
  # store information which observations were switched with each other
  label = text
  
  return(list(data = data_switched, label = label))
  
}


switchEH = function(n, beta, data_clean, p, k = 10){
  
  # switch observations that are most different according to first (global) eigenvector score
  # n: number of observations
  # beta: percentage of switched observation
  # data_clean: uncontaminated data
  # p: number of variables of data_clean
  # k: number of observations around switched observations that 
  #    are removed from switching procedure
  
  # robust PCA and first eigenvector score
  pca_rob = rrcov::PcaCov(data_clean[, 1:p], cov.control = CovControlMcd(), scale = TRUE)
  score1 = pca_rob$scores[, 1]
  
  # number of swaps
  n_swaps = round(n*beta/2)
  
  ind_toclose = c()  # collects observations that are not allowed to be switched
  ind_swapped = rep(NA, n)
  for(i in 1:n_swaps){
    
    if(i == 1) {
      vals = score1
    } else  {
      vals = score1
      vals[ind_toclose] = NA  # remove observations from switching process
    }
    
    # exchange observations with highest and lowest score values
    a = which.max(vals)
    b = which.min(vals)
    
    if(length(a) == 0 | length(b) == 0) break
    ind_swapped[a] = i
    ind_swapped[b] = i
    
    # swapping
    tmp <- data_clean[a, 1:p]
    data_clean[a, 1:p] <- data_clean[b, 1:p]
    data_clean[b, 1:p] <- tmp
    
    # removing close outliers using spatial coordinates
    xa = data_clean[a, "x"]
    ya = data_clean[a, "y"]
    xb = data_clean[b, "x"]
    yb = data_clean[b, "y"]
    
    # remove observations too close to swapped observations
    q <- data_clean[ ,c("x", "y")][c(a, b),]
    nn <- dbscan::kNN(data_clean[, c("x", "y")], k = k, query = q)
    tmp = unique(c(nn$id))
    ind_toclose = c(ind_toclose, tmp)
  }
  
  # create return values
  label = ind_swapped
  data_switched = cbind(data_clean, "out" = !is.na(ind_swapped))
  
  # check if there were enough switches (possible problem for high beta)
  if( sum(!is.na(ind_swapped))/dim(data_clean)[1] < (beta - 0.005))  {
    stop(paste("Not enough switches!", sum(!is.na(ind_swapped))/dim(data_clean)[1]))
  }
  
  return(list(data = data_switched, label = label))
}



contaminated_data_randomfields = function(seed, coords, param, beta = 0.05, type, k=10){
  # data construction with Random Fields
  
  # seed: for reprocibility
  # coords: values for both coordinates
  # param: a list consist of the parameters for bivariate matern model
  # beta: contamination level
  # type: type of switching
  # k: for type = 'EH' only.
  
  # set seeds
  set.seed(seed) # general

  # define random fields model

  # parsimonious Mathern model

  # make grid and simulate clean data
  model = GeoSim(coordx = coords, corrmodel = "Bi_matern", param = param)
  data_all = cbind(t(model$data), model$coordx, model$coordy)
  colnames(data_all) = c('variable1', 'variable2', 'x', 'y')
  
  # contaminate data
  n = dim(data_all)[1]
  data_contam = data_all
  
  
  if(beta > 0) {
    if(type == "Random"){
      switched = switchRandom(n, beta, data_all, p=2)
      label = switched$label
      data_contam = switched$data
    }
    if(type == "EH"){
      switched = switchEH(n, beta, data_all, p=2, k=k)
      label = switched$label
      data_contam = switched$data
    } 
  } else {
    data_contam = cbind(data_contam, "out" = rep(FALSE, dim(data_contam)[1]))
  }
  
  # return
  return(as.data.frame(data_contam))
}


###################################################
### EVALUATION FUNCTIONS 
#########################################################
classification = function(data, outliers){
  # function calculates false negative rate (FNR) and false positive rate (FPR)
  
  # outliers: indices of flagged outliers as numbers
  # data: needs to include "out" column indicating constructed outliers
  
  n = dim(data)[1]
  
  data_res = cbind(data, "detected" = (1:n) %in% outliers)
  tab = prop.table(table(data_res[, c("out", "detected")]), margin = 1)
  
  if(sum(dim(tab) == c(2,2)) == 2) {
    FNR = tab[2,1]
    FPR  = tab[1,2]
  }
  
  if(sum(dim(tab) == c(1,2)) == 2){ # only one true value ( = no outliers)
    FNR = 0
    FPR = tab[1,2]
  }
  
  if(sum(dim(tab) == c(2,1)) == 2){  # only one kind of classification given (no outliers)
    FPR = 0
    FNR = tab[2,1]
  }
  
  if(sum(dim(tab) == c(1,1)) == 2){ #  no outliers and non detected
    FPR = 0
    FNR = 0
  }
  
  return(c(FNR, FPR))
}

eval_runtimes = function(t1, t2, unit = "secs"){
  # returns time difference between t1 and t2 in given unit as numeric value
  origin = Sys.time()
  diff = as.POSIXct(t2, origin = origin)  - as.POSIXct(t1, origin = origin) 
  z = as.numeric(diff, units = unit)
  return(z)
}


F_measure  = function(fnr, fpr, beta, n){
  # calculates F1 score based on fnr and fpr
  
  number_out = 2*round(n*beta/2)
  number_in = n - number_out
  fn = fnr*number_out
  tp = number_out - fn
  fp = fpr*number_in
  tn = number_in - fp
  
  F1_score = 2*tp/(2*tp + fp + fn)
  
  return(F1_score)
}

##################### heatmap for contaminated data ####################
# library(ggplot2)
# library(viridis)

# Model PARAMETER SETTINGS:
# p = 2
# beta = 0.15
# type_par = c("Random", 'EH')
# grid.design = expand.grid(seq(0, 5, 0.25), seq(0, 5, 0.25))
# nu = 2.5
# rho = 0.7
# 
# param = list(mean_1 = 0, mean_2 = 0, smooth_1 = nu, smooth_2 = nu, smooth_12 = nu,
#                    scale_1 = 1, scale_2 = 1, scale_12 = 1,
#                    sill_1 = 1, sill_2 = 1, nugget_1 = 0, nugget_2 = 0, pcol = rho)
# 
# 
# grid.EH = contaminated_data_randomfields(seed = 1, coords = grid.design, param = param, beta=beta, type = 'EH', k =10)
# grid.Random = contaminated_data_randomfields(seed = 1, coords = grid.design, param = param, beta=beta, type = 'Random')
# 
# 
# # plot
# df = data.frame(
#   x1 = c(grid.Random$x, grid.EH$x),
#   x2 = c(grid.Random$y, grid.EH$y),
#   y = c(grid.Random$variable1, grid.EH$variable1),
#   swapping = rep(c('Random', "EH"), each = nrow(grid.design))
#   #nu = rep(c('v=1.5', 'v=2.5'), each = nrow(grid.design))
# )
# 
# ggplot(df, aes(x=x1, y=x2, fill=y)) +
#   geom_tile() +  
#   scale_fill_viridis(option = 'C',direction = -1) +
#   labs(fill='', x=expression(x[1]), y=expression(x[2])) +
#   facet_grid(~ swapping) + 
#   theme_bw() +
#   theme(axis.title = element_text(size = 16),
#         axis.text = element_text(size = 16),
#         panel.background = element_rect(fill = "white"),
#         plot.background = element_rect(fill = "white"),
#         axis.line = element_line(color = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         )
