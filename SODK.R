library(reticulate) # for interface to 'Python'
library(Matrix)

# NOTE: The proposed ODK model is based on the Python package 'gpflow'. Thus, 
#       the user should install the Python environment on a PC, and then install
#       the Python package 'gpflow'.

#       Before carrying out ODK, the user should activate the Python environment 
#       in R by using R package 'reticulate'. 


# NOTE: User should specify the path of Python environment
python_path = "D:/Anaconda3/python.exe" 
use_python(python_path)

# check the Python environment and the Python package 'gpflow'
if(!py_module_available('gpflow')){
  stop("Error message: the GPflow module should be available. \n")
}

### Main Function for ODK #######################################################################
# main function
SODK = function(X, y, cov='Matern5_2', iso=TRUE, noise=FALSE, omega=30, phi=1, phi_estim = FALSE, max_iter=25, change_between_iter = 1e-8, u_threshold = 1e-2, threshold=FALSE){
  # X: design matrix 
  # y: response
  # cov: covariance functions: 
  #     'sq_exp': square exponential covariance
  #     'gauss': Gaussian covariance
  #     'exp': exponential covariance
  #     'Matern3_2': Matern class with smooth parameter nu=3/2
  #     'Matern3_2': Matern class with smooth parameter nu=5/2
  # noise: a logical scalar: TRUE for regression kriging; FALSE for Interpolation kriging
  # iso: a logical scalar: TRUE for isotropic covariance while FALSE for anisotropic covariance.
  # scale_y: logical scalar: if TRUE, the response y is scaled as y/sd(y);
  #                         otherwise, use y
  # omega: the regularization parameter. We suggest omega = 30.
  # phi, phi_estim: if phi_estim is TRUE, phi is used as the initial value of dispersion parameter;
  #                 if phi_estim is FALSE, phi is fixed.
  #                 Particularly, phi=0 gives the conventional kriging.
  # u_threshold: the thresholding rule for obtaining the zero solution for u.
  #              If phi*u < u_threshold * (sd(y))^2, we then set u=0.
  # threshold: a logical scalar: TRUE for thresholding during the iterations; 
  #                             FALSE for thresholding as last iteration.
  #                             (Empitically, it is found that they are similar.)
  #         we consider K + phi*U + jitter * I instead, where jitter is a small positive value.
  # max_iter, change_between_iter: stopping rule
  
  # return: result list
  
  # import the gpflow from Python.
  gpflow = reticulate::import('gpflow')
  
  X= as.matrix(X)
  y = as.vector(y)
  n = dim(X)[1]
  d = dim(X)[2]
  y_scale= y/sd(y) # scale y to have unit variance
  u0 = rep(0, n)
  
  # initialize covariance
  if(iso) {lengthscales = 1} else {lengthscales = as.vector(rep(1, d))}
  
  if(cov == 'sq_exp') kernel = gpflow$kernels$SquaredExponential()
  if(cov == 'gauss') kernel = gpflow$kernels$RBF(lengthscales=lengthscales)
  if(cov == 'exp') kernel = gpflow$kernels$Exponential(lengthscales=lengthscales)
  if(cov == 'Matern3_2') kernel = gpflow$kernels$Matern32(lengthscales=lengthscales)
  if(cov == 'Matern5_2') kernel = gpflow$kernels$Matern52(lengthscales=lengthscales)
  
  #################### if phi=0, it gives a CK ####################
  
  if(phi == 0) {
    u1 = u0
    if(noise) noise_var = 0.1
    if(!noise) noise_var = 1e-5
    
    model = gpflow$models$GPR(
      data=list(X, as.matrix(y_scale)),
      kernel=kernel,
      mean_function=gpflow$mean_functions$Constant(1),
      likelihood=gpflow$likelihoods$Gaussian(variance=noise_var)
    )
    # make interpolation or regression
    if(!noise) gpflow$set_trainable(model$likelihood, FALSE)
    opt = gpflow$optimizers$Scipy()
    opt$minimize(model$training_loss, model$trainable_variables)
    noise_var = model$likelihood$variance$numpy()[1]
    if(!noise) {cat('SK for interpolation is done. \n')} else {cat('SK for regression is done. \n')}
  }
  
  ######### if phi > 0 and noise=FALSE, it gives SODK for interpolation #######
  
  # initialize un
  u1 = rep(1, n)
  u_gap = c()
  theta_gap = c()
  
  if(phi > 0 & !noise) {
    noise_var = 1e-5
    theta1 = unlist(lapply(kernel$parameters, function(x) x$numpy()))
    iter = 1
    while (iter <= 15 | (mean((u1-u0)^2)/(mean(u0^2) + 1e-8) > change_between_iter & iter < max_iter + 1)) {
      iter = iter + 1
      model = gpflow$models$GPR(
        data=list(X, as.matrix(y_scale)),
        kernel=kernel,
        mean_function=gpflow$mean_functions$Constant(1),
        likelihood=gpflow$likelihoods$Gaussian(variance = noise_var + phi * as.matrix(u1)) # 1e-5 is a jitter term for stable computation
      ) 
      
      gpflow$set_trainable(model$likelihood, FALSE)
      opt = gpflow$optimizers$Scipy()
      opt$minimize(model$training_loss, model$trainable_variables)
      
      # update theta
      theta0 = theta1
      theta1 = unlist(lapply(model$kernel$parameters, function(x) x$numpy()))
      
      # update phi (if phi_estim=TRUE)
      ## Here, we only update phi at first 10 iterations to reduce computation time
      if(phi_estim) {
        mu = model$mean_function(X)$numpy()[, 1]
        K = model$kernel$K(X, X)$numpy()
        opt = optim(par=phi, fn=phi_loss, gr=gradient_phi, method='L-BFGS-B', lower = 0.1, upper = 10, X=X, y=y_scale, mu=mu, K=K, u=u1, noise_var = 1e-5)
        phi = opt$par
        }
      
      # update u
      f_hat = model$predict_y(X)[[1]]$numpy()[, 1]
      epsilon_hat = y_scale - f_hat
      u0 = u1
      u1 = (2 - omega) / 4 + sqrt(
        (2 - omega) ** 2 + 8 * omega * (epsilon_hat ** 2) / phi) / 4
      
      # save the change of un and theta between two iterations
      theta_gap = c(theta_gap, sqrt(sum((theta0 - theta1)**2)))
      u_gap = c(u_gap, sqrt(sum((phi*u1 - phi*u0)^2)))
      
      # threshold
      if(threshold & iter >= 15){
        u1[phi * u1 < u_threshold] = 0
      }
    }
    cat('SODK for interpolation is done. \n')
  }
  
  ######### if phi > 0 and noise==TRUE, it gives SODK for regression #########
  
  if(phi >0 & noise){
    noise_var = 0.01
    iter = 1
    while (iter <= 15 | (mean((u1-u0)^2)/(mean(u0^2) + 1e-8) > change_between_iter & iter < max_iter + 1)) {
      iter = iter + 1
      # update theta
      model = gpflow$models$GPR(
        data=list(X, as.matrix(y_scale)),
        kernel=kernel + gpflow$kernels$White(noise_var),
        mean_function=gpflow$mean_functions$Constant(1),
        likelihood=gpflow$likelihoods$Gaussian(variance = 1e-5 + phi* as.matrix(u1)) # 1e-5 is a jitter term for stable computation
      ) 
      theta1 = unlist(lapply(model$kernel$parameters, function(x) x$numpy()))
      gpflow$set_trainable(model$likelihood, FALSE)
      opt = gpflow$optimizers$Scipy()
      opt$minimize(model$training_loss, model$trainable_variables)
      
      noise_var = model$kernel$kernels[1]$variance$numpy()[1]
      # cat('The ', iter, 'th: noise_var = ', noise_var, '\n')
      
      theta0 = theta1
      theta1 = unlist(lapply(model$kernel$parameters, function(x) x$numpy()))

      # update phi (if phi_estim=TRUE)
      # Here, we only update phi at first 10 iterations to reduce computation time
      if(phi_estim & iter < 10){
        mu = model$mean_function(X)$numpy()[, 1]
        K = model$kernel$K(X, X)$numpy()
        opt = optim(par=phi, fn=phi_loss, gr=gradient_phi, method='L-BFGS-B', lower = 0.1, upper = 10, X=X, y=y_scale, mu=mu, K=K, u=u1, noise_var = noise_var)
        phi = opt$par
      }

      # if(phi_estim & iter < 15){
      #   model_phi = gpflow$models$GPR(
      #     data=list(X, as.matrix(y_scale/sqrt(phi))),
      #     kernel=kernel + gpflow$kernels$White(0.1),
      #     mean_function=gpflow$mean_functions$Constant(1),
      #     likelihood=gpflow$likelihoods$Gaussian(variance=1e-5 + as.matrix(u1))
      #   )
      #   gpflow$set_trainable(model_phi$likelihood, FALSE)
      #   opt = gpflow$optimizers$Scipy()
      #   opt$minimize(model_phi$training_loss, model_phi$trainable_variables)
      #   tau = model$kernel$kernels[1]$variance$numpy()[1]
      #   phi = as.numeric(noise_var/tau)
      # }
      
      cat('The ', iter, 'th optimization: ', noise_var, '\n')
      cat('phi =  ', phi, '\n')
      cat('noise variance =  ', noise_var, '\n')
      
      # update u
      tf = import('tensorflow')
      mu = model$mean_function(X)$numpy()[, 1]
      K = model$kernel$K(X)$numpy()
      Q_tilde = K +  phi * diag(u1, n, n) + diag(noise_var, n, n) 
      invQ_tilde = tf$linalg$inv(Q_tilde)$numpy()
      f_tilde = K %*% invQ_tilde %*% (y_scale - mu)
    
      Q = K + diag(noise_var, n, n)
      invQ = tf$linalg$inv(Q)$numpy()
      f_hat = K %*% invQ %*% (y_scale - mu)
      epsilon_hat = (f_hat - f_tilde)[, 1]
      u0 = u1
      u1 = (2 - omega) / 4 + sqrt(
        (2 - omega) ** 2 + 8 * omega * (epsilon_hat ** 2) / phi) / 4
      
      # save the change of un and theta between two iterations
      theta_gap = c(theta_gap, sqrt(sum((theta0 - theta1)**2)))
      u_gap = c(u_gap, sqrt(sum((phi*u1 - phi*u0)^2)))
      
      # threshold
      if(threshold & iter >= 15){
        u1[phi * u1 < u_threshold] = 0
      }
    }
    cat('SODK for regression is done. \n')
  }
  if(phi > 0) out = which(phi * u1 > u_threshold)
  if(phi == 0) out=NA
  return(list(model=model, X=X, y=y, noise=noise, noise_var = noise_var, u=u1, phi=phi, out = out, u_gap = u_gap, theta_gap=theta_gap))
}

# prediction function
SODK.predict = function(object, Xnew){
  # object: result of ODK
  # Xnew: new inputs
  # return: prediction at Xnew
  tf = import('tensorflow')
  gpflow = import('gpflow')
  model = object$model
  Xnew = as.matrix(Xnew)
  X = model$data[0]
  y_scale = model$data[1]
  y = object$y
  u = object$u
  phi = object$phi
  noise = object$noise
  noise_var = object$noise_var
  
  mu = model$mean_function(X)
  K_new = model$kernel$K(X=X, X2=Xnew)
  if(noise) { K = model$kernel$K(X, X); likelihood_variance = noise_var + phi*u}
  if(!noise) {K = model$kernel$K(X, X); likelihood_variance = 1e-5 + phi*u}
  Q = gpflow$utilities$add_noise_cov(K, likelihood_variance = likelihood_variance)
  invQ = tf$linalg$inv(Q)
  alpha = tf$linalg$matmul(invQ, y_scale - mu)
  f_pred = model$mean_function(Xnew) + tf$linalg$matmul(K_new, alpha, transpose_a=TRUE)
  y_pred = sd(y) * f_pred$numpy()[, 1] 
  return(y_pred)
}



########## subfunctions ########
phi_loss = function(phi, X, y, mu, K, u, noise_var) {
  # phi = exp(phi)
  n = dim(X)[1]
  Q = K + noise_var*diag(1, n, n) + phi* diag(u, n, n)
  invQ = solve(Q)
  chol_decomp = chol(Q)
  f1 = 2 * sum(log(diag(chol_decomp)))
  f2 = (t(y - mu) %*% invQ %*% (y - mu))[1, 1]
  return(0.5*(f1 + f2))
}

gradient_phi = function(phi, X, y, mu, K, u, noise_var) {
  n = dim(X)[1]
  Un = diag(u, n, n)
  Q = K + noise_var*diag(1, n, n) + phi* Un
  invQ = solve(Q)
  alpha = invQ %*% (y - mu)
  df = - 0.5 * sum(diag((alpha %*% t(alpha) - invQ) %*% Un))
  return(df)
}

# noise_loss = function(noise_var, X, y, mu, K, u, phi) {
#   n = dim(X)[1]
#   Q = K + noise_var*diag(1, n, n) + phi * diag(u, n, n) 
#   chol_decomp = chol(Q)
#   f1 = 2 * sum(log(diag(chol_decomp)))
#   invQ = solve(Q)
#   f2 = (t(y - mu) %*% invQ %*% (y - mu))[1, 1]
#   return(0.5*(f1 + f2))
# }
# 
# gradient_noise = function(noise_var, X, y, mu, K, u, phi) {
#   n = dim(X)[1]
#   Q = K + noise_var*diag(1, n, n) + phi * diag(u, n, n) 
#   invQ = solve(Q)
#   alpha = invQ %*% (y - mu)
#   df = - 0.5 * sum(diag(alpha %*% t(alpha) - invQ))
#   return(df)
# }

# sqrt(mean((y_test - y_pred)**2))
