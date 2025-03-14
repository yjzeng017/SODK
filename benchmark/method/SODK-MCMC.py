
import gpflow
import tensorflow as tf
import tensorflow_probability as tfp
from gpflow.ci_utils import reduce_in_tests
import numpy as np
import random


f64 = gpflow.utilities.to_default_float
int32 = gpflow.utilities.to_default_int

def ODKMCMC_py(X, y, kernel, omega, jitter=1e-4, num_sample=5000, num_burnin = 2000):
  n = X.shape[0]
  d = X.shape[1]
  
  gpflow.config.set_default_jitter(jitter)
  mean_function = gpflow.mean_functions.Constant()
  random.seed(123)
  # variance_priors = []
  # for i in range(n):
  #   variance_priors = variance_priors + gpflow.Parameter(1.0, transform=gpflow.utilities.positive(),
  #                       prior=tfp.distributions.Gamma(f64(0.8), f64(0.8)))
  #
  likelihood = gpflow.likelihoods.Gaussian(variance=np.ones(shape=(n, 1)))
  model = gpflow.models.GPR(data=(X, y), kernel=kernel, mean_function=mean_function,likelihood=likelihood)
  model.kernel.lengthscales.prior = tfp.distributions.Gamma(f64(1), f64(1.0))
  model.kernel.variance.prior = tfp.distributions.Gamma(f64(1.0), f64(1.0))
  model.mean_function.c.prior = tfp.distributions.Normal(f64(0.0), f64(10.0))
  #
  model.phi = gpflow.Parameter(1, 
  transform=gpflow.utilities.positive(), 
  prior=tfp.distributions.Gamma(f64(1/2), f64(1/2)))
  
  model.likelihood.variance.prior = tfp.distributions.Gamma(f64(1/omega), f64(1/(model.phi*omega)))
  

  num_burnin_steps = reduce_in_tests(int32(num_burnin))
  num_samples = reduce_in_tests(int32(num_sample))
  
  # Note that here we need model.trainable_parameters, not trainable_variables - only parameters can have priors!
  hmc_helper = gpflow.optimizers.SamplingHelper(
      model.log_posterior_density, model.trainable_parameters
  )
  
  hmc = tfp.mcmc.HamiltonianMonteCarlo(
      target_log_prob_fn=hmc_helper.target_log_prob_fn,
      num_leapfrog_steps=10,
      step_size=0.01,
  )
  
  adaptive_hmc = tfp.mcmc.SimpleStepSizeAdaptation(
      hmc,
      num_adaptation_steps=10,
      target_accept_prob=f64(0.75),
      adaptation_rate=0.1,
  )
  
  
  @tf.function
  def run_chain_fn():
      return tfp.mcmc.sample_chain(
          num_results=num_samples,
          num_burnin_steps=num_burnin_steps,
          current_state=hmc_helper.current_state,
          kernel=adaptive_hmc,
          trace_fn=lambda _, pkr: pkr.inner_results.is_accepted,
      )
  
  
  samples, _ = run_chain_fn()
  
  # parameter_samples = hmc_helper.convert_to_constrained_values(samples)
  # 
  # param_to_name = {
  #   name: param
  #   for name, param in gpflow.utilities.parameter_dict(model).items()
  #   }
  return model

