# Gibbs sampling for Bayesian lasso

gibbs_blasso = function(X,y,priors, I, w_init, sigma_init){
  
  # priors should be a list containing s, a_sigma, b_sigma, s_0^2
  s = priors$s
  a_sigma = priors$a_sigma
  b_sigma = priors$b_sigma
  s2_0 = priors$s2_0
  
  # I is number of mcmc iterations
  I = I
  
  #Data and sizes
  N = dim(X)[1]
  D = dim(X)[2]
  bX = cbind(rep(1,N),X)
  
  # Output
  w_mat = matrix(0, I, D+1)
  sigma_mat = matrix(0, I, 1)
  tau_mat = matrix(0, I, D)
  
  # Initialize
  w = w_init
  sigma = sigma_init
  tau = matrix(0,D,1)
  
  # Iterate over mcmc
  for (i in 1:I){
    
    # Sample tau
    tau = 1/rinvgauss(D, mean=sqrt(sigma/(s^2*w^2)), shape=1/s^2)
    
    # Sample w
    S = solve(t(bX)%*%bX+ diag(1/c(s2_0,tau)))
    mu = S%*%t(bX)%*%y
    w = t(rmvnorm(1, mean = mu, sigma = S))
    
    # Sample sigma
    ahat = a_sigma + N/2 + (D+1)/2
    bhat = b_sigma + 1/2*sum((y-bX%*%w)^2) + 1/2* sum(w^2/c(s2_0,tau))
    sigma = 1/rgamma(1, shape = ahat, rate = bhat)
    
    #Store output
    w_mat[i,] = w
    sigma_mat[i] = sigma
    tau_mat[i,] = tau
    
    if(i %% 100==0){
      print(paste("Number of iterations completed:", i))
    }
    
  }
  return(list(w=w_mat, sigma=sigma_mat, tau=tau_mat))
}