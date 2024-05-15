# Gibbs sampling for Bayesian lasso

vblasso = function(X,y,priors, I = 100, thresh = 0.001, init){
  
  # priors should be a list containing s, a_sigma, b_sigma, s_0^2
  s = priors$s
  a_sigma = priors$a_sigma
  b_sigma = priors$b_sigma
  s2_0 = priors$s2_0
  
  # I is number of mcmc iterations
  I = I
  # The threshold value
  thresh = thresh
  
  #Data and sizes
  N = dim(X)[1]
  D = dim(X)[2]
  bX = cbind(rep(1,N),X)
  
  # Output
  elbo = c()
  
  # Initialize: init is a list containing the initialization of the variational parameters
  mu_w = matrix(init$mu_w, D+1,1)
  Sigma_w = init$Sigma_w
  a_sigma_vb = init$a_sigma
  b_sigma_vb = init$b_sigma
  mu_tau = init$mu_tau
  
  # Initialize difference between elbo and number of iterations
  diff = Inf
  i = 1
  
  # Iterate
  while((i<=I)&(diff>thresh)){
    
    # VB sigma
    a_sigma_vb = a_sigma + N/2 + (D+1)/2
    b_sigma_vb = c(b_sigma + 1/2*(sum(y^2) - 2*t(y)%*%bX%*%mu_w + 
                                  sum((Sigma_w + mu_w%*%t(mu_w))*(t(bX)%*%bX+ diag(c(1/s2_0,mu_tau)))))) 
    
    # VB tau
    mu_tau =1/sqrt(s^2*(diag(Sigma_w)[-1] + mu_w[-1]^2)*a_sigma_vb/b_sigma_vb)
    
    # VB w
    Sigma_w = b_sigma_vb/a_sigma_vb*solve(t(bX)%*%bX+  diag(c(1/s2_0,mu_tau)))
    mu_w = a_sigma_vb/b_sigma_vb*Sigma_w%*%t(bX)%*%y
    
    # Compute elbo
    elbo[i] = -1/2*a_sigma_vb/b_sigma_vb*(sum(y^2) - 2*t(y)%*%bX%*%mu_w + 
                                           sum((Sigma_w + mu_w%*%t(mu_w))*(t(bX)%*%bX+ diag(c(1/s2_0,mu_tau))))) +  
      1/2*determinant(Sigma_w, logarithm = TRUE)$modulus - D*log(s) -1/s^2*sum(1/mu_tau) - b_sigma/b_sigma_vb*a_sigma_vb - 
      a_sigma_vb*log(b_sigma_vb) 
    
    
    if(i %% 10==0){
      print(paste("Number of iterations completed:", i, " ELBO=", elbo[i]))
    }
    
    # Update difference and epoch
    if(i>1){diff = elbo[i]- elbo[i-1]}
    i = i + 1
    
  }
  return(list(mu_w=mu_w, Sigma_w=Sigma_w, a_sigma = a_sigma_vb, b_sigma = b_sigma_vb, mu_tau = mu_tau, elbo = elbo))
}
