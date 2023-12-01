# Telesca model
library(splines)
library(MASS)

# it returns Omega and P
omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

# n_knots_h -> number of internal knots for warping functions
# n_knots_m -> number of internal knots for common shape function

# n_obs is the number of observations for each patient
telesca_model = function(niter = 1000, nburn = 1000, m_c0, sigma_c0, a_c, b_c, m_a0, sigma_a0, a_a, b_a,
                         a_lambda, b_lambda, a_phi, b_phi, a_eps, b_eps, n_knots_h, n_knots_m,
                         y, n_obs){
  order = 4 # equal to 4 for cubic splines
  p = n_knots_m + order
  q = n_knots_h + order
  Omega = omega_P(p)
  P = omega_P(q)
  
  n_patients = dim(y)[2] # y has patients in the columns
  
  time = matrix(NA, nrow = dim(y)[1], ncol=n_patients)
  for (i in 1:n_patients){
    time[1:n_obs[i],i] = seq(0,1, length.out = n_obs[i])
  }
  
  # knots and B-spline matrix for common shape function m
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  B_m = bs(time, knots = knots_m, intercept=T)
  
  # knots and B-spline matrix for warping functions h_i
  knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
  B_h = bs(time, knots = knots_h, intercept=T)
  
  # Upsilon
  nu <- c(rep(0,order),knots_h,rep(1,order))
  Upsilon <- (nu[order] - nu[1])/(order-1)
  for(i in 1:(n_knots_h+order-1)){
    Upsilon[i+1] <- (nu[i+order] - nu[i+1])/(order-1) + Upsilon[i]
  }
  
  beta_0 = rep(0, p) # prior mean of beta
  
  # MCMC acceptance rates
  accepts = numeric(n_patients*(q-2))
  
  nrun = niter + nburn
  
  # matrices and vectors in which we save all the iterations of the MCMC
  a_save = matrix(NA, nrow=nrun, ncol = n_patients)
  a_save[1,] = rep(m_a0, n_patients)
  
  c_save = matrix(NA, nrow=nrun, ncol = n_patients)
  c_save[1,] = rep(m_c0, n_patients)
  
  
  c_0_save = numeric(nrun)
  c_0_save[1] = m_c0
  
  sigmac_save = numeric(nrun) # it is sigma^2
  sigmac_save[1] = sigma_c0
  
  a_0_save = numeric(nrun)
  a_0_save[1] = m_a0
  
  sigmaa_save = numeric(nrun) # it is sigma^2
  sigmaa_save[1] = sigma_a0
  
  beta_save = matrix(NA, nrow=nrun, ncol = p)
  beta_save[1, ] = beta_0
  
  lambda_save = numeric(nrun)
  lambda_save[1] = 1
  
  # tensor with all phi_i
  phi_save = array(rep(NA, nrun *q*n_patients), dim=c(nrun , q , n_patients))
  for (i in 1:n_patients){
    phi_save[1, ,i] = Upsilon
  }
  # phi[a,b,c] -> a=run, b=component of vector phi_i, c=patient i
  
  # save first value of phi_save
  
  sigma_phi_save = numeric(nrun)
  sigma_phi_save[1] =1
  
  sigma_eps_save = numeric(nrun)
  sigma_eps_save[1] =1
  
  # MCMC
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(iter in 2:nrun){
    setTxtProgressBar(bar, iter)
  }
  
  
  for(iter in 2:nrun){
    # update of phi_i
    for (i in 1:n_patients){
    
    h_i = b_splines(q, time[1:n_obs[i],i], phi_save[iter-1,,i], transpose=T)    
    
    # A.2.1 beta
    # vec(A) = as.numeric(t(A))
    X = matrix(NA, nrow=max(n_obs), ncol=p)
    # X = matrix(NA, nrow=n_patients, ncol=p) # ncol = dim beta
    for (pat in 1:n_patients){
      h_pat = b_splines(q, time[1:n_obs[pat],pat], phi_save[iter-1,,pat], transpose=T)  
      X[pat, ] = as.numeric(t(a_save[iter-1, pat]*b_splines(p,h_pat)))
    }
    C = matrix(NA, nrow=n_patients, ncol=max(n_obs))
    for(pat in 1:n_patients){
      C[pat,1:n_obs[pat]] = rep(c_save[iter-1, pat], n_obs[pat])
    }
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = inv_sigma_beta + (1/sigma_eps_save[iter-1])*t(X)%*%X
    V_beta = inv(inv_V_beta)
    m_beta = (1/sigma_eps_save[iter-1])*V_beta%*%(t(X)%*%(t(y)-C))
    
    #A.2.2 a_0, c_0
    sigma_star = 1/(1/sigma_a0 + n_patients/sigmaa_save[iter-1])
    a_star = sigma_star*(sum(a_save[iter-1, ])/sigmaa_save[iter-1] + m_a0/sigma_a0)
    a_0_save[iter] = rnorm(n=1, mean=a_star, sd=sqrt(sigma_star))
      
    sigma_star = 1/(1/sigma_c0 + n_patients/sigmac_save[iter-1])
    c_star = sigma_star*(sum(c_save[iter-1, ])/sigmac_save[iter-1] + m_c0/sigma_c0)
    c_0_save[iter] = rnorm(n=1, mean=c_star, sd=sqrt(sigma_star))
    
      
    # A.2.3
    # c_i, a_i
    sigma_ca = matrix(0, nrow = 2, ncol=2)
    sigma_ca[1,1] = sigmac_save[iter-1]
    sigma_ca[2,2] = sigmaa_save[iter-1]
    W = cbind(rep(1, n_obs[i]), b_splines(p, h_i, beta_save[iter,], transpose = F))
    inv_sigma_l = inv(sigma_ca) + (1/sigma_eps_save[iter-1])*t(W)%*%W
    sigma_l = inv(inv_sigma_l)
    mu_l = sigma_l%*%(inv(sigma_ca)%*%c(c_0_save[iter], a_0_save[iter]) + 
                        (1/sigma_eps_save[iter-1])*t(W)%*%y[1:n_obs[i],i])
    
    vec = mvrnorm(n=1, mu = mu_l, Sigma = sigma_l)
    c_save[iter, i] = vec[1]
    a_save[iter, i] = vec[2]
    
    # A.2.4 sigma eps
    a_star = a_eps + 0.5*sum(n_obs)
    m_tilde = matrix(NA, nrow=max(n_obs), ncol=n_patients)
    for (pat in 1:n_patients){
      h_pat = b_splines(q, time[1:n_obs[pat],pat], phi_save[iter-1,,pat], transpose=T)  
      m_tilde[,pat] = c_save[iter-1, pat]*rep(1, n_obs[pat]) + 
        a_save[iter-1, pat]*b_splines(p, h_pat, beta_save[iter-1,], transpose = F)
    }
    sum_b = 0
    for(pat in 1:n_patients){
      sum_b = sum_b + t(y[1:n_obs[pat],]-m_tilde[1:n_obs[pat],])%*%(y[1:n_obs[pat],]-m_tilde[1:n_obs[pat],])
    }
    b_star = b_eps + 0.5*sum_b
    sigma_eps_save[iter] = 1/rgamma(a_star, 1/b_star)
   
    # A.2.5
    # sigma_c
    a_star = a_c + 0.5*n_patients
    b_star = b_c + 0.5*sum((c_save[iter, ]-c_0_save[iter])^2)
    sigmac_save[iter] = 1/rgamma(a_star, 1/b_star)
    
    # sigma_a
    a_star = a_a + 0.5*n_patients
    b_star = b_a + 0.5*sum((a_save[iter, ]-a_0_save[iter])^2)
    sigmaa_save[iter] = 1/rgamma(a_star, 1/b_star)
    
    #lambda
    # k = p
    a_star = a_lambda + 0.5*p
    b_star = b_lambda + 0.5*t(beta_save[iter,])%*%Omega%*%beta_save[iter,]
    lambda_save[iter] = 1/rgamma(a_star, 1/b_star)
    
    #sigma_phi
    a_star = a_phi + 0.5*n_patients*q
    prod = 0
    for(pat in 1:n_patients){
      prod = prod + t(phi_save[iter, ,pat] - Upsilon)%*%P%*%(phi_save[iter, ,pat] - Upsilon)
    }
    b_star = b_phi + 0.5*prod
    sigma_phi_save[iter] = 1/rgamma(a_star, 1/b_star)
    
  }
  }
  for(i in 1:n_patients){
    tmp_phi <- phi_save[, ,i] # patient i
    
    current_llik <- dmnorm(y = y_list[[i]], mu =  H_list[[i]] %*% beta_save[it - 1,],
                           prec = 1 / sig2_save[it - 1] * diag_m, log = TRUE, unnorm = TRUE)
    current_lprior <- dmnorm(y = phi[[i]], mu = Upsilon,
                             prec = 1 / lam2_save[it - 1] * Q, log = TRUE, unnorm = TRUE)
    for(j in 2:(q-1)){
      tmp_phi[j] <- runif(1, min = max(tmp_phi[j] - tune, tmp_phi[j-1]),
                          max = min(tmp_phi[j] + tune, tmp_phi[j+1]))
      tmp_Hp <- bs(Hq %*% tmp_phi, knots = knot_loc_p, intercept = TRUE)
      #tmp_Hp <- cbs(Hq %*% tmp_phi, int_p)
      cand_llik <- dmnorm(y = y_list[[i]], mu =  tmp_Hp %*% beta_save[it - 1,],
                          prec = 1 / sig2_save[it - 1] * diag_m, log = TRUE, unnorm = TRUE)
      cand_lprior <- dmnorm(y = tmp_phi, mu = Upsilon,
                            prec = 1 / lam2_save[it - 1] * Q, log = TRUE, unnorm = TRUE)
      lratio <- cand_llik + cand_lprior - current_llik - current_lprior
      
      if(log(runif(1)) < lratio){
        current_llik <- cand_llik
        current_lprior <- cand_lprior
        accepts[(i-1)*(q-2)+(j-1)] <- accepts[(i-1)*(q-2)+(j-1)] + 1
      }
      else{
        tmp_phi[j] <- phi[[i]][j]
      }
    }
    phi[[i]] <- tmp_phi
  }
  
  



}

# Upsilon = c(1,1)
# phi_save[1,1,1] = 1
# phi_save[1,1,2] = 1
# sum(phi_save[1, ,] - Upsilon)
# phi_save[1,,]- Upsilon

b_splines = function(n_knots, time_vec, coeff_vec=NULL, transpose=F){
  knots = seq(0, 1, length.out = n_knots+2)[-c(1,n_knots+2)]
  B = bs(time_vec, knots = knots, intercept=T)
  if(!is.null(coeff_vec)){
  if(transpose){
    return(t(B)%*%coeff_vec)
  }
  else {return(B%*%coeff_vec)}}
  else{
    return(t(B))
  }
}
