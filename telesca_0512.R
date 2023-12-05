# Telesca model
library(splines)
library(MASS)
library(fda)
library(pracma)
y = CanadianWeather$dailyAv[,,1]


# it returns Omega and P
omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

# n_knots_h -> number of internal knots for warping functions
# n_knots_m -> number of internal knots for common shape function

# niter = 1000
# nburn = 1000
# m_c0 = 0
# sigma_c0 = 1
# a_c=1
# b_c=1
# m_a0=1
# sigma_a0 = 1
# a_a=1
# b_a=1
# a_lambda=1
# b_lambda=1
# a_phi=1
# b_phi=1
# a_eps=1
# b_eps=1 
# n_knots_h = 300
# n_knots_m = 10
# y = y
# n_obs=rep(365,35)

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
  knots_m = seq(0, 1, length.out = n_knots_m+2-order)[-c(1,n_knots_m+2-order)]
  B_m = bs(time, knots = knots_m, intercept=T)
  
  # knots and B-spline matrix for warping functions h_i
  knots_h = seq(0, 1, length.out = n_knots_h+2-order)[-c(1,n_knots_h+2-order)]
  B_h = bs(time, knots = knots_h, intercept=T) #B_h has nrows = length of time, ncols = length(knots_h)+order
  
  # Upsilon
  nu <- c(rep(0,order),seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)],rep(1,order))
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
  # for(iter in 2:nrun){
  #   setTxtProgressBar(bar, iter)
  # }
  # 
  #iter = 2
  #i = 1
  for(iter in 2:nrun){
    setTxtProgressBar(bar, iter)
    #A.2.2 a_0, c_0
    sigma_star = 1/(1/sigma_a0 + n_patients/sigmaa_save[iter-1])
    a_star = sigma_star*(sum(a_save[iter-1, ])/sigmaa_save[iter-1] + m_a0/sigma_a0)
    a_0_save[iter] = rnorm(n=1, mean=a_star, sd=sqrt(sigma_star))
    
    sigma_star = 1/(1/sigma_c0 + n_patients/sigmac_save[iter-1])
    c_star = sigma_star*(sum(c_save[iter-1, ])/sigmac_save[iter-1] + m_c0/sigma_c0)
    c_0_save[iter] = rnorm(n=1, mean=c_star, sd=sqrt(sigma_star))
    
    # A.2.1 beta
    # vec(A) = as.numeric(t(A))
    # X = matrix(NA, nrow=max(n_obs), ncol=p)
    # X = matrix(NA, nrow=n_patients, ncol=p) # ncol = dim beta
    # for (pat in 1:n_patients){
    #   h_pat = b_splines(q, time[1:n_obs[pat],pat], phi_save[iter-1,,pat])  
    #   X[pat, ] = as.numeric(t(a_save[iter-1, pat]*b_splines(p,h_pat)))
    # }
    
    X = matrix(1, nrow = max(n_obs), ncol=p)
    C = matrix(NA, nrow=n_patients, ncol=max(n_obs))
    for(pat in 1:n_patients){
      C[pat,1:n_obs[pat]] = rep(c_save[iter-1, pat], n_obs[pat])
    }
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = inv_sigma_beta + (1/sigma_eps_save[iter-1])*t(X)%*%X
    V_beta = inv(inv_V_beta)
   # print(V_beta)
    m_beta = rep(1, p)
    beta_save[iter, ] =  mvrnorm(n=1, mu = m_beta, Sigma = V_beta)
    # m_beta = NULL
    # for(pat in 1:n_patients){
    #   for(k in 1:p){
    #     
    #   }
    #   s = (1/sigma_eps_save[iter-1])
    #   m_beta = c(m_beta, s)
    #     (1/sigma_eps_save[iter-1])*V_beta%*%(t(X)%*%(y-t(C)))  
    # }
    # 
    #V-beta pxp
    #X n_obs x p
    # C n_obs x n_pat
    
    #pxn_pat
    
    # A.2.4 sigma eps
    a_star = a_eps + 0.5*sum(n_obs)
    m_tilde = matrix(NA, nrow=max(n_obs), ncol=n_patients)
    for (pat in 1:n_patients){
      h_pat = b_splines(q, time[1:n_obs[pat],pat], phi_save[iter-1,,pat])  
      m_tilde[1:n_obs[pat],pat] = c_save[iter-1, pat]*rep(1, n_obs[pat]) + 
        a_save[iter-1, pat]*b_splines(p, h_pat, beta_save[iter-1,])
    }
    sum_b = 0
    for(pat in 1:n_patients){
      sum_b = sum_b + t(y[1:n_obs[pat],pat]-m_tilde[1:n_obs[pat],pat])%*%(y[1:n_obs[pat],pat]-m_tilde[1:n_obs[pat],pat])
    }
    b_star = as.numeric(b_eps + 0.5*sum_b)
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, 1/b_star)
    
    # A.2.5
    # sigma_c
    a_star = a_c + 0.5*n_patients
    b_star = b_c + 0.5*sum((c_save[iter-1, ]-c_0_save[iter])^2)
    sigmac_save[iter] = 1/rgamma(n=1,a_star, 1/b_star)
    
    # sigma_a
    a_star = a_a + 0.5*n_patients
    b_star = b_a + 0.5*sum((a_save[iter-1, ]-a_0_save[iter])^2)
    sigmaa_save[iter] = 1/rgamma(n=1,a_star, 1/b_star)
    #lambda
    # k = p
    a_star = a_lambda + 0.5*p
    b_star = b_lambda + 0.5*t(beta_save[iter,])%*%Omega%*%beta_save[iter,]
    lambda_save[iter] = 1/rgamma(n=1,a_star, 1/b_star)
    
    #sigma_phi
    a_star = a_phi + 0.5*n_patients*q
    prod = 0
    for(pat in 1:n_patients){
      prod = prod + t(phi_save[iter-1, ,pat] - Upsilon)%*%P%*%(phi_save[iter-1, ,pat] - Upsilon)
    }
    b_star = b_phi + 0.5*prod
    sigma_phi_save[iter] = 1/rgamma(n=1, a_star, 1/b_star)
    
    for (i in 1:n_patients){
      # Update of phi_i
      # phi[a,b,c] -> a=run, b=component of vector phi_i, c=patient i
      # phi_save[iter, 1, i]=1
      # for(b in 2:q){
      #   phi_save[iter,b ,i] =  phi_save[iter, b-1, i] + 0.01
      # }
    h_i = b_splines(q, time[1:n_obs[i],i], phi_save[iter-1,,i])    
    
    phi_save[iter, 1, i] = 0
    phi_save[iter, q, i] = 1
    for (b in 2:(q-1)){
      # #B_h has nrows = length of time, ncols = length(knots_h)+order = q
      # warped_t = B_h[,b]*phi_save[iter-1, b,i]
      # prod1 = b_splines(n_knots_m, warped_t, beta_save[iter,])*a_save[iter-1, i]
      # add1 = y[1:n_obs[i], i] - c_save[iter-1,i]*t(rep(1, n_obs[i])) - prod1
      # # da capire righe/colonne
      # term1 = t(add1)%*%add1/sigma_eps_save[iter]
      # term2 = (phi_save[iter-1, b, i]-Upsilon[b])^2*P[b,b]/sigma_phi_save[iter]
      # exp(-0.5*(term1+term2))
      
      
      
      phi_old=phi_save[iter-1, b,i] #x
      Upsilon = Upsilon[b]
      beta = beta_save[iter,]
      sigma_eps = sigma_eps_save[iter]
      a_i = a_save[iter-1, i]
      c_i = c_save[iter-1,i]
      P = P[b,b]
      sigma_phi = sigma_phi_save[iter]
      #n_knots_m = n_knots_m
      B_h = B_h[,b]
      y_i = y[1:n_obs[i], i]
      n_obs = n_obs[i]
      pi_phi_old = pi_phi(phi_old, beta, sigma_eps, a_i, c_i, P, sigma_phi, Upsilon, n_knots_m, B_h, y_i, n_obs_i)
      phi_new=runif(n=1, min=phi_save[iter,b-1,i] , max=1+phi_save[iter,b-1,i])   
      pi_phi_new = pi_phi(phi_new, beta, sigma_eps, a_i, c_i, P, sigma_phi, Upsilon, n_knots_m, B_h, y_i, n_obs_i)
      alpha = min(1, (pi_phi_new*dunif(x=phi_old, min=phi_save[iter,b-1,i], max=1+phi_save[iter,b-1,i]))/
                    (pi_phi_old*dunif(x=phi_new, min=phi_save[iter,b-1,i], max=1+phi_save[iter,b-1,i])))
      u = runif(n=1, min=0, max=1)
      if(u<alpha){
        #accept
        phi_save[iter, b, i] = phi_new
      }
      else{ phi_save[iter, b, i] = phi_old}
    # b_h is the j-th row of B_h (j of phi_[iter, j, i])
    # phi is phi[iter-1, j, i]
    
    
    
    # A.2.3
    # c_i, a_i
    sigma_ca = matrix(0, nrow = 2, ncol=2)
    sigma_ca[1,1] = sigmac_save[iter-1]
    sigma_ca[2,2] = sigmaa_save[iter-1]
    W = cbind(rep(1, n_obs[i]), b_splines(p, h_i, beta_save[iter,]))
    inv_sigma_l = inv(sigma_ca) + (1/sigma_eps_save[iter-1])*t(W)%*%W
    sigma_l = inv(inv_sigma_l)
    mu_l = sigma_l%*%(inv(sigma_ca)%*%c(c_0_save[iter], a_0_save[iter]) + 
                        (1/sigma_eps_save[iter-1])*t(W)%*%y[1:n_obs[i],i])
    
    vec = mvrnorm(n=1, mu = mu_l, Sigma = sigma_l)
    c_save[iter, i] = vec[1]
    a_save[iter, i] = vec[2]
    }
print(iter)  }
  
  
return (list (a = a_save, c=c_save, c_0 = c_0_save, sigmac = sigmac_save, a_0 = a_0_save, 
              sigmaa = sigmaa_save, beta = beta_save, lambda = lambda_save, phi = phi_save, 
              sigma_phi = sigma_phi_save, sigma_eps = sigma_eps_save))
}

# Upsilon = c(1,1)
# phi_save[1,1,1] = 1
# phi_save[1,1,2] = 1
# sum(phi_save[1, ,] - Upsilon)
# phi_save[1,,]- Upsilon

b_splines = function(n_knots, time_vec, coeff_vec=NULL){
  n_knots = n_knots-4
  knots = seq(0, 1, length.out = n_knots+2)[-c(1,n_knots+2)]
  B = bs(time_vec, knots = knots, intercept=T)
  if(!is.null(coeff_vec)){
  return(B%*%coeff_vec)
}
  else{
    return(B)
  }
}
# warped_t = B_h[,b]*phi_save[iter-1, b,i]
# B = bs(seq(0,1, len=100), knots = c(0,0.25,0.5,0.75,1), intercept=T)
# dim(B)

#y_i are data of patient i (vector)
# b_h is the j-th row of B_h (j of phi_[iter, j, i])
# phi is phi[iter-1, j, i]
# inputs

pi_phi = function (phi, beta, sigma_eps, a_i, c_i, P, sigma_phi, Upsilon, n_knots_m, B_h, y_i, n_obs_i){
  #B_h has nrows = length of time, ncols = length(knots_h)+order = q
  warped_t = B_h*phi
  prod1 = b_splines(n_knots_m, warped_t, beta)*a_i
  add1 = y_i - c_i*rep(1, n_obs) - prod1
  # da capire righe/colonne
  term1 = t(add1)%*%add1/sigma_eps
  term2 = (phi-Upsilon)^2*P/sigma_phi
  return(exp(-0.5*(term1+term2)))
}

res=telesca_model(niter = 1000, nburn = 1000, m_c0 = 0, sigma_c0 = 1, a_c=1, b_c=1, m_a0=1, 
                         sigma_a0 = 1, a_a=1, b_a=1,
                         a_lambda=1, b_lambda=1, a_phi=1, b_phi=1, a_eps=1, b_eps=1, 
                         n_knots_h = 300, n_knots_m =4 ,
                         y = y, n_obs=c(rep(365,34),306)))
