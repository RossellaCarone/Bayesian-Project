
# Telesca model

library(splines)
library(MASS)
library(fda)
library(pracma)
library(ddalpha)
library(emdbook)
library(RcppArmadillo)
library(splines2)
library(pbapply)
library(parallel)
cores = detectCores()
cores
#cl=makeCluster(parallel::detectCores()/2)
#clusterExport(cl=cl,list('x1'))

rm(list=ls())
graphics.off()

Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends="RcppArmadillo")
Rcpp::cppFunction("arma::vec Arma_mvrnorm(const arma::vec& mean, const arma::mat& cov) {
  return arma::mvnrnd(mean, cov);}", depends="RcppArmadillo")
Rcpp::cppFunction("double  Arma_rnorm(const double &mu, const double &sigma){
      return arma::randn(arma::distr_param(mu,sigma));} ", depends="RcppArmadillo") #variance = sigma^2
Rcpp::cppFunction("double  Arma_runif(const double &a, const double &b){
      return arma::randu(arma::distr_param(a,b));} ", depends="RcppArmadillo")


# it returns Omega and P
omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  return (K)
}

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

b_splines_new = function(n_knots, time_vec, coeff_vec=NULL){
  n_knots = n_knots-4
  knots = seq(0, 1, length.out = n_knots+2)[-c(1,n_knots+2)]
  B = bSpline(time_vec, knots = knots, intercept=T)
  if(!is.null(coeff_vec)){
    return(B%*%coeff_vec)
  }
  else{
    return(B)
  }
}


pi_phi = function (phi, beta, sigma_eps, a_i, c_i, P, sigma_phi, Upsilon, p, B_h, y_i, n_obs_i){
  warped_t = B_h*phi
  prod1 = b_splines(p, warped_t, beta)*a_i
  add1 = y_i - rep(c_i, n_obs_i) - prod1
  term1 = crossprod(add1)/sigma_eps
  term2 = (phi-Upsilon)^2*P/sigma_phi
  return(-0.5*(term1+term2))
}


# Inputs:
# da completare
# n_knots_h -> number of internal knots for warping functions
# n_knots_m -> number of internal knots for common shape function


# Outputs:
# da completare

#0.02 solo gibbs
telesca_model = function(niter = 1000, nburn = 1000, m_c0, sigma_c0, a_c, b_c, m_a0, sigma_a0, a_a, b_a,
                         a_lambda, b_lambda, a_phi, b_phi, a_eps, b_eps, n_knots_h, n_knots_m,
                         y){
  
  n_obs_max  = dim(y)[1] 
  n_patients = dim(y)[2] # y has patients in the columns
  y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))
  n_obs=sapply(y_list,length)
  y=unlist(y_list)
  y=as.numeric(y)
  order = 4 # equal to 4 for cubic splines
  p = n_knots_m + order
  q = n_knots_h + order
  Omega = omega_P(p)
  inv_Omega = armaInv(Omega)
  P = omega_P(q)
  inv_P = armaInv(P)
  
  time = matrix(NA, nrow = n_obs_max, ncol=n_patients)
  for (i in 1:n_patients){
    time[1:n_obs[i],i] = seq(0,1, length.out = n_obs[i])
  }
  
  # knots and B-spline matrix for common shape function m
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  #B_m = bs(time, knots = knots_m, intercept=T)
  
  # knots and B-spline matrix for warping functions h_i
  knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
  #B_h = bs(time, knots = knots_h, intercept=T) #B_h has nrows = length of time, ncols = length(knots_h)+order
  
  # Upsilon
  nu <- c(rep(0,order),seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)],rep(1,order))
  Upsilon <- (nu[order] - nu[1])/(order-1)
  for(i in 1:(n_knots_h+order-1)){
    Upsilon[i+1] <- (nu[i+order] - nu[i+1])/(order-1) + Upsilon[i]
  }
  
  beta_0 = rep(0, p) # prior mean of beta
  
  # MCMC acceptance rates
  nrun = niter + nburn
  accepts = matrix(0, nrow=n_patients, ncol=q-2)
  
  # matrices and vectors in which we save all the iterations of the MCMC
  a_save = matrix(NA, nrow=nrun, ncol = n_patients)
  a_save[1,] = rnorm(n=1, mean=m_a0, sd=sqrt(sigma_a0))
  
  c_save = matrix(NA, nrow=nrun, ncol = n_patients)
  c_save[1,] = rnorm(n=1, mean=m_c0, sd=sqrt(sigma_c0))
  
  
  c_0_save = numeric(nrun)
  c_0_save[1] = rnorm(n=1, mean=m_c0, sd=sqrt(sigma_c0))
  
  sigmac_save = numeric(nrun) # it is sigma^2
  sigmac_save[1] = 1/rgamma(n=1,a_c, b_c)
  
  a_0_save = numeric(nrun)
  a_0_save[1] = rnorm(n=1, mean=m_a0, sd=sqrt(sigma_a0))
  
  sigmaa_save = numeric(nrun) # it is sigma^2
  sigmaa_save[1] = 1/rgamma(n=1,a_a, b_a)
  
  lambda_save = numeric(nrun)
  lambda_save[1] = 1/rgamma(n=1,a_lambda, b_lambda)
  
  beta_save = matrix(NA, nrow=nrun, ncol = p)
  beta_save[1, ] = Arma_mvrnorm(beta_0, lambda_save[1] * inv_Omega)
  
  h_save = array(rep(NA, nrun *n_obs_max*n_patients), dim=c(nrun , n_obs_max , n_patients))
  
  # tensor with all phi_i
  phi_save = array(rep(NA, nrun *q*n_patients), dim=c(nrun , q , n_patients))
  for (i in 1:n_patients){
    phi_save[1, ,i] = Upsilon
  }
  # phi[a,b,c] -> a=run, b=component of vector phi_i, c=patient
  
  phi <- lapply(1:n_patients, function(i) Upsilon)  # serve a salvare solo l'ultima
  
  sigma_phi_save = numeric(nrun)
  sigma_phi_save[1] =1
  
  sigma_eps_save = numeric(nrun)
  sigma_eps_save[1] =1
  
  
  phi_save = array(rep(NA, nrun *q*n_patients), dim=c(nrun , q , n_patients))
  for (i in 1:n_patients){
    phi_save[1, ,i] = Upsilon
  }
  # phi[a,b,c] -> a=run, b=component of vector phi_i, c=patient
  
  B_h = list()
  for(i in 1:n_patients){
    B_h[[i]]=bSpline(time[1:n_obs[i],i], knots = knots_h, intercept=T)
  }
  # MCMC
  
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
 
  for(iter in 2:nrun){
    
    setTxtProgressBar(bar, iter)
    
    # Update of phi_i (IT DOESN'T WORK!!)
    for (i in 1:n_patients){
      
      h_i = b_splines(q, time[1:n_obs[i],i], phi_save[iter-1,,i])    
      h_save[iter-1,1:n_obs[i],i] = h_i[1:n_obs[i]]
      
      phi_save[iter, 1, i] = 0
      phi_save[iter, q, i] = 1
      
      
      for (b in 2:(q-1)){
        
        pi_phi_old = pi_phi(phi_save[iter-1, b,i] , beta_save[iter-1,],  sigma_eps_save[iter-1], a_save[iter-1, i],
                            c_save[iter-1,i], P[b,b], sigma_phi_save[iter-1], Upsilon[b], p,  B_h[[i]][,b],
                            as.numeric(y_list[[i]]), n_obs[i])
  
        
        phi_old = phi_save[iter-1, b,i]
        phi_new=Arma_runif(phi_save[iter,b-1,i], phi_save[iter-1,b+1,i])  
        
        pi_phi_new = pi_phi(phi_new, beta_save[iter-1,], sigma_eps_save[iter-1],a_save[iter-1, i], 
                            c_save[iter-1,i], P[b,b], sigma_phi_save[iter-1], Upsilon[b],p,  B_h[[i]][,b],
                            as.numeric(y_list[[i]]), n_obs[i])
       
        alpha = min(0, (pi_phi_new-pi_phi_old))
        
        u = Arma_runif(0, 1)
        
        if(u<exp(alpha)){
          #accept
          phi_save[iter, b, i] = phi_new
          accepts[i, b-1]= accepts[i, b-1]+1
        } else{ 
          phi_save[iter, b, i] = phi_old
        }
        
      }
    }
    
    
    
    
    # A.2.1 beta
   
    h_pat_list = NULL
    X = NULL
    C = NULL
    #h_pat_list = lapply(1:n_patients, function(pat) b_splines(q, time[1:n_obs[pat],pat], phi_save[iter,,pat]))
    #start.time <- Sys.time()
    for (pat in 1:n_patients){
       h_pat_list[[pat]] = b_splines(q, time[1:n_obs[pat],pat], phi_save[iter,,pat])
       X= rbind(X,a_save[iter-1, pat]*b_splines(p,h_pat_list[[pat]]))
       C = c(C,rep(c_save[iter-1, pat], n_obs[pat]))
     }
    #end.time <- Sys.time()
    #time.taken <- round(end.time - start.time,5)
    #print(paste('time for: ',time.taken))
    
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = inv_sigma_beta + (1/sigma_eps_save[iter-1])*crossprod(X)
    
    V_beta = armaInv(inv_V_beta)
    m_beta =(1/sigma_eps_save[iter-1])*V_beta%*%(crossprod(X, y-C))
    beta_save[iter, ] =  Arma_mvrnorm(m_beta, V_beta)
    
    
    #A.2.2 
    
    #a_0, c_0
    
    sigma_star = 1/(1/sigma_a0 + n_patients/sigmaa_save[iter-1])
    a_star = sigma_star*(sum(a_save[iter-1, ])/sigmaa_save[iter-1] + m_a0/sigma_a0)
    a_0_save[iter] = Arma_rnorm(a_star, sqrt(sigma_star))
    
    sigma_star = 1/(1/sigma_c0 + n_patients/sigmac_save[iter-1])
    c_star = sigma_star*(sum(c_save[iter-1, ])/sigmac_save[iter-1] + m_c0/sigma_c0)
    c_0_save[iter] = Arma_rnorm(c_star, sqrt(sigma_star))
    
    # A.2.3
    
    # c_i, a_i
    
    for (i in 1:n_patients){
      
      sigma_ca = matrix(0, nrow = 2, ncol=2)
      sigma_ca[1,1] = sigmac_save[iter-1]
      sigma_ca[2,2] = sigmaa_save[iter-1]
      W = cbind(rep(1, n_obs[i]), b_splines(p, h_pat_list[[i]], beta_save[iter,]))
      inv_sigma_l = armaInv(sigma_ca) + (1/sigma_eps_save[iter-1])*crossprod(W)
      sigma_l = armaInv(inv_sigma_l)
      mu_l = sigma_l%*%(armaInv(sigma_ca)%*%c(c_0_save[iter], a_0_save[iter]) + 
                          (1/sigma_eps_save[iter-1])*crossprod(W,y_list[[i]]))
      vec = Arma_mvrnorm(mu_l, sigma_l)
      c_save[iter, i] = vec[1]
      a_save[iter, i] = vec[2]
      
    }

    
    # A.2.4 
    
    #sigma eps
    a_star = a_eps + 0.5*sum(n_obs)
    m_tilde = matrix(NA, nrow=max(n_obs), ncol=n_patients)
    for (pat in 1:n_patients){
      m_tilde[1:n_obs[pat],pat] = c_save[iter, pat]*rep(1, n_obs[pat]) + 
        a_save[iter, pat]*b_splines(p, h_pat_list[[pat]], beta_save[iter,])
    }
    sum_b = 0
    for(pat in 1:n_patients){
      sum_b = sum_b + crossprod(y_list[[pat]]-m_tilde[1:n_obs[pat],pat])
    }
    b_star = as.numeric(b_eps + 0.5*sum_b)
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    # A.2.5
    
    # sigma_c
    a_star = a_c + 0.5*n_patients
    b_star = b_c + 0.5*sum((c_save[iter, ]-c_0_save[iter])^2)
    sigmac_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    # sigma_a
    a_star = a_a + 0.5*n_patients
    b_star = b_a + 0.5*sum((a_save[iter, ]-a_0_save[iter])^2)
    sigmaa_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    #lambda
    a_star = a_lambda + 0.5*p
    b_star = b_lambda + 0.5*t(beta_save[iter,])%*%Omega%*%beta_save[iter,]
    lambda_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    #sigma_phi
    a_star = a_phi + 0.5*n_patients*q
    prod = 0
    for(pat in 1:n_patients){
      prod = prod + t(phi_save[iter, ,pat] - Upsilon)%*%P%*%(phi_save[iter, ,pat] - Upsilon)
    }
    b_star = b_phi + 0.5*prod
    sigma_phi_save[iter] = 1/rgamma(n=1, a_star, b_star)
    
    
  }

  
    close(bar)
    accepts <- accepts / nrun
    
    a_post <- apply(a_save[-c(1:nburn),], 2, mean)
    c_post <- apply(c_save[-c(1:nburn),], 2, mean)
    c_0_post <- mean(c_0_save[-c(1:nburn)])
    sigmac_post<-mean(sigmac_save[-c(1:nburn)])
    a_0_post <- mean(a_0_save[-c(1:nburn)])
    sigmaa_post<-mean(sigmaa_save[-c(1:nburn)])
    beta_post <- apply(beta_save[-c(1:nburn),], 2, mean)
    lambda_post <- mean(lambda_save[-c(1:nburn)])
    phi_post<-apply(phi_save[-c(1:nburn),,], c(2, 3), mean)
    sigma_phi_post<-mean(sigma_phi_save[-c(1:nburn)])
    sigma_eps_post<-mean(sigma_eps_save[-c(1:nburn)])

    
    for(i in 1:n_patients){
      h_i = b_splines(q, time[1:n_obs[i],i], phi_save[iter,,i])    
      h_save[iter,1:n_obs[i],i] = h_i[1:n_obs[i]]
    }
    h_post = apply(h_save[-c(1:nburn),,], c(2, 3), mean)
    y_star= matrix(NA, nrow=max(n_obs), ncol=n_patients)
    
    for(i in 1:n_patients){
      warped_time = seq(min(h_post[1:n_obs[i],i]), max(h_post[1:n_obs[i],i]) , length.out=n_obs[i])
      y_star[1:n_obs[i],i] =     splinefun(h_post[1:n_obs[i],i], y_list[[i]]) (warped_time)
    }
  
   
return (list (post=list(a = a_post, c=c_post, c_0 = c_0_post, sigmac = sigmac_post, a_0 = a_0_post, 
              sigmaa = sigmaa_post, beta = beta_post, lambda = lambda_post, phi = phi_post, 
              sigma_phi = sigma_phi_post, sigma_eps = sigma_eps_post, accepts = accepts ), 
              full=list(a_s=a_save, c_s = c_save, c0_s=c_0_save, sigmac_s=sigmac_save, 
                   a0_s=a_0_save, sigmaa_s=sigmaa_save, beta_s=beta_save,
                   lambda_s=lambda_save, phi_s=phi_save, sigma_phi_s=sigma_phi_save,
                   sigma_eps_s=sigma_eps_save),
              y_star=y_star))
}


#Now we test the model
m_c0=50
sigma_c0=10
a_c=10
b_c=1
m_a0=50
sigma_a0=10
a_a=1
b_a=1
a_lambda=10
b_lambda=1
a_phi=10
b_phi=1
a_eps=1
b_eps=1
n_knots_h=7
n_knots_m=27
y = growth$hgtm[-c(1,2,3,4),]  
#y = CanadianWeather$dailyAv[,,1]

nburn = 1000
niter = 1000
telesca_output <-telesca_model(niter = niter, nburn = nburn, m_c0, sigma_c0, a_c, b_c, m_a0, sigma_a0, a_a, b_a,
                                         a_lambda, b_lambda, a_phi, b_phi, a_eps, b_eps, n_knots_h, n_knots_m,
                                         y)

#0.00017
#0.00060

acceptance_output = telesca_output$post$accepts
x11()
hist(acceptance_output, main='Acceptance rate')

a_output=telesca_output$full$a_s
c_output=telesca_output$full$c_s

c0_output = telesca_output$full$c0_s
sigma_c_output = telesca_output$full$sigmac_s

a0_output = telesca_output$full$a0_s
sigma_a_output = telesca_output$full$sigmaa_s

beta_output=telesca_output$full$beta_s

lambda_output = telesca_output$full$lambda_s

phi_output = telesca_output$full$phi_s
sigma_phi_output = telesca_output$full$sigma_phi_s

sigma_eps_output = telesca_output$full$sigma_eps_s

x11()
par(mfrow=c(2,1))
matplot(a_output, type='l', main='a_i', xlab='Iteration', xlim=c(nburn+1, niter+nburn))
matplot(c_output, type='l', main='c_i', xlab='Iteration', xlim=c(nburn+1, niter+nburn))

x11()
par(mfrow=c(2,1))
matplot(a0_output, type='l', main='a0', xlab='Iteration', xlim=c(nburn+1, niter+nburn))
matplot(c0_output, type='l', main='c0', xlab='Iteration', xlim=c(nburn+1, niter+nburn))

x11()
par(mfrow=c(2,1))
matplot(sigma_a_output, type='l', main='sigma^2_a0', xlab='Iteration', xlim=c(nburn+1, niter+nburn))
matplot(sigma_c_output, type='l', main='sigma^2_c0', xlab='Iteration', xlim=c(nburn+1, niter+nburn))


x11()
matplot(beta_output, xlim=c(nburn+1, niter+nburn),type='l', ylim=range(beta_output[(nburn+1):(niter+nburn),]),
        xlab = 'Iteration', main='beta')

x11()
matplot(lambda_output, xlim=c(nburn+1, niter+nburn),type='l', ylim=range(lambda_output[(nburn+1):(niter+nburn)]),
        xlab = 'Iteration', main='lambda')

x11()
par(mfrow=c(2,1))
matplot(phi_output[,,1], xlim=c(nburn+1, niter+nburn),type='l',xlab = 'Iteration', main='phi of patient 1')
matplot(phi_output[,,5], xlim=c(nburn+1, niter+nburn),type='l',xlab = 'Iteration', main='phi of patient 5')

x11()
par(mfrow=c(2,1))
matplot(sigma_phi_output, xlim=c(nburn+1, niter+nburn),type='l',xlab = 'Iteration', 
        main='sigma^2_phi', ylim=range(sigma_phi_output[(nburn+1):(niter+nburn)]))
matplot(sigma_eps_output, xlim=c(nburn+1, niter+nburn),type='l',xlab = 'Iteration', 
        main='sigma^2_eps')

#graphics.off()
n_patients=39

n_times = 27
time = matrix(NA, nrow = n_times, ncol=n_patients)
for (i in 1:n_patients){
  time[1:n_times,i] = seq(0,1, length.out = n_times)
}

time_grid = seq(0,1, length.out=1000)

knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
#B_h = bSpline(1:n_times, knots = knots_h, intercept=T)

h_i = matrix(NA, nrow=length(time_grid), ncol=n_patients )
for (i in 1:n_patients){
  #h_i[,i] = b_splines(n_knots_h+4, 1:27, coeff_vec=spero_funzioni$post$phi[,i])
  h_i[,i] = b_splines(n_knots_h+4, time_grid, coeff_vec=telesca_output$post$phi[,i])
}

#m = b_splines(n_knots_m+4, 1:27, coeff_vec=spero_funzioni$post$beta)
m = b_splines(n_knots_m+4, time_grid, coeff_vec=telesca_output$post$beta)

m_i = matrix(NA, nrow=length(time_grid), ncol=n_patients )
for (i in 1:n_patients){
  m_i[,i] = telesca_output$post$c[i] + telesca_output$post$a[i]*m
}

m_tilde = matrix(NA, nrow=length(time_grid), ncol=n_patients )
for (i in 1:n_patients){
  m_tilde[,i] = telesca_output$post$c[i] + telesca_output$post$a[i]*b_splines(n_knots_m+4, h_i[,i], coeff_vec=telesca_output$post$beta)
}

x11()
par(mfrow=c(1,2))
matplot(y, type='l', main='data')
matplot(m_tilde, type='l', main='m_tilde')
x11()
par(mfrow=c(1,2))
matplot(m_i, type='l', main='m_i')
matplot(m, type='l', main='m')
x11()
par(mfrow=c(1,2))
matplot(y, type='l', main='y')
matplot(telesca_output$y_star, type='l', main='y*')


#graphics.off()

#### data knee ####

#Now we test the model with the final data
#inv_x = armaInv(x)
# P = omega_P(5)
# inv_P = solve(P)
# inv2_P = armaInv(P)

m_c0=30
sigma_c0=10
a_c=1e-4
b_c=1
m_a0=30
sigma_a0=10
a_a=1e-4
b_a=1
a_lambda=1e-2
b_lambda=1
a_phi=1e-4
b_phi=1
a_eps=1e-4
b_eps=1
n_knots_h=7
n_knots_m=600 #600
y = read.csv('data.csv')
y = y[,1:7] # i consider only 7 patients

spero_funzioni<-telesca_model(niter = 200, nburn = 300, m_c0, sigma_c0, a_c, b_c, m_a0, sigma_a0, a_a, b_a,
                              a_lambda, b_lambda, a_phi, b_phi, a_eps, b_eps, n_knots_h, n_knots_m,
                              y)
# circa 2-3 sec a iterazione
# 23:18

beta_full=spero_funzioni$full$beta
#acceptance = spero_funzioni$post$accepts
a=spero_funzioni$full$a_s
c=spero_funzioni$full$c_s

x11()
par(mfrow=c(2,1))
matplot(a, type='l')
matplot(c, type='l')
x11()
matplot(beta_full, type='l')

n_obs_max  = dim(y)[1] 
n_patients = dim(y)[2] # y has patients in the columns
y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))
n_obs=sapply(y_list,length)

time = matrix(NA, nrow = n_obs_max, ncol=n_patients)
for (i in 1:n_patients){
  time[1:n_obs[i],i] = seq(0,1, length.out = n_obs[i])
}


knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
#B_h = bs(1:27, knots = knots_h, intercept=T)

h_i = matrix(NA, nrow=n_obs_max, ncol=n_patients )
for (i in 1:n_patients){
  h_i[1:n_obs[i],i] = b_splines(n_knots_h+4, time[1:n_obs[i],i], coeff_vec=spero_funzioni$post$phi[,i])
}

#m = b_splines(n_knots_m+4, time[1:n_obs[i],i], coeff_vec=spero_funzioni$post$beta)

m_i = matrix(NA, nrow=n_obs_max, ncol=n_patients )
for (i in 1:n_patients){
  m_i[1:n_obs[i],i] = spero_funzioni$post$c[i] + spero_funzioni$post$a[i]*b_splines(n_knots_m+4, time[1:n_obs[i],i], coeff_vec=spero_funzioni$post$beta) 
}

m_tilde = matrix(NA, nrow=n_obs_max, ncol=n_patients )
for (i in 1:n_patients){
  m_tilde[1:n_obs[i],i] = spero_funzioni$post$c[i] + spero_funzioni$post$a[i]*b_splines(n_knots_m+4, h_i[1:n_obs[i],i], coeff_vec=spero_funzioni$post$beta)
}


x11()
par(mfrow=c(1,2))
matplot(y, type='l')
matplot(m_tilde, type='l')
x11()
#par(mfrow=c(1,2))
matplot(m_i, type='l')
#matplot(m, type='l')
phi_full=spero_funzioni$full$phi
phi_1=phi_full[,,1]
x11()
matplot(phi_1, type='l')


