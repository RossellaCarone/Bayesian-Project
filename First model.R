rm(list=ls())
#### Libraries ####
library(splines)
library(MASS)
library(fda)
library(ddalpha)
library(emdbook)
library(RcppArmadillo)
library(splines2)
library(mvnfast)
library(coda)
library(LearnBayes)
library(lattice)


#### Rcpp functions ####
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends="RcppArmadillo")
Rcpp::cppFunction("arma::vec Arma_mvrnorm(const arma::vec& mean, const arma::mat& cov) {
  return arma::mvnrnd(mean, cov);}", depends="RcppArmadillo")
Rcpp::cppFunction("double  Arma_rnorm(const double &mu, const double &sigma){
      return arma::randn(arma::distr_param(mu,sigma));} ", depends="RcppArmadillo") #variance = sigma^2
Rcpp::cppFunction("double  Arma_runif(const double &a, const double &b){
      return arma::randu(arma::distr_param(a,b));} ", depends="RcppArmadillo")


#### Telesca model function + auxiliary functions ####
pi_phi = function (phi, beta, sigma_eps, a_i, c_i, P, sigma_phi, Upsilon, p, B_h, y_i, n_obs_i){
  warped_t = B_h%*%phi
  prod1 = b_splines(p, warped_t, beta)*a_i
  add1 = y_i - rep(c_i, n_obs_i) - prod1
  term1 = crossprod(add1)/sigma_eps
  term2 = t(phi-Upsilon)%*%P%*%(phi-Upsilon)/sigma_phi
  return(-0.5*(term1+term2))
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

# it returns Omega and P
omega_P <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diag=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(rep(2,dim-1),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  return (K)
}

stack_Matrix<-function(x_list){
  n=length(x_list)
  X=NULL
  for (i in 1:n){
    X= rbind(X,x_list[[i]])
  }
  return (X)
}


interp_spline <- function(x, y, nout = length(y)) {
  ind_out <- seq(min(x), max(x), len = nout)
  spfit <- splinefun(x, y)
  return(spfit(ind_out))
}

telesca_model <- function(y, niter = 100, nburn = 100, n_knots_h = 10, n_knots_m = 30,m_c0=0, sigma_c0=1, a_c=0, b_c=1, m_a0=1, sigma_a0=1, a_a=1, b_a=1,
                          a_eps = .1, b_eps = .1, a_lambda = .1, b_lambda = .1, a_phi = .1, b_phi = .1
){
  
  
  n_obs_max  = dim(y)[1] 
  n_patients = dim(y)[2] # y has patients in the columns
  y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))  # na.omit?
  n_obs=sapply(y_list,length)
  y=unlist(y_list)
  y=as.numeric(y)
  order = 4 # equal to 4 for cubic splines
  p = n_knots_m + order
  q = n_knots_h + order
  Omega = omega_P(p)
  P = omega_P(q)
  inv_P = armaInv(P)
  time= lapply(1:n_patients, function(i)  seq(0,1, length.out = n_obs[i]))
  
  # knots and Bspline matrix for warping functions h_i
  knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]
  Bh <- lapply(1:n_patients, function(i) bs(time[[i]], knots = knots_h, intercept = T))
  
  # knots and Bspline matrix for common shape function m
  knots_m = seq(0, 1, length.out = n_knots_m+2)[-c(1,n_knots_m+2)]
  Bm <- lapply(1:n_patients, function(i) bs(time[[i]], knots = knots_m, intercept = T))
  
  # Upsilon
  nu <- c(rep(0,order),seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)],rep(1,order))
  Upsilon <- (nu[order] - nu[1])/(order-1)
  for(i in 1:(n_knots_h+order-1)){
    Upsilon[i+1] <- (nu[i+order] - nu[i+1])/(order-1) + Upsilon[i]
  }
  
  bigP <- Matrix::bdiag(replicate(n_patients, P, simplify = FALSE))
  beta_0 <- rep(0,p)
  
  tune <- .005
  accepts = matrix(0, nrow=n_patients, ncol=q-2)
  
  
  #----- Save Structures -----#
  
  nrun <- nburn + niter
  
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
  
  beta_save <- matrix(NA, nrow = nrun, ncol = p)
  beta_save[1,] <- beta_0
  
  sigma_eps_save <- numeric(nrun)
  sigma_eps_save[1] <- 1
  lambda_save <- numeric(nrun)
  lambda_save[1] <- 1
  sigma_phi_save <- numeric(nrun)
  sigma_phi_save[1] <- 1
  
  phi <- lapply(1:n_patients, function(i) Upsilon)
  
  # tensor with all phi_i
  phi_save = array(rep(NA, nrun *q*n_patients), dim=c(nrun , q , n_patients))
  for (i in 1:n_patients){
    phi_save[1, ,i] = Upsilon
  }
  # phi[a,b,c] -> a=run, b=component of vector phi_i, c=patient
  
  h <- lapply(1:n_patients, function(i) time[[i]])
  h_save <- lapply(1:n_patients, function(i){
    out <- matrix(NA, nrow = nrun, ncol = n_obs[i])
    out[1,] <- time[[i]]
    return(out)
  })
  
  
  
  X <- stack_Matrix(Bm)
  
  diag_m <- lapply(1:n_patients, function(i) diag(n_obs[i]))
  
  
  
  #----- MCMC Loop -----#
  
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(iter in 2:nrun){
    
    setTxtProgressBar(bar, iter)
    
    
    #-- Update Phi --#
    
    # tmp_Bm<-list()
    # for(i in 1:n_patients){
    #   tmp_phi <- phi[[i]]
    #   current_llik <- mvnfast::dmvn(X = as.numeric(y_list[[i]]), mu =  as.numeric(c_save[iter-1,i]+a_save[iter-1,i]*Bm[[i]] %*% beta_save[iter - 1,]),
    #                        sigma = diag_m[[i]]* sigma_eps_save[iter - 1] , log = TRUE)
    #   current_lprior <- mvnfast::dmvn(X = phi[[i]], mu = Upsilon,
    #                          sigma = sigma_phi_save[iter - 1] * inv_P, log = TRUE)
    #   for(j in 2:(q-1)){
    #     tmp_phi[j] <- Arma_runif(max(tmp_phi[j] - tune, tmp_phi[j-1]),min(tmp_phi[j] + tune, tmp_phi[j+1]))
    #     tmp_Bm[[i]] <- bs(Bh[[i]] %*% tmp_phi, knots = knots_m, intercept = TRUE)
    # 
    #     cand_llik <- mvnfast::dmvn(X = as.numeric(y_list[[i]]), mu =  as.numeric(c_save[iter-1,i]+a_save[iter-1,i]*tmp_Bm[[i]] %*% beta_save[iter - 1,]),
    #                       sigma = sigma_eps_save[iter - 1] * diag_m[[i]], log = TRUE)
    #     cand_lprior <- mvnfast::dmvn(X = tmp_phi, mu = Upsilon,
    #                         sigma = sigma_phi_save[iter - 1] * inv_P, log = TRUE)
    #     lratio <- cand_llik + cand_lprior - current_llik - current_lprior
    # 
    #     if(log(Arma_runif(0, 1)) < lratio){
    #       current_llik <- cand_llik
    #       current_lprior <- cand_lprior
    #       accepts[i, j-1]= accepts[i, j-1]+1
    #     }
    #     else{
    #       tmp_phi[j] <- phi[[i]][j]
    #     }
    #   }
    #   phi[[i]] <- tmp_phi
    #   phi_save[iter,,i]=tmp_phi
    # 
    # }
    
    for (i in 1:n_patients){


      phi_old = phi[[i]]
      phi_new=phi_old

      for (b in 2:(q-1)){


        pi_phi_old = pi_phi(phi_old , beta_save[iter-1,],  sigma_eps_save[iter-1], a_save[iter-1, i],
                            c_save[iter-1,i], P, sigma_phi_save[iter-1], Upsilon, p,  Bh[[i]],
                            as.numeric(y_list[[i]]), n_obs[i])

        phi_new[b] <- Arma_runif(max(phi_new[b] - tune, phi_new[b-1]),min(phi_new[b] + tune, phi_new[b+1]))

        pi_phi_new = pi_phi(phi_new, beta_save[iter-1,], sigma_eps_save[iter-1],a_save[iter-1, i],
                            c_save[iter-1,i], P, sigma_phi_save[iter-1], Upsilon,p,  Bh[[i]],
                            as.numeric(y_list[[i]]), n_obs[i])

        alpha = min(0, (pi_phi_new-pi_phi_old))

        u = Arma_runif(0, 1)

        if(u<exp(alpha)){
          #accept
          accepts[i, b-1]= accepts[i, b-1]+1
        } else{
          phi_new[b]= phi_old[b]
        }
      }
      phi[[i]]=phi_new
      phi_save[iter,,i]=phi_new
    }



    #-- Update h and h_save --#
    
    h <- lapply(1:n_patients, function(i) as.numeric(Bh[[i]] %*% phi[[i]]))
    for(i in 1:n_patients){
      h_save[[i]][iter,] <- h[[i]]
    }
    
    #-- Update Bm --#
    
    Bm <- lapply(1:n_patients, function(i) bs(h[[i]], knots = knots_m, intercept = TRUE))
    
    #-- Update X --#
    
    X <- stack_Matrix(Bm)
    
    #-- Update Beta --#
    
    inv_sigma_beta = Omega/lambda_save[iter-1]
    inv_V_beta = inv_sigma_beta + (1/sigma_eps_save[iter-1])*crossprod(X)
    
    V_beta = armaInv(inv_V_beta)
    m_beta =V_beta%*%(t(X) %*% (1/sigma_eps_save[iter-1]*y) + inv_sigma_beta %*% beta_0 )
    beta_save[iter, ] =  Arma_mvrnorm(m_beta, V_beta)
    
    
    ##-- Update a_0,c_0 --#
    
    sigma_star = 1/(1/sigma_a0 + n_patients/sigmaa_save[iter-1])
    a_star = sigma_star*(sum(a_save[iter-1, ])/sigmaa_save[iter-1] + m_a0/sigma_a0)
    a_0_save[iter] = rnorm(1,a_star, sqrt(sigma_star))
    
    sigma_star = 1/(1/sigma_c0 + n_patients/sigmac_save[iter-1])
    c_star = sigma_star*(sum(c_save[iter-1, ])/sigmac_save[iter-1] + m_c0/sigma_c0)
    c_0_save[iter] =rnorm(1,c_star, sqrt(sigma_star))
    
    
    # #-- Update c_i, a_i --#
    
    for (i in 1:n_patients){
      
      sigma_ca = matrix(0, nrow = 2, ncol=2)
      sigma_ca[1,1] = sigmac_save[iter-1]
      sigma_ca[2,2] = sigmaa_save[iter-1]
      inv_sigma_ca = armaInv(sigma_ca)
      W = cbind(rep(1, n_obs[i]), Bm[[i]]%*% beta_save[iter,])
      inv_sigma_l = inv_sigma_ca + (1/sigma_eps_save[iter-1])*crossprod(W)
      sigma_l = armaInv(inv_sigma_l)
      mu_l = sigma_l%*%(inv_sigma_ca%*%c(c_0_save[iter], a_0_save[iter]) + 
                          (1/sigma_eps_save[iter-1])*crossprod(W,y_list[[i]]))
      vec =Arma_mvrnorm(mu_l, sigma_l)
      c_save[iter, i] = vec[1]
      a_save[iter, i] = vec[2]
      
    }
    
    
    #-- Update sigma eps --#
    
    a_star = a_eps + 0.5*sum(n_obs)
    m_tilde = list()
    for (pat in 1:n_patients){
      m_tilde[[pat]] = c_save[iter, pat]*rep(1, n_obs[pat]) + 
        a_save[iter, pat]*Bm[[pat]]%*%beta_save[iter,]
    }
    sum_b = 0
    for(pat in 1:n_patients){
      sum_b = sum_b + crossprod(y_list[[pat]]-m_tilde[[pat]])
    }
    b_star = as.numeric(b_eps + 0.5*sum_b)
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    
    #-- Update sigma_c --#
    
    a_star = a_c + 0.5*n_patients
    b_star = b_c + 0.5*sum((c_save[iter, ]-c_0_save[iter])^2)
    sigmac_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    #-- Update sigma_a --#
    
    a_star = a_a + 0.5*n_patients
    b_star = b_a + 0.5*sum((a_save[iter, ]-a_0_save[iter])^2)
    sigmaa_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    #-- Update lambda --#
    
    a_star=0.5*p+a_lambda
    b_star=0.5 * as.numeric(t(beta_save[iter,]) %*% Omega %*% (beta_save[iter,])) +  b_lambda
    lambda_save[iter]<-1/rgamma(n=1, shape=a_star, rate=b_star)
    
    #-- Update sigma_phi --#
    
    a_star = a_phi + 0.5*q*n_patients
    b_star = b_phi + 0.5*as.numeric(t(unlist(phi)- rep(Upsilon, n_patients)) %*% bigP %*% (unlist(phi)- rep(Upsilon, n_patients)))
    sigma_phi_save[iter] = 1/rgamma(n=1, shape=a_star, rate=b_star)
    
    
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
  h_p <- lapply(h_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  Bm_post <- lapply(h_p, function(w)  bs(w, knots = knots_m, intercept = T))
  y_p <- lapply(1:n_patients, function(i) as.numeric(c_post[i]+a_post[i]*Bm_post[[i]] %*% beta_post))
  y_reg <- lapply(1:n_patients, function(i) interp_spline(h_p[[i]], y_list[[i]]))
  
  # ---List conversion--- #
  
  h_post=matrix(NA,nrow=n_obs_max,ncol=n_patients)
  
  for(i in 1:n_patients){
    h_post[1:n_obs[i],i]=h_p[[i]]
  }
  
  y_post=matrix(NA,nrow=n_obs_max,ncol=n_patients)
  
  for(i in 1:n_patients){
    y_post[1:n_obs[i],i]=y_p[[i]]
  }
  
  y_star=matrix(NA,nrow=n_obs_max,ncol=n_patients)
  
  for(i in 1:n_patients){
    y_star[1:n_obs[i],i]=y_reg[[i]]
  }
  
  return (list (post=list(a = a_post, c=c_post, c_0 = c_0_post, sigmac = sigmac_post, a_0 = a_0_post, 
                          sigmaa = sigmaa_post, beta = beta_post, lambda = lambda_post, phi=phi_post,
                          sigma_phi = sigma_phi_post, sigma_eps = sigma_eps_post, accepts = accepts,h=h_post ,y=y_post), 
                full=list(a_s=a_save, c_s = c_save, c0_s=c_0_save, sigmac_s=sigmac_save, 
                          a0_s=a_0_save, sigmaa_s=sigmaa_save, beta_s=beta_save,
                          lambda_s=lambda_save,phi_s=phi_save, sigma_phi_s=sigma_phi_save,
                          sigma_eps_s=sigma_eps_save),
                y_star=y_star))
}

#### Run the function ####
m_c0=0
sigma_c0=10
m_a0=5 
sigma_a0=10
a_c=250
b_c=600
a_a=25
b_a=100
a_lambda=1200
b_lambda=3500
a_phi=200
b_phi=600
a_eps=3000
b_eps=5000
n_knots_h=3
n_knots_m=20

# graphics.off()
# a_star = 2.5 + 22
# x11()
# invgamma = 1/rgamma(n=1000, 100, 1000)
# hist(invgamma)
# median(invgamma)
# mean(invgamma)

#y= read.csv('Data_Undersampled_surgery.csv') #qua caricate il dataset con l'undersampling

nburn = 50000
niter = 15000

# 1 ora 14%
#19:11-22:01
# 15 min 10%
start.time <- Sys.time()
telesca_output = telesca_model(y=y, niter = niter, nburn = nburn, n_knots_h = n_knots_h, 
                             n_knots_m = n_knots_m,
                             m_c0=m_c0, sigma_c0=sigma_c0, a_c=a_c, b_c=b_c, m_a0=m_a0, 
                             sigma_a0=sigma_a0, a_a=a_a, b_a=b_a,
                             a_eps = a_eps, b_eps = b_eps, a_lambda = a_lambda, b_lambda = b_lambda, 
                             a_phi = a_phi, b_phi =b_phi)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,5)
print(time.taken)
# 1 ora per 2000 iterazioni con 10 pazienti, 8 nodi h, 40 nodi m (dati senza undersampling)
# 2.8 minuti per 2000 iterazioni con 10 pazienti, 8 nodi h, 40 nodi m, con undersamplin (1 ogni 4)
# 12.8 minuti per 2000 iterazioni con 32 pazienti, 8 nodi h, 40 nodi m, con undersampling (1 ogni 4)
# 40 min 32 pazienti 8000 iterazioni

# load("telesca_surgery.RData")
# load("Telesca_fisio_def.Rdata")
# load("runnata_sani_22_12.Rdata")
# load("telesca_complete_Simo_v2.Rdata")

# save the output of the model
acceptance_output = telesca_output$post$accepts
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


#### Acceptance rate + trace plots ####
x11()
hist(acceptance_output, main='Acceptance rate')
x11()
par(mfrow=c(2,1))
matplot(a_output[(nburn+1):(niter+nburn),], type='l', main='a_i', xlab='Iteration')
matplot(c_output[(nburn+1):(niter+nburn),], type='l', main='c_i', xlab='Iteration')

x11()
par(mfrow=c(2,1))
matplot(a0_output[(nburn+1):(niter+nburn)], type='l', main='a0', xlab='Iteration')
matplot(c0_output[(nburn+1):(niter+nburn)], type='l', main='c0', xlab='Iteration')

x11()
par(mfrow=c(2,1))
matplot(sigma_a_output[(nburn+1):(niter+nburn)], type='l', main='sigma^2_a0', xlab='Iteration')
matplot(sigma_c_output[(nburn+1):(niter+nburn)], type='l', main='sigma^2_c0', xlab='Iteration')

x11()
matplot(beta_output[(nburn+1):(niter+nburn),],type='l', xlab = 'Iteration', main='beta')

x11()
matplot(lambda_output[(nburn+1):(niter+nburn)],type='l', xlab = 'Iteration', main='lambda')

x11()
par(mfrow=c(2,1))
matplot(phi_output[(nburn+1):(niter+nburn),,9],type='l',xlab = 'Iteration', main='phi of patient 9')
matplot(phi_output[(nburn+1):(niter+nburn),,16],type='l',xlab = 'Iteration', main='phi of patient 16')


x11()
par(mfrow=c(2,1))
matplot(sigma_phi_output[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='sigma^2_phi')
matplot(sigma_eps_output[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='sigma^2_eps')


#graphics.off()



#### Trace plots with burn-in ####
x11()
par(mfrow=c(2,1))
matplot(a_output, type='l', main='a_i', xlab='Iteration')
matplot(c_output, type='l', main='c_i', xlab='Iteration')

x11()
par(mfrow=c(2,1))
matplot(a0_output, type='l', main='a0', xlab='Iteration')
matplot(c0_output, type='l', main='c0', xlab='Iteration')

x11()
par(mfrow=c(2,1))
matplot(sigma_a_output, type='l', main='sigma^2_a0', xlab='Iteration')
matplot(sigma_c_output, type='l', main='sigma^2_c0', xlab='Iteration')

x11()
matplot(beta_output,type='l', xlab = 'Iteration', main='beta')

x11()
matplot(beta_output[,15],type='l', xlab = 'Iteration', main='beta')

x11()
matplot(lambda_output,type='l', xlab = 'Iteration', main='lambda')


x11()
par(mfrow=c(2,1))
matplot(phi_output[,,9],type='l',xlab = 'Iteration', main='phi of patient 9')
matplot(phi_output[,,16],type='l',xlab = 'Iteration', main='phi of patient 16')


x11()
par(mfrow=c(2,1))
matplot(sigma_phi_output, type='l',xlab = 'Iteration', 
        main='sigma^2_phi')
matplot(sigma_eps_output, type='l',xlab = 'Iteration', 
        main='sigma^2_eps')



# fancier traceplot
prova = mcmc(a0_output)
x11()
xyplot(ts(prova))

#graphics.off()




#### Plot of the curves ####
n_obs_max  = dim(y)[1] 
n_patients = dim(y)[2] # y has patients in the columns
y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))
n_obs=sapply(y_list,length)

time = matrix(NA, nrow = n_obs_max, ncol=n_patients)
for (i in 1:n_patients){
  time[1:n_obs[i],i] = seq(0,1, length.out = n_obs[i])
}

knots_h = seq(0, 1, length.out = n_knots_h+2)[-c(1,n_knots_h+2)]

time_grid = seq(0,1, length.out=1000)

h_i = matrix(NA, nrow=length(time_grid), ncol=n_patients )
for (i in 1:n_patients){
  h_i[,i] = b_splines(n_knots_h+4, time_grid, coeff_vec=telesca_output$post$phi[,i])
}

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




col_pat = rainbow(n_patients)
ylim_plot = c(0,90) #set ylim for the plot

x11()
par(mfrow=c(1,3))
matplot(y, ylim=ylim_plot, type='l', main='y', col = col_pat)
plot(seq(0,1, len=n_obs[1]),y[1:n_obs[1],1], type='l', main='y', ylim=ylim_plot,
    col=col_pat[1], xlab='time')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),y[1:n_obs[i],i], type='l', col=col_pat[i])
}
plot(seq(0,1, len=n_obs[1]),telesca_output$y_star[1:n_obs[1],1], type='l', main='y*',
     ylim=ylim_plot, col=col_pat[1], xlab='time')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),telesca_output$y_star[1:n_obs[i],i], type='l', col=col_pat[i])
}

x11()
plot(seq(0,1, len=n_obs[1]),y[1:n_obs[1],1], type='l', main='Observed curves', ylim=ylim_plot,
     col=col_pat[1], xlab='t', ylab='y(t)', cex.main=1.5, cex.lab=1.5, cex.axis=1.5)
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),y[1:n_obs[i],i], type='l', col=col_pat[i])
}

x11()
plot(seq(0,1, len=n_obs[1]),telesca_output$y_star[1:n_obs[1],1], type='l', main='Aligned curves',
     ylim=ylim_plot, col=col_pat[1], xlab='t', ylab='y*(t)')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),telesca_output$y_star[1:n_obs[i],i], type='l', col=col_pat[i])
}


#graphics.off()
x11()
matplot(time_grid, h_i, type='l', main='h_i', xlab = 't')


#### If you want to compute the execution time of a chunk ####
start.time <- Sys.time()
end.time <- Sys.time()
time.taken <- round(end.time - start.time,5)
print(time.taken)



#### Diagnostics  - create MCMC objects ####
a_mcmc = mcmc(a_output[-c(1:nburn),])
c_mcmc = mcmc(c_output[-c(1:nburn),])
a0_mcmc = mcmc(a0_output[-c(1:nburn)])
c0_mcmc = mcmc(c0_output[-c(1:nburn)])
sigma_a_mcmc = mcmc(sigma_a_output[-c(1:nburn)])
sigma_c_mcmc = mcmc(sigma_c_output[-c(1:nburn)])
beta_mcmc = mcmc(beta_output[-c(1:nburn),])
lambda_mcmc = mcmc(lambda_output[-c(1:nburn)])
sigma_phi_mcmc = mcmc(sigma_phi_output[-c(1:nburn)])
sigma_eps_mcmc = mcmc(sigma_eps_output[-c(1:nburn)])

#### Diagnostics - a ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(a_mcmc)
s$statistics[,4] #MC standard errors

# ESS
ess_a = effectiveSize(a_mcmc)
median(ess_a)

# Autocorrelation plot
x11()
#acf(a_mcmc, lag.max=500)

# graphics.off()




#### Diagnostics - c ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(c_mcmc)
s$statistics[,4] #MC standard errors

# ESS
ess_c = effectiveSize(c_mcmc)
median(ess_c)

# Autocorrelation plot
x11()
#acf(c_mcmc, lag.max=500)

# graphics.off()

#### Diagnostics - a0 ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(a0_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(a0_mcmc)

#Autocorrelation plot
x11()
acf(a0_mcmc, lag.max=500)


#### Diagnostics - c0 ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(c0_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(c0_mcmc)

#Autocorrelation plot
x11()
acf(c0_mcmc, lag.max=500)





#### Diagnostics - sigma_a ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(sigma_a_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(sigma_a_mcmc)

#Autocorrelation plot
x11()
acf(sigma_a_mcmc, lag.max=500)









#### Diagnostics - sigma_c ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(sigma_c_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(sigma_c_mcmc)

#Autocorrelation plot
x11()
acf(sigma_c_mcmc, lag.max=500)


#### Diagnostics - beta ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(beta_mcmc)
s$statistics[,4] #MC standard errors

# ESS
ess_beta = effectiveSize(beta_mcmc)
median(ess_beta)

# Autocorrelation plot
# x11()
# autocorr.plot(sigma_eps_mcmc, lag.max=500)

x11()
acf(beta_mcmc[,22], lag.max=500)

# graphics.off()



#### Diagnostics - lambda ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(lambda_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(lambda_mcmc)

#Autocorrelation plot
x11()
acf(lambda_mcmc, lag.max=500)


#### Diagnostics - sigma phi ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(sigma_phi_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(sigma_phi_mcmc)


#Autocorrelation plot
x11()
acf(sigma_phi_mcmc, lag.max=500)

#### Diagnostics - sigma eps ####
# batchSE(prova, batchSize = 50) # da capire
s=summary(sigma_eps_mcmc)
s$statistics[4] #MC standard error

# ESS
effectiveSize(sigma_eps_mcmc)


#Autocorrelation plot
x11()
acf(sigma_eps_mcmc, lag.max=500)


#### Diagnostics - phi ####
ess_phi = matrix(NA, nrow = dim(y)[2], ncol=n_knots_h+2)
sd_error = matrix(NA, nrow = dim(y)[2], ncol=n_knots_h+2)
for (i in 1:dim(y)[2]){
  phi_temp = mcmc(phi_output[-c(1:nburn),-c(1,7),i])
  ess_phi[i,] = effectiveSize(phi_temp)
  sd_error[i,] = summary(phi_temp)$statistics[,4]
}

median(ess_phi)
median(sd_error)

