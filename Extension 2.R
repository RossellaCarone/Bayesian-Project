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

#### Rcpp functions and auxiliary functions ####
Rcpp::cppFunction("arma::mat armaInv(const arma::mat & x) { return arma::inv(x); }", depends="RcppArmadillo")
Rcpp::cppFunction("arma::vec Arma_mvrnorm(const arma::vec& mean, const arma::mat& cov) {
  return arma::mvnrnd(mean, cov);}", depends="RcppArmadillo")
Rcpp::cppFunction("double  Arma_rnorm(const double &mu, const double &sigma){
      return arma::randn(arma::distr_param(mu,sigma));} ", depends="RcppArmadillo") #variance = sigma^2
Rcpp::cppFunction("double  Arma_runif(const double &a, const double &b){
      return arma::randu(arma::distr_param(a,b));} ", depends="RcppArmadillo")



pi_phi = function (phi, theta,sigma_eps, P, sigma_phi, Upsilon, p, B_h, y_i){
  warped_t = B_h%*%phi
  mtilde_i = as.numeric(b_splines(p, warped_t, theta))
  add1 = y_i - mtilde_i
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


interp_spline <- function(x, y, nout = length(y)) {
  ind_out <- seq(min(x), max(x), len = nout)
  spfit <- splinefun(x, y)
  return(spfit(ind_out))
}


#### Extension 2 ####
extension2 <- function(y, niter = 1000, nburn = 1000, n_knots_h = 4, n_knots_m = 4,
                       a_eps = .5, b_eps = .2, a_phi = .1, b_phi = .1,v_psi=5,a1=1.5,a2=1.5,a_sigma=.5,b_sigma=0.25,
                       epsilon=1e-2,alpha_0=1.5,alpha_1=5*1e-4
){ 
  
  prop = 1.00
  num=0
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
  P = omega_P(q)
  inv_P = armaInv(P)
  k_tilde<-floor(log(p)*4)
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
  
  tune <- .005
  accepts = matrix(0, nrow=n_patients, ncol=q-2)

  
  #----- Save Structures -----#
  
  nrun <- nburn + niter
  k_star = k_tilde
  k_save=numeric(nrun)
  k_save[1]=k_star
  
  psi_save=list()
  psi_save[[1]]=matrix(rgamma(p*k_star,v_psi/2,v_psi/2),p,k_star)
  
  delta_save=list()
  delta_save[[1]]=rep(0,k_star)
  
  delta_save[[1]][1]=rgamma(1,a1,1)
  delta_save[[1]][2:k_star]=rgamma(k_star-1,a2,1)
  
  tau_save=list()
  tau_save[[1]]=cumprod(delta_save[[1]])
  
  eta_save=list()
  eta_save[[1]]=matrix(0,n_patients,k_star)
  for (i in 1:n_patients){
    eta_save[[1]][i,]=Arma_mvrnorm(rep(0,k_star),diag(k_star))
  }
  
  sigma_lambda_save=list()
  sigma_lambda_save[[1]]=1/rgamma(p,a_sigma,b_sigma)
  
  Ptht= psi_save[[1]] * matrix(rep(tau_save[[1]],p),p,k_star,byrow=T)
  
  lambda_save=list()
  lambda_save[[1]]=matrix(0,p,k_star)
  
  #Recall that: theta_i = Lambda_matrix%*%eta_i + gamma_i
  
  Sigma=diag(sigma_lambda_save[[1]])
  Sigma=(Sigma + t(Sigma))/2
  Sigmainv=diag(1/sigma_lambda_save[[1]])
  Sigmainv=(Sigmainv + t(Sigmainv))/2

  
  theta_save=list()  
  theta_save[[1]]=matrix(0,n_patients,p)
  for (i in 1:n_patients){
    theta_save[[1]][i,]=Arma_mvrnorm(lambda_save[[1]]%*%eta_save[[1]][i,],Sigma)
  }
  
  sigma_eps_save <- numeric(nrun)
  sigma_eps_save[1] <- 1 
  
  
  # warping parameters
  
  sigma_phi_save <- numeric(nrun)
  sigma_phi_save[1] <- 1
  
  phi <- lapply(1:n_patients, function(i) Upsilon)
  
  # tensor with all phi_i
  phi_save = array(rep(NA, nrun *q*n_patients), dim=c(nrun , q , n_patients))
  for (i in 1:n_patients){
    phi_save[1, ,i] = Upsilon
  }
  
  h <- lapply(1:n_patients, function(i) time[[i]])
  h_save <- lapply(1:n_patients, function(i){
    out <- matrix(NA, nrow = nrun, ncol = n_obs[i])
    out[1,] <- time[[i]]
    return(out)
  })

  diag_m <- lapply(1:n_patients, function(i) diag(n_obs[i]))
  
  
  #----- MCMC Loop -----#
  
  bar <- txtProgressBar(min = 2, max = nrun, style = 3)
  for(iter in 2:nrun){
    
    setTxtProgressBar(bar, iter)
    
    #--------------------------------SMOOTHING----------------------------------
    
    #-- Update lambda --#
    
    lambda_save[[iter]]=matrix(0,p,k_star)
    
    for(j in 1:p){
      
      Vlam = diag(Ptht[j,])+(1/sigma_lambda_save[[iter-1]][j])*t(eta_save[[iter-1]])%*%eta_save[[iter-1]]
      Vlam = (Vlam + t(Vlam))/2
      Vlam= armaInv(Vlam)
      Sigma_lambdaj =(Vlam + t(Vlam))/2
      Elam= Sigma_lambdaj %*% t(eta_save[[iter-1]])%*% theta_save[[iter-1]][,j]/sigma_lambda_save[[iter-1]][j]
      lambda_save[[iter]][j,]=Arma_mvrnorm(Elam,Sigma_lambdaj)
      
    }
    
    
    #-- Update Psi --#
    
    psi_save[[iter]]=matrix(0,p,k_star)
    
    for(j in 1:p){
      for(h in 1:k_star){
        psi_save[[iter]][j,h]=rgamma(1,(v_psi+1)/2, (tau_save[[iter-1]][h]*(lambda_save[[iter]][j,h]^2))/2+v_psi/2)
      }
    }
    
    
    #-- Update Delta --#
    
    delta_save[[iter]]=rep(0,k_star)
    term=0
    for (l in 2:k_star){
      tau_l1=prod(delta_save[[iter-1]][2:l])
      term=term+tau_l1*psi_save[[iter]][,l]%*%(lambda_save[[iter]][,l]^2)
    }
    delta_save[[iter]][1]=rgamma(1,a1+p*k_star/2,1+0.5*term)
    
    for (h in 2:k_star){
      term=0
      for (l in h:k_star){
        d=c(delta_save[[iter]][1:(h-1)],delta_save[[iter-1]][h:l])
        d=d[-h]
        tau_lh=prod(d)
        term=term+tau_lh*psi_save[[iter]][,l]%*%(lambda_save[[iter]][,l]^2)
      }
      delta_save[[iter]][h]=rgamma(1,a2+p*(k_star-h+1)/2,1+0.5*term)
    }
    
    #-- Update Tau --#
    
    tau_save[[iter]]=cumprod(delta_save[[iter]])
    
    #--Update Ptht --#
    
    Ptht=psi_save[[iter]] * matrix(rep(tau_save[[iter]],p),p,k_star,byrow=T)
    
    #-- Update sigma lambda --#
    
    sigma_lambda_save[[iter]]=rep(0,p)
    
    
    for (j in 1:p){
      term=0
      for (i in 1:n_patients)  {
        term=term + (theta_save[[iter-1]][i,j]-t(lambda_save[[iter]][j,])%*%eta_save[[iter-1]][i,])^2
      }
      
      sigma_lambda_save[[iter]][j]=1/rgamma(1,n_patients/2 + a_sigma,b_sigma+0.5*term)
    }
    
    
    #--Update Sigma --#
    
    Sigma=diag(sigma_lambda_save[[iter]])
    Sigma=(Sigma + t(Sigma))/2
    Sigmainv=diag(1/sigma_lambda_save[[iter]])
    Sigmainv=(Sigmainv + t(Sigmainv))/2
    
    #-- Update eta --#
    
    eta_save[[iter]]=matrix(0,n_patients,k_star)
    
    
    for(i in 1:n_patients){
      
      term1=Bm[[i]]%*%lambda_save[[iter]]
      term2=(sigma_eps_save[iter-1])*diag(n_obs[i])+Bm[[i]]%*%Sigma%*%t(Bm[[i]])
      term2=(term2+t(term2))/2
      term2=armaInv(term2)
      term2=(term2+t(term2))/2
      A =diag(k_star)+ t(term1)%*%term2%*%term1
      C=t(term1)%*%term2%*%y_list[[i]]
      inv_A=armaInv(A)
      inv_A=(inv_A+t(inv_A))/2
      eta_save[[iter]][i,]=Arma_mvrnorm(inv_A%*%C,inv_A)
    }
    
    #-- Update theta --#
    
    theta_save[[iter]]=matrix(0,n_patients,p)
    
    for (i in 1:n_patients){
      term1=t(Bm[[i]])%*%Bm[[i]]/sigma_eps_save[iter-1]+Sigmainv
      term1=(term1+t(term1))/2
      term1=armaInv(term1)
      term1=(term1+t(term1))/2
      term2=t(Bm[[i]])%*%y_list[[i]]/sigma_eps_save[iter-1]+Sigmainv%*%lambda_save[[iter]]%*%eta_save[[iter]][i,]
      theta_save[[iter]][i,]=Arma_mvrnorm(term1%*%term2,term1)
    }
    
    
    #-- Update k_star--#
    
    prob=1/exp(alpha_0+alpha_1*(iter-1))
    uu=Arma_runif(0, 1)
    lind=(colSums(abs(lambda_save[[iter]])<epsilon))/p
    vec=lind>=prop
    num=sum(vec)
    
    if(uu<min(prob,1)){
      if((iter-1)>20 & num==0 & all(lind<0.995)){
        k_star=k_star+1
        lambda_save[[iter]]=cbind(lambda_save[[iter]],rep(0,p))
        eta_save[[iter]]=cbind(eta_save[[iter]],Arma_mvrnorm(rep(0,n_patients),diag(n_patients)))
        psi_save[[iter]]=cbind(psi_save[[iter]],rgamma(p,v_psi/2,v_psi/2))
        delta_save[[iter]]=c(delta_save[[iter]],rgamma(1,a2,1))
        tau_save[[iter]]=cumprod(delta_save[[iter]])
        Ptht=psi_save[[iter]] * matrix(rep(tau_save[[iter]],p),p,k_star,byrow=T)
      }
      else{
        if(num>0){
          nonred <- setdiff(1:k_star, which(vec))
          if(num==k_star){
            k_star=1
            lambda_save[[iter]]=rep(0,p)
            eta_save[[iter]]=Arma_mvrnorm(rep(0,n_patients),diag(n_patients))
            psi_save[[iter]]=matrix(rgamma(p*k_star,v_psi/2,v_psi/2),p,k_star)
            delta_save[[iter]]=rgamma(1,a1,1)
            tau_save[[iter]]=delta_save[[iter]]
            Ptht=psi_save[[iter]] * matrix(rep(tau_save[[iter]],p),p,k_star,byrow=T)
            theta_save[[iter]]=matrix(0,n_patients,p)
            for (i in 1:n_patients){
              theta_save[[1]][i,]=Arma_mvrnorm(lambda_save[[iter]]%*%eta_save[[iter]][i,],Sigma)
            }
          }
          else{
            k_star = max(k_star - num,1)
            lambda_save[[iter]]=lambda_save[[iter]][,nonred]
            psi_save[[iter]]=psi_save[[iter]][,nonred]
            eta_save[[iter]]=eta_save[[iter]][,nonred]
            delta_save[[iter]]=delta_save[[iter]][nonred]
            tau_save[[iter]]=cumprod(delta_save[[iter]])
            Ptht=psi_save[[iter]] * matrix(rep(tau_save[[iter]],p),p,k_star,byrow=T)
          }
        }
      }
    }
    k_save[iter]=k_star
    
    #-------------------------------WARPING-------------------------------------
    
    #-- Update Phi --#
    
    for (i in 1:n_patients){
      
      
      phi_old = phi[[i]]
      phi_new=phi_old
      
      for (b in 2:(q-1)){
        
        
        pi_phi_old = pi_phi(phi_old ,theta_save[[iter]][i,],sigma_eps_save[iter-1], 
                            P, sigma_phi_save[iter-1], Upsilon, p,  Bh[[i]],
                            as.numeric(y_list[[i]]))
        
        phi_new[b] <- Arma_runif(max(phi_new[b] - tune, phi_new[b-1]),min(phi_new[b] + tune, phi_new[b+1]))
        
        pi_phi_new = pi_phi(phi_new,theta_save[[iter]][i,], sigma_eps_save[iter-1],
                            P, sigma_phi_save[iter-1], Upsilon,p,  Bh[[i]],
                            as.numeric(y_list[[i]]))
        
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
    
    #-- Update sigma eps --#
    
    a_star = a_eps + 0.5*sum(n_obs)
    m_tilde = list()
    for (pat in 1:n_patients){
      m_tilde[[pat]]=Bm[[pat]]%*%theta_save[[iter]][pat,]
    }
    sum_b = 0
    for(pat in 1:n_patients){
      sum_b = sum_b + crossprod(y_list[[pat]]-m_tilde[[pat]])
    }
    b_star = as.numeric(b_eps + 0.5*sum_b)
    sigma_eps_save[iter] = 1/rgamma(n=1,a_star, b_star)
    
    
    
    #-- Update sigma phi --#
    
    a_star = a_phi + 0.5*q*n_patients
    b_star = b_phi + 0.5*as.numeric(t(unlist(phi)- rep(Upsilon, n_patients)) %*% bigP %*% (unlist(phi)- rep(Upsilon, n_patients)))
    sigma_phi_save[iter] = 1/rgamma(n=1, shape=a_star, rate=b_star)
    
  }
  
  
  #-------------------------------RETURN----------------------------------------
  
  
  close(bar)
  accepts <- accepts / nrun
  theta_post<-Reduce('+',theta_save[-c(1:nburn)])/niter
  phi_post<-apply(phi_save[-c(1:nburn),,], c(2, 3), mean)
  sigma_phi_post<-mean(sigma_phi_save[-c(1:nburn)])
  sigma_eps_post<-mean(sigma_eps_save[-c(1:nburn)])
  sigma_lambda_post<-Reduce('+',sigma_lambda_save[-c(1:nburn)])/niter
  h_p <- lapply(h_save, function(w) apply(w[-c(1:nburn),], 2, mean))
  Bm_post <- lapply(h_p, function(w)  bs(w, knots = knots_m, intercept = T))
  y_p <- lapply(1:n_patients, function(i) as.numeric(Bm_post[[i]] %*% theta_post[i,]))
  m_i <- lapply(1:n_patients, function(i)  as.numeric(b_splines(p,time[[i]],theta_post[i,])))
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
  
  mi=matrix(NA,nrow=n_obs_max,ncol=n_patients)
  
  for(i in 1:n_patients){
    mi[1:n_obs[i],i]=m_i[[i]]
  }
  
  y_star=matrix(NA,nrow=n_obs_max,ncol=n_patients)
  
  for(i in 1:n_patients){
    y_star[1:n_obs[i],i]=y_reg[[i]]
  }
  
  return (list (post=list(phi=phi_post,sigmaj = sigma_lambda_post,sigma_phi = sigma_phi_post, sigma_eps = sigma_eps_post, accepts = accepts,h=h_post ,mi_tilde=y_post,mi=mi), 
                full=list(k_s=k_save,lambda_s=lambda_save,psi_s=psi_save,eta_s=eta_save,theta_s=theta_save,delta_s=delta_save,tau_s=tau_save,sigma_j_s=sigma_lambda_save,phi_s=phi_save, sigma_phi_s=sigma_phi_save, sigma_eps_s=sigma_eps_save),
                y_star=y_star))
}

# library(fdasrvf)
# y<-growth_vel$f[,1:20]
# result<-extension2(y)

x11()
par(mfrow=c(1,4))
matplot(y, type='l', main='y')
matplot(result$post$mi_tilde, type='l', main='mi tilde')
matplot(result$post$mi, type='l', main='mi')
matplot(result$y_star, type='l', main='y*')


#### Run the function ####
a_phi=200
b_phi=600
a_eps=3000
b_eps=5000
n_knots_h=3
n_knots_m=36

#y= read.csv('Data_Undersampled_healthy.csv') #qua caricate il dataset con l'undersampling
#y = read.csv('healthy_new_undersampling.csv')

nburn = 50000
niter = 15000
start.time <- Sys.time()
result = extension2(y=y, n_knots_h = n_knots_h, n_knots_m = n_knots_m, a_eps = a_eps, b_eps = b_eps,
                    a_phi = a_phi, b_phi = b_phi, nburn = nburn, niter =niter)
end.time <- Sys.time()
time.taken <- round(end.time - start.time,5)
print(time.taken)



#### Plot curves ####
# mettere cex.lab = 1.2 e cex.axis=1.2 per ingrandire i valori sugli assi!!
n_obs_max  = dim(y)[1] 
n_patients = dim(y)[2] # y has patients in the columns
y_list = lapply(seq_len(ncol(y)), function(i) na.omit(y[, i]))
n_obs=sapply(y_list,length)

col_pat = rainbow(n_patients)
ylim_plot = c(0,90) #set ylim for the plot

x11()
plot(seq(0,1, len=n_obs[1]),y[1:n_obs[1],1], type='l', main='y', ylim=ylim_plot,
     col=col_pat[1], xlab='time')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),y[1:n_obs[i],i], type='l', col=col_pat[i])
}

x11()
plot(seq(0,1, len=n_obs[1]),result$y_star[1:n_obs[1],1], type='l', main='y*',
     ylim=ylim_plot, col=col_pat[1], xlab='time')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]),result$y_star[1:n_obs[i],i], type='l', col=col_pat[i])
}


x11()
plot(seq(0,1, len=n_obs[1]),result$post$mi[1:n_obs[1],1], type='l', main='m_i', ylim=ylim_plot,
     col=col_pat[1], xlab='time')
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]), result$post$mi[1:n_obs[i],i], type='l', col=col_pat[i])
}

x11()
plot(seq(0,1, len=n_obs[1]),result$post$mi_tilde[1:n_obs[1],1], type='l', ylab='m_i[h_i(t)]', ylim=ylim_plot,
     col=col_pat[1], xlab='t', cex.axis=1.2, cex.lab=1.2)
for(i in 2:n_patients){
  points(seq(0,1, len=n_obs[i]), result$post$mi_tilde[1:n_obs[i],i], type='l', col=col_pat[i])
}



#### Save structures ####
acceptance_output = result$post$accepts
k_output=result$full$k_s
lambda_output=result$full$lambda_s
psi_output = result$full$psi_s
eta_output = result$full$eta_s
theta_output = result$full$theta_s
delta_output=result$full$delta_s
tau_output = result$full$tau_s
phi_output = result$full$phi_s
sigma_phi_output = result$full$sigma_phi_s
sigma_eps_output = result$full$sigma_eps_s
sigmaj_output = result$full$sigma_j_s

##### Traceplots + acceptance rate ####
#Acceptance rate 
x11()
hist(acceptance_output, xlab = 'Rate', main='')

x11()
matplot(1:15000,k_ouput[],'l',main='k')

x11()
matplot(1:15000,k_output[(nburn+1):(niter+nburn)], 'l', ylab='', xlab='Iteration')


#Phi
x11()
matplot(phi_output[(nburn+1):(niter+nburn),,1],type='l',xlab = 'Iteration', main='', ylab='')
x11()
matplot(phi_output[(nburn+1):(niter+nburn),,12],type='l',xlab = 'Iteration', main='', ylab='')


#Sigma phi e sigma eps 
x11()
matplot(sigma_phi_output[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='', ylab='')
x11()
matplot(sigma_eps_output[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='', ylab='')


#Delta 
delta_out=sapply(delta_output, '[', 2)  #selezionare la componente desiderata 
x11()
matplot(delta_out[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='', ylab='')

#Tau
tau_out=sapply(tau_output, '[', 2)  #selezionare la componente desiderata 
x11()
matplot(tau_out[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='', ylab='')

# sigma_j 
x11()
sigmaj_out=sapply(sigmaj_output, '[', 3)
matplot(sigmaj_out[(nburn+1):(niter+nburn)], type='l',xlab = 'Iteration', 
        main='', ylab='')


#Theta (n patients X p)
#theta_out=sapply(theta_output, [, 1)  #(questo seleziona l' elemento 1 da tutte le matrici) 
#Per un paziente io devo selezionare tutta una riga--> devo avere p vettori 
paz=7  #Scegliere il paziente 
p=n_knots_m+4   #mettere il p usato 
ind=(paz-1)*p+1#primo elemento paziente paz
pos=8      #pos va da 0 a p-1 (0--> prima componente, p-1-->ultima componente) 
ind=ind+pos
theta_out=sapply(theta_output, '[',66) #theta del paziente 2, terzo p
x11()
matplot(theta_out[(nburn+1):(niter+nburn)],type='l',main='', ylab='', xlab='Iteration')

# p = 3
# n_pat = 2
# k = 1

#Lambda (p X k)
#65000, ognuna 40x17
#lambda_output[[2]]
#head(l_out)
#,3 significa terza riga, prima colonna della matrice
l_out=sapply(lambda_output, '[',3)   #selezionare la componente (il massimo è max(k)*p)
x11()
matplot(l_out[(nburn+1):(niter+nburn)],type='l', ylab='', main='', xlab='Iteration')

#Psi (p X k)
psi_out=sapply(psi_output, '[',3)   #selezionare la componente (il massimo è max(k)*p)
x11()
matplot(psi_out[(nburn+1):(niter+nburn)],type='l',main='', ylab='')

#Eta (n_patients X k)
#Per un paziente  devo selezionare tutta una riga--> devo avere k vettori 
paz=7  #Scegliere il paziente 
k_act=14  #mettere il k usato (cambia!!!) 
ind=(paz-1)*k_act+1#primo elemento paziente paz
pos=2        #pos va da 0 a k-1 (0--> prima componente, k-1-->ultima componente) 
ind=ind+pos
eta_out=sapply(eta_output, '[',2) #paziente 2, prima componente di k
x11()
matplot(eta_out[(nburn+1):(niter+nburn)],type='l',main='', ylab='')


#graphics.off()


#### Diagnostics  - create MCMC objects + ACF ####
sigma_phi_mcmc = mcmc(sigma_phi_output[-c(1:nburn)])
sigma_eps_mcmc = mcmc(sigma_eps_output[-c(1:nburn)])
lambda_mcmc = mcmc(l_out[-c(1:nburn)])
psi_mcmc = mcmc(psi_out[-c(1:nburn)])
tau_mcmc = mcmc(tau_out[-c(1:nburn)])
eta_mcmc = mcmc(eta_out[-c(1:nburn)])
delta_mcmc = mcmc(delta_out[-c(1:nburn)])
theta_mcmc = mcmc(theta_out[-c(1:nburn)])
sigmaj_mcmc = mcmc(sigmaj_out[-c(1:nburn)])

# acf
x11()
acf(sigma_phi_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(sigma_eps_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(lambda_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(psi_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(tau_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(eta_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(delta_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(theta_mcmc, lag.max=500, main = '', ylab = '')
x11()
acf(sigmaj_mcmc, lag.max=500, main = '', ylab = '')
