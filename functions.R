####### Packages #######

library(readxl)
library(bayesforecast)
library(forecast)
library(scoringRules)

####### Functions #######

loglik_arma_pq = function(y,mu0,phi,theta,sigma,yn,an,p,q){
  lk=vector(mode(0))
  if(p==1 && q==0) lk = -0.5*log(sigma^2)-0.5*log(2*pi)-1/(2*sigma^2)*(y-mu0-phi*yn)^2
  if(q==1 && p==0) lk = -0.5*log(sigma^2)-0.5*log(2*pi)-1/(2*sigma^2)*(y-mu0+theta*an)^2
  if(p==1 && q==1) lk = -0.5*log(sigma^2)-0.5*log(2*pi)-1/(2*sigma^2)*(y-mu0-phi*yn+theta*an)^2
  if(p==1 && q>1){
    for(j in 1:M){
      lk[j] = -0.5*log(sigma[j]^2)-0.5*log(2*pi)-1/(2*sigma[j]^2)*(y-mu0[j]-phi[j]*yn+theta[j,]%*%an[1:q])^2
    }
  }
  if(q==1 && p>1){
    for(j in 1:M){
      lk[j] = -0.5*log(sigma[j]^2)-0.5*log(2*pi)-1/(2*sigma[j]^2)*(y-mu0[j]-phi[j,]%*%yn[1:p]+theta[j]*an)^2
    }
  }
  if(p>1 && q==0){
    for(j in 1:M){
      lk[j] = -0.5*log(sigma[j]^2)-0.5*log(2*pi)-1/(2*sigma[j]^2)*(y-mu0[j]-phi[j,]%*%yn[1:p])^2
    }
  }
  if(q>1 && p==0){
    for(j in 1:M){
      lk[j] = -0.5*log(sigma[j]^2)-0.5*log(2*pi)-1/(2*sigma[j]^2)*(y-mu0[j]+theta[j,]%*%an[1:q])^2
    }
  }
  if(p>1 && q>1){
    for(j in 1:M){
      lk[j] = -0.5*log(sigma[j]^2)-0.5*log(2*pi)-1/(2*sigma[j]^2)*(y-mu0[j]-phi[j,]%*%yn[1:p]+theta[j,]%*%an[1:q])^2
    }
  }
  lk
}

marg.veros = function(logl.th0,d0){ # d0: numero de parametros
  
  l.bar0 = mean(logl.th0)
  s2.0 = var(logl.th0)
  
  lamb = 0.98
  l0.hat = max(l.bar0+s2.0,logl.th0)
  alpha.0 = d0/2
  logpy0 = l0.hat+alpha.0*log(1-lamb)
  return(logpy0)
}
