##### Time series Simulation #####
######### rolling origin #########

set.seed(1234)

m = 100 # simulated series
n = 100 # number of observations
p1 = 1 # Model 0
q1 = 0
p2 = 0 # Model 1
q2 = 1
ytt = matrix(0,n,m)
f1=f2=r=pp=qq=PP=QQ=mar_lik_p=mar_lik_q=vector(mode(0))
pietraP=pietraQ=giniP=giniQ=log_likeli_q=log_likeli_p=vector(mode(0))

km = n-floor(2*n/3)+1
mm = km-1
PPp = QQq = matrix(0,mm,m)
logscorepp = logscoreqq = vector(mode(0))
logscoreP = logscoreQ = matrix(0,mm,m)

iteraciones = 15000
M=iteraciones/2
posterior_arma1 = array(data = NA, dim = c(iteraciones/2,p1+q1+3,n))
posterior_arma2 = array(data = NA, dim = c(iteraciones/2,p2+q2+3,n))

## Generated data with an ARMA(1,1)
p_s = 1
q_s = 1

for(j in 1:m){
ytt[,j] = arima.sim(n = n, list(ar = 0.2, ma = 0.5, sd = sqrt(0.16)))

## Generation of the time series
for(i in floor(2*n/3):(n-1)){
  yt = ytt[1:i,j]
  m1 = arima(yt,order=c(p1,0,q1))
  m2 = arima(yt,order=c(p2,0,q2))
  sfarma1 = stan_sarima(ts = yt,order = c(p1,0,q1),chains = 1,iter = iteraciones)
  sfarma2 = stan_sarima(ts = yt,order = c(p2,0,q2),chains = 1,iter = iteraciones)
  
  post_arma1 = extract_stan(sfarma1)
  posterior_arma1[,1,i] = post_arma1$mu0
  if(p1!=0) posterior_arma1[,2:(1+p1),i] = post_arma1$ar
  if(q1!=0) posterior_arma1[,(2+p1):(1+p1+q1),i] = post_arma1$ma
  posterior_arma1[,2+p1+q1,i] = post_arma1$sigma0
  posterior_arma1[,p1+q1+3,i] = post_arma1$loglik
  
  post_arma2 = extract_stan(sfarma2)
  posterior_arma2[,1,i] = post_arma2$mu0
  if(p2!=0) posterior_arma2[,2:(1+p2),i] = post_arma2$ar
  if(q2!=0) posterior_arma2[,(2+p2):(1+p2+q2),i] = post_arma2$ma
  posterior_arma2[,2+p2+q2,i] = post_arma2$sigma0
  posterior_arma2[,p2+q2+3,i] = post_arma2$loglik
  
  f1[i] = ytt[i+1,j]
  an1 = residuals(sfarma1)
  qqm = loglik_arma_pq(as.double(f1[i]),posterior_arma1[,1,i],posterior_arma1[,2:(1+p1),i],posterior_arma1[,(2+p1):(1+p1+q1),i],posterior_arma1[,2+p1+q1,i],yt[i:(i-p1+1)],an1[i:(i-q1+1)],p1,q1)
  qq[i] = mean(exp(qqm))
  logscoreqq[i] = mean(-qqm)
  mar_lik_q[i] = marg.veros(posterior_arma1[,p1+q1+3,i],p1+q1+2)
  
  an2 = residuals(sfarma2)
  ppm = loglik_arma_pq(as.double(f1[i]),posterior_arma2[,1,i],posterior_arma2[,2:(1+p2),i],posterior_arma2[,(2+p2):(1+p2+q2),i],posterior_arma2[,2+p2+q2,i],yt[i:(i-p2+1)],an2[i:(i-q2+1)],p2,q2)
  pp[i] = mean(exp(ppm))
  logscorepp[i] = mean(-ppm)
  mar_lik_p[i] = marg.veros(posterior_arma2[,p2+q2+3,i],p2+q2+2)
  r[i] = pp[i]/qq[i]
}

nn = n-1
pp = pp[floor(2*n/3):nn]
qq = qq[floor(2*n/3):nn]
r = r[floor(2*n/3):nn]
f1 = f1[floor(2*n/3):nn]
p_a = pp/sum(pp)
q_a = qq/sum(qq)
aux1 = cbind(r,f1,p_a,q_a)
aux1 = aux1[order(aux1[,3]), ]
aux2 = aux1[order(aux1[,4]), ]

logscoreP[,j] = logscorepp[floor(2*n/3):nn]
logscoreQ[,j] = logscoreqq[floor(2*n/3):nn]

for (k in 1:mm){
  PP[k] = sum(aux1[1:k,3])
  QQ[k] = sum(aux2[1:k,4])
}

PPp[,j] = PP/PP[mm]
QQq[,j] = QQ/QQ[mm]

# Pietra and Gini index
pietraP[j] = max(1:mm/mm-PPp[,j])
pietraQ[j] = max(1:mm/mm-QQq[,j])

ident = 1:mm/mm
giniP[j] = 1/2-1/2*sum((PPp[2:mm,j]+PPp[1:(mm-1),j])*(ident[2:mm]-ident[1:(mm-1)]))
giniQ[j] = 1/2-1/2*sum((QQq[2:mm,j]+QQq[1:(mm-1),j])*(ident[2:mm]-ident[1:(mm-1)]))

# Marginal likelihood
mar_lik_p = mar_lik_p[floor(2*n/3):nn]
mar_lik_q = mar_lik_q[floor(2*n/3):nn]

log_likeli_q[j] = median(mar_lik_q)
log_likeli_p[j] = median(mar_lik_p)
}

n_puntos = nrow(PPp)
n_replicas = ncol(QQq)
x = y = seq(0, 1, length.out = n_puntos)
mediaP = rowMeans(PPp)
mediaQ = rowMeans(QQq)

ic_infP = apply(PPp, 1, quantile, probs = 0.025)
ic_supP = apply(PPp, 1, quantile, probs = 0.975)
ic_infQ = apply(QQq, 1, quantile, probs = 0.025)
ic_supQ = apply(QQq, 1, quantile, probs = 0.975)

plot(x, mediaP, type = "l", lwd = 2, col = "blue",
     ylim = range(c(ic_infP, ic_supP)),
     xlab = "x", ylab = "")
lines(x,y,type = "l", col=1,lwd=3)
polygon(c(x, rev(x)), c(ic_supP, rev(ic_infP)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(x, mediaP, lwd = 2, col = "blue")

lines(x, mediaQ, type = "l", lwd = 2, col = "red",
      ylim = range(c(ic_infQ, ic_supQ)),
      xlab = "Población acumulada", ylab = "Variable acumulada")
polygon(c(x, rev(x)), c(ic_supQ, rev(ic_infQ)),
        col = rgb(1, 0, 1, 0.2), border = NA)
lines(x, mediaQ, lwd = 2, col = "red")
legend("bottomright",legend=c("M0: AR(1)","M1: MA(1)"),col=c(2,4),lty =c(1,1),lwd = rep(2,2),bty = "n")

plot(mediaQ,mediaP,type = "l", lwd = 3, col = "brown",
     ylim = range(c(0,1)),
     xlab = "x", ylab = "")
lines(x,y,type = "l", col=1,lwd=3,lty=1)
