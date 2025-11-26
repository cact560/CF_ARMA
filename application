##### Application #####

data_gas <- read_excel("data_gas.xlsx")

data_gas = ts(data_gas$`United States Natural Gas Industrial Price (Dollars per Thousand Cubic Feet)`,start=c(2010,1),frequency=12)
plot(data_gas,ylab="")

auto.sarima(data_gas) # AR(1)

## Concentration function
p1=1 # model proposed by auto.sarima (an AR(1)),
q1=0
p2=2 # other model
q2=0

n = length(data_gas)
m = n-floor(2*n/3)+1

yt = data_gas

f1=f1p=f2p=r=pp=qq=PP=QQ=mar_lik_p=mar_lik_q=vector(mode(0))
pietraP=pietraQ=giniP=giniQ=log_likeli_q=log_likeli_p=vector(mode(0))
logscorepp = logscoreqq = vector(mode(0))
crps_m1 <- crps_m2 <- numeric(n)

iteraciones = 15000
M=iteraciones/2
posterior_arma1 = array(data = NA, dim = c(iteraciones/2,p1+q1+3,n))
posterior_arma2 = array(data = NA, dim = c(iteraciones/2,p2+q2+3,n))

## Generation of the series
for(i in floor(2*n/3):(n-1)){
  yyt=yt[1:i]
  m1=arima(yyt,order=c(p1,0,q1))
  m2=arima(yyt,order=c(p2,0,q2))
  sfarma1 = stan_sarima(ts = yyt,order = c(p1,0,q1),chains = 1,iter = iteraciones)
  sfarma2 = stan_sarima(ts = yyt,order = c(p2,0,q2),chains = 1,iter = iteraciones)

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
  f1p[i] = forecast(sfarma1,h = 1)$mean
  f2p[i] = forecast(sfarma2,h = 1)$mean
  f1[i]=yt[i+1]

  an1 = residuals(sfarma1)
  qqm = loglik_arma_pq(as.double(f1[i]),posterior_arma1[,1,i],posterior_arma1[,2:(1+p1),i],posterior_arma1[,(2+p1):(1+p1+q1),i],posterior_arma1[,2+p1+q1,i],yyt[i:(i-p1+1)],an1[i:(i-q1+1)],p1,q1)
  qq[i]= mean(exp(qqm))
  logscoreqq[i]=mean(-qqm)
  crps_m1[i] <- crps_sample(f1[i],f1p[i])
  mar_lik_q[i] = marg.veros(posterior_arma1[,p1+q1+3,i],p1+q1+2)

  an2 = residuals(sfarma2)
  ppm = loglik_arma_pq(as.double(f1[i]),posterior_arma2[,1,i],posterior_arma2[,2:(1+p2),i],posterior_arma2[,(2+p2):(1+p2+q2),i],posterior_arma2[,2+p2+q2,i],yyt[i:(i-p2+1)],an2[i:(i-q2+1)],p2,q2)
  pp[i]= mean(exp(ppm))
  logscorepp[i]=mean(-ppm)
  crps_m2[i] <- crps_sample(f1[i],f2p[i])
  mar_lik_p[i] = marg.veros(posterior_arma2[,p2+q2+3,i],p2+q2+2)
  r[i] = pp[i]/qq[i]
}

nn=n-1
pp = pp[floor(2*n/3):nn]
qq = qq[floor(2*n/3):nn]
r = r[floor(2*n/3):nn]
f1 = f1[floor(2*n/3):nn]
p_a = pp/sum(pp)
q_a = qq/sum(qq)
aux1 = cbind(r,f1,p_a,q_a)
aux1 <- aux1[order(aux1[,3]), ]
aux2 <- aux1[order(aux1[,4]), ]

mm=m-1
for (j in 1:mm){
  PP[j] =sum(aux1[1:j,3])
  QQ[j] =sum(aux2[1:j,4])
}

PPp = PP/PP[mm]
QQq = QQ/QQ[mm]

par(mfrow=c(1,1))
x=y = seq(0,1,length = 100)
plot(QQq,PPp,type = "l", lwd = 3, col = "brown",ylim = range(c(0,1)),xlab = "x", ylab = "",main = "(a)")
lines(x,y,type = "l", col=1,lwd=3,lty=1)

# Pietra and Gini index
pietraPQ = max(QQq-PPp)
giniPQ = 1-sum((PPp[2:mm]+PPp[1:(mm-1)])*(QQq[2:mm]-QQq[1:(mm-1)]))

# logscore
logscoreP = mean(logscorepp[floor(2*n/3):nn])
logscoreQ = mean(logscoreqq[floor(2*n/3):nn])

# --- Average CRPS over forecast period ---
mean_crps_m1 <- mean(crps_m1[floor(2 * n / 3):(n - 1)], na.rm = TRUE)
mean_crps_m2 <- mean(crps_m2[floor(2 * n / 3):(n - 1)], na.rm = TRUE)
