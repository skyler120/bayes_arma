rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
initialize = rep(0,p+q+1)
num_series = 25
vs = gen_series(num_series, rp,rq)
# saving vs, do for each rp and rq
saveRDS(vs,file='other_approaches_series_21') # change file name!!!

################# Maximum Likelihood Estimation #########################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  ev = matrix(-Inf,maxp+1,maxq+1)
  print(i)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = arima(x, order=c(p,0,q), include.mean=F, method="ML")
      ev[p+1,q+1] = a$loglik
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  
  results55[[i]] = c(mlp, mlq)
}
saveRDS(results55,file='MLE_approach_21') # change file name!!!



################# Cross Validation (LOOCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  ev = matrix(-Inf,maxp+1,maxq+1)
  print(i)
  x = vs[[i]]$series
  train = x[1:floor(length(x)*0.8)]
  test = x[floor(length(x)*0.8):length(x)]
  #y = vs[[i]]$forc
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = arima(train, order=c(p,0,q), include.mean=F, method="ML")
      fax = forecast(a, h=length(test))
      ev[p+1,q+1] = sqrt(sum((test-fax$mean)^2  ) / length(test))
    }
  }
  m = which(ev == min(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  
  results55[[i]] = c(mlp, mlq)
}
saveRDS(results55,file='LOOCV_approach_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=10, max.q = 10, allowmean = F, approximation = F, ic="aic", stepwise = T)
  aicp = a$arma[1]
  aicq = a$arma[2]
  
  results55[[i]] = c(aicp, aicq)
}
saveRDS(results55,file='aic_approach_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=10, max.q = 10, allowmean = F, approximation = F, ic="aicc", stepwise = T)
  aiccp = a$arma[1]
  aiccq = a$arma[2]
  
  results55[[i]] = c(aicp, aicq)
}
saveRDS(results55,file='aicc_approach_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=10, max.q = 10, allowmean = F, approximation = F, ic="bic", stepwise = T)
  bicp = a$arma[1]
  bicq = a$arma[2]
  
  results55[[i]] = c(aicp, aicq)
}
saveRDS(results55,file='aicc_approach_21') # change file name!!!
