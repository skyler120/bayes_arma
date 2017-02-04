################# Scaling the number of data points ############################
#This file creates nums number of series of size sampls and evaluates all methods
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
initialize = rep(0,p+q+1)

#Change this to set the size of the time series
sampls = c(10,seq(100,1000,10), seq(1000,10000,1000), seq(10000,50000, 20000))
nums = 10

vs = gen_var_samp_series(nums, sampls, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='scaling_sampsize_method_series_21') # change file name!!!

################# Bayes ARMA #############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  ev = matrix(-Inf,maxp+1,maxq+1)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  #run our algorithm
  prev_pdm = -Inf
  best_params = NULL
  for(p in 0:maxp){
    for(q in 0:maxq){
      bayes_ar_res = bayes_arima(x,y,p,q, initialize)
      ev[p+1,q+1] = bayes_ar_res$pdm
      if(ev[p+1, q+1]> prev_pdm){
        prev_pdm = ev[p+1, q+1]
        best_params = bayes_ar_res$transformed_params
      }
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  bp = m[1,1] - 1
  bq = m[1,2] - 1
  lag_terms = get_coeffs(best_params, bp, bq)
  results55[[i]] = c(bp, bq, fitted_acc(x,y,bp,bq,rp,rq))
}
saveRDS(results55,file='BARMA_scaling_size_21') # change file name!!!

# saving results55, do for each rp and rq
saveRDS(results55,file='scaling_sampsize_method_results_21') # change file name!!!

################# Maximum Likelihood Estimation #########################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  mlos = mle_arma(x)
  results55[[i]] = c(mlos[1], mlos[2], fitted_acc(x,y,mlos[1],mlos[2],rp,rq))
}
saveRDS(results55,file='MLE_scaling_size_21') # change file name!!!


################# Cross Validation (OFCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  cvos = ofcv_arma(x)
  results55[[i]] = c(cvos[1], cvos[2], fitted_acc(x,y,cvos[1],cvos[2],rp,rq))
}
saveRDS(results55,file='OFCV_scaling_size_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aic", stepwise = F)
  aicp = a$arma[1]
  aicq = a$arma[2]
  
  results55[[i]] = c(aicp, aicq, fitted_acc(x,y,aicp,aicq,rp,rq))
}
saveRDS(results55,file='aic_scaling_size_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aicc", stepwise = F)
  aiccp = a$arma[1]
  aiccq = a$arma[2]
  
  results55[[i]] = c(aiccp, aiccq, fitted_acc(x,y,aiccp,aiccq,rp,rq))
}
saveRDS(results55,file='aicc_scaling_size_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  print(i)
  x = vs[[i]]$series
  a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="bic", stepwise = F)
  bicp = a$arma[1]
  bicq = a$arma[2]
  
  results55[[i]] = c(bicp, bicq, fitted_acc(x,y,bicp,bicq,rp,rq))
}
saveRDS(results55,file='aicc_scaling_size_21') # change file name!!!