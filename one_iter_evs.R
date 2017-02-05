################# Basic comparison of algorithm for each of different approaches #############################
#This file creates 25 series of orders rp and rq and outputs 
#the orders, training rmse, and forecast rmse for each method
setwd("~/Desktop/bayes_arma")
source("all_code_get_results.R")
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
num_series = 1
samp_size = 125
nois = 1
vs = gen_series(num_series,samp_size,nois, rp,rq)
# saving vs, do for each rp and rq
saveRDS(vs,file='ev_other_approaches_series_21') # change file name!!!


################# Bayes ARMA #############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  ev = matrix(-Inf,maxp+1,maxq+1)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  #run our algorithm
  prev_pdm = -Inf
  best_params = NULL
  for(p in 0:maxp){
    for(q in 0:maxq){
      initialize = rep(0,p+q+1)
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
  results55[[i]] = list(mat = ev, res = c(proc.time()[3] - pt, bp, bq, fitted_acc(x,y,bp,bq,rp,rq), lag_terms))
}
saveRDS(results55,file='ev_BARMA_approach_21') # change file name!!!

################# Maximum Likelihood Estimation #########################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  mlos = mle_arma(x)
  results55[[i]] = list(mat = mlos$mat, res = c(proc.time()[3] - pt, mlos$ords, fitted_acc(x,y,mlos$ords[1],mlos$ords[2],rp,rq), mlos$coeffs))
}
saveRDS(results55,file='ev_MLE_approach_21') # change file name!!!


################# Cross Validation (OFCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  x = vs[[i]]$series
  cvos = ofcv_arma(x)
  results55[[i]] = list(mat = cvos$mat, res = c(proc.time()[3] - pt, cvos$ords, fitted_acc(x,y,cvos$ords[1],cvos$ords[2],rp,rq), cvos$coeffs))
}
saveRDS(results55,file='ev_OFCV_approach_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  x = vs[[i]]$series
  aicos = aic_arma(x)
  #a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aic", stepwise = F)
  #aicp = a$arma[1]
  #aicq = a$arma[2]
  
  results55[[i]] = list(mat = aicos$mat, res = c(proc.time()[3] - pt, aicos$ords[1], aicos$ords[2], fitted_acc(x,y,aicos$ords[1],aicos$ords[2],rp,rq), aicos$coeffs))
}
saveRDS(results55,file='ev_aic_approach_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  x = vs[[i]]$series
  aiccos = aicc_arma(x)
  #a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aicc", stepwise = F)
  #aiccp = a$arma[1]
  #aiccq = a$arma[2]
  
  results55[[i]] = list(mat = aiccos$mat, res = c(proc.time()[3] - pt, aiccos$ords[1], aiccos$ords[2], fitted_acc(x,y,aiccos$ords[1],aiccos$ords[2],rp,rq), aiccos$coeffs))
}
saveRDS(results55,file='ev_aicc_approach_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  print(i)
  x = vs[[i]]$series
  bicos = bic_arma(x)
  #a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="bic", stepwise = F)
  #bicp = a$arma[1]
  #bicq = a$arma[2]
  
  results55[[i]] = list(mat = bicos$mat, res = c(proc.time()[3] - pt, bicos$ords[1], bicos$ords[2], fitted_acc(x,y,bicos$ords[1],bicos$ords[2],rp,rq), bicos$coeffs))
}
saveRDS(results55,file='ev_bic_approach_21') # change file name!!!

