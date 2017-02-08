################# Scaling the number of data points ############################
#This file creates nums number of series of size sampls and evaluates all methods
setwd("~/Desktop/bayes_arma")
source("all_code_get_results.R")
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;

#Change this to set the size of the time series
sampls = c(125, 625, 1875)
nums = 10

v_temp = gen_series(nums, sampls[length(sampls)], 1, rp,rq)
vs <- vector("list", nums)
for(i in 1:nums){
  v_samps <- vector("list", length(sampls))
  v_samps[[length(sampls)]] = v_temp[[i]]
  for(samp in 1:(length(sampls)-1)){
        n = sampls[samp]
        train = floor(0.8*n)
        v_samps[[samp]] = list(series = v_temp[[i]]$series[1:train], forc = v_temp[[i]]$forc[1:(n-train)], p = v_temp[[i]]$p, q = v_temp[[i]]$q)
  }
  vs[[i]] = v_samps
}
# saving vs, do for each rp and rq
saveRDS(vs,file='scaling2_sampsize_method_series_21') # change file name!!!

################# Bayes ARMA #############################
results55 <- vector("list", length(sampls))
for(i in 1:nums){
  print(sampls[i])
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    print(j)
    pt <- proc.time()[3]
    ev = matrix(-Inf,maxp+1,maxq+1)
    x = vs[[i]][[j]]$series
    y = vs[[i]][[j]]$forc
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
    results5[[j]] = c(proc.time()[3] - pt, bp, bq, fitted_acc(x,y,bp,bq,rp,rq), lag_terms)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='BARMA2_scaling_size_21') # change file name!!!

# saving results55, do for each rp and rq
#saveRDS(results55,file='scaling_sampsize_method_results_21') # change file name!!!

################# Maximum Likelihood Estimation #########################
results55 <- vector("list", length(sampls))
for(i in 1:nums){
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]]$series
    y = vs[[i]][[j]]$forc
    mlos = mle_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, mlos$ords, fitted_acc(x,y,mlos$ords[1],mlos$ords[2],rp,rq), mlos$coeffs)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='MLE2_scaling_size_21') # change file name!!!


################# Cross Validation (OFCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", length(sampls))
for(i in 1:nums){
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]]$series
    cvos = ofcv_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, cvos$ords, fitted_acc(x,y,cvos$ords[1],cvos$ords[2],rp,rq), cvos$coeffs)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='OFCV2_scaling_size_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", length(sampls))
for(i in 1:nums){
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    aicp = a$arma[1]
    aicq = a$arma[2]
    results5[[j]] = c(proc.time()[3] - pt, aicp, aicq, fitted_acc(x,y,aicp,aicq,rp,rq), a$sigma2, a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aic2_scaling_size_21') # change file name!!!

results55 <- vector("list", length(sampls))
for(i in 1:nums){
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aicc", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    aiccp = a$arma[1]
    aiccq = a$arma[2]
    results5[[j]] = c(proc.time()[3] - pt, aiccp, aiccq, fitted_acc(x,y,aiccp,aiccq,rp,rq), a$sigma2, a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aicc2_scaling_size_21') # change file name!!!

results55 <- vector("list", length(sampls))
for(i in 1:nums){
  results5 <- vector("list", length(nums))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="bic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    bicp = a$arma[1]
    bicq = a$arma[2]
    results5[[j]] = c(proc.time()[3] - pt, bicp, bicq, fitted_acc(x,y,bicp,bicq,rp,rq), a$sigma2, a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='bic2_scaling_size_21') # change file name!!!