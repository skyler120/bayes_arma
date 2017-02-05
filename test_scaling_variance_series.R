################# Scaling the variance of the arima process ############################
#This file creates nums number of series of size with varying amounts of noise in noiss

rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
noiss = 1:10
nums = 10
vs = gen_var_noise_series(nums, noiss, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='scaling_variance_method_series_21') # change file name!!!

################# Bayes ARMA #############################
results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    ev = matrix(-Inf,maxp+1,maxq+1)
    print(i)
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
    results5[[j]] = c(proc.time()[3] - pt, bp, bq, fitted_acc(x,y,bp,bq,rp,rq), lag_terms)
  }
  results55[[i]] = results5
}
# saving results55, do for each rp and rq
saveRDS(results55,file='ABARMA_scaling_variance_results_21') # change file name!!!


################# Maximum Likelihood Estimation #########################
results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]]$series
    y = vs[[i]]$forc
    mlos = mle_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, mlos$ords, fitted_acc(x,y,mlos$ords[1],mlos$ords[2],rp,rq), mlos$coeffs)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='MLE_scaling_var_21') # change file name!!!


################# Cross Validation (OFCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]]$series
    cvos = ofcv_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, cvos$ords, fitted_acc(x,y,cvos$ords[1],cvos$ords[2],rp,rq), cvos$coeffs)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='OFCV_scaling_var_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aic", stepwise = F)
    aicp = a$arma[1]
    aicq = a$arma[2]
    
    results5[[j]] = c(proc.time()[3] - pt, aicp, aicq, fitted_acc(x,y,aicp,aicq,rp,rq), a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aic_scaling_var_21') # change file name!!!

results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="aicc", stepwise = F)
    aiccp = a$arma[1]
    aiccq = a$arma[2]
    
    results5[[j]] = c(proc.time()[3] - pt, aiccp, aiccq, fitted_acc(x,y,aiccp,aiccq,rp,rq), a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aicc_scaling_var_21') # change file name!!!

results55 <- vector("list", length(noiss))
for(i in 1:length(noiss)){
  results5 <- vector("list", length(nums))
  for(j in 1:length(nums)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]]$series
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, approximation = F, ic="bic", stepwise = F)
    bicp = a$arma[1]
    bicq = a$arma[2]
    
    results5[[j]] = c(proc.time()[3] - ptbicp, bicq, fitted_acc(x,y,bicp,bicq,rp,rq), a$coef)
  }
  results55[[i]] = results5
}
saveRDS(results55,file='bic_scaling_var_21') # change file name!!!