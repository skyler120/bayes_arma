################# Scaling the variance of the arima process ############################
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
initialize = rep(0,p+q+1)
noiss = c(1:10, seq(10,100,10))
vs = gen_var_noise_series(noiss, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='scaling_variance_method_series_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  ev = matrix(-Inf,maxp+1,maxq+1)
  print(i)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  #run our algorithm
  for(p in 0:maxp){
    for(q in 0:maxq){
      ev[p+1,q+1] = bayes_arima(x,y,p,q, initialize)
      if(is.nan(ev[p+1, q+1])){
        ev[p+1, q+1] = -Inf
      }
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  bp = m[1,1] - 1
  bq = m[1,2] - 1
  
  arma_fits = fitted_acc(x,y,bp, bq, rp,rq)
  ap = arma_fits[1]
  aq = arma_fits[2]
  results55[[i]] = c(rp-ap, rp-bp, rq-aq, rq-bq, rp + rq - (ap + aq), rp + rq - (bp+bq), arma_fits[3:length(arma_fits)])
}

# saving results55, do for each rp and rq
saveRDS(results55,file='scaling_variance_method_results_21') # change file name!!!