################# Scaling the number of data points ############################
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
initialize = rep(0,p+q+1)
sampls = c(10,seq(1,1000,10), seq(1000,10000,1000), seq(10000,50000, 20000))
vs = gen_var_samp_series(sampls, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='scaling_sampsize_method_series_21') # change file name!!!

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
saveRDS(results55,file='scaling_sampsize_method_results_21') # change file name!!!