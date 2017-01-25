################# Changing intiialization of optim ############################
rp = 2; rq = 1;  #change these to test different series
num_initial = 25
maxp = 10; maxq = 10;
######### need to look at how to change initialization just set to random?
num_series = 1
vs = gen_series(num_series, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='init_test_method_series_21') # change file name!!!

results55 <- vector("list", length(vs))
for(i in 1:num_initial){
  ev = matrix(-Inf,maxp+1,maxq+1)
  print(i)
  x = vs[[1]]$series
  y = vs[[1]]$forc
  #run our algorithm
  for(p in 0:maxp){
    for(q in 0:maxq){
      initial = runif(p+q+1,-3,3)
      ev[p+1,q+1] = bayes_arima(x,y,p,q, initial)
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
saveRDS(results55,file='change_init_pm3_results_21') # change file name!!!