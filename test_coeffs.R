setwd("~/Desktop/bayes_arma")
source("all_code_get_results.R")
rp = 1; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
num_series = 10
samp_size = 125
nois = 1
vs = gen_series(num_series,samp_size,nois, rp,rq)
# saving vs, do for each rp and rq
#vs = readRDS("init_test_method_series_21")
saveRDS(vs,file='test_coeffs_11') # change file name!!!

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
saveRDS(results55,file='BARMA_test_coeff_21')