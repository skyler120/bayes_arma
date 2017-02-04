################# Changing intiialization of optim ############################
rp = 2; rq = 1;  #change these to test different series
num_initial = 25
maxp = 10; maxq = 10;
######### need to look at how to change initialization just set to random?
num_series = 25
vs = gen_series(num_series, rp,rq)

# saving vs, do for each rp and rq
saveRDS(vs,file='init_test_method_series_21') # change file name!!!

results55 <- vector("list", length(vs))
for(j in 1:num_series){
  results_os <- vector("list", num_initial)
  for(i in 1:num_initial){
    ev = matrix(-Inf,maxp+1,maxq+1)
    print(i)
    x = vs[[1]]$series
    y = vs[[1]]$forc
    #run our algorithm  
    prev_pdm = -Inf
    best_params = NULL
    for(p in 0:maxp){
      for(q in 0:maxq){
        initial = runif(p+q+1,-3,3)
        bayes_ar_res = bayes_arima(x,y,p,q, initial)
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
    results_os[[i]] = list(orders = c(bp,bq), lags = lag_terms)
  }
  results55[[j]] = results_os
}

# saving results55, do for each rp and rq
saveRDS(results55,file='change_init_pm3_results_21') # change file name!!!