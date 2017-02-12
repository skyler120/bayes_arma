################# Stepwise implementation of ABARMA #############################
setwd("~/Desktop/Cornell/ORIE 6741/bayes_arma")
#setwd("~/Desktop/bayes_arma")
source("all_code_get_results.R")
#initialize
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
num_series = 1
samp_size = 125
nois = 1
vs = gen_series(num_series,samp_size,nois, rp,rq)
#vs = readRDS('other_approaches_series_21')
##################Stepwise bayesian ############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  ev = matrix(-Inf,maxp+1,maxq+1)
  x = vs[[i]]$series
  y = vs[[i]]$forc
  
  #step 1 initialization
  p = 2; q = 2; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 0; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 1; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 0; q = 1; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  m = which(ev == max(ev), arr.ind = TRUE)
  curr = c(m[1,1]-1, m[1,2]-1)
  curr_val = ev[curr[1]+1, curr[2]+2]
  
  #step 2
  count = 0
  params_list = numeric(16)
  while(count < (length(params_list)/2)){
    print(curr)
    params_list = c(curr[1]+1, curr[2], curr[1]-1, curr[2], curr[1], curr[2]+1, curr[1], curr[2]-1,
                    c(curr[1]+1, curr[2]+1, curr[1]-1, curr[2]+1, curr[1]+1, curr[2]-1, curr[1]-1, curr[2]-1))
    count = 0;
    for(j in 1:(length(params_list)/2)){
      p = params_list[(2*j-1)]; q = params_list[(2*j)];
      if(p>0 && q>0){
        if(ev[p+1,q+1] == -Inf) {ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm}
        if(ev[p+1, q+1]> curr_val){
          curr = c(p,q)
          curr_val = ev[p+1,q+1]
          break
        }
        else{count = count+1}
      }
      else{count = count+1}
    }
  }
  results55[[i]] = c(proc.time()[3] - pt, curr)
}
saveRDS(results55,file='BARMA_stepwise_21') # change file name!!!

