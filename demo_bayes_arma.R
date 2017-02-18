########################################################################################
#                                                                                      #
#                                                                                      #
#                 Modify these lines to run the algorithm                              #
#                                                                                      #
setwd("~/bayes_arma") #set the working directory - default will be home directory      #
#                                                                                      #
#                                                                                      #
#                                                                                      #
#you may need to install the package forecast                                          #
source("all_code_get_results.R")                                                       #
#                                                                                      #
#                                                                                      #                            
#set this variable with the univariate time series you want to fit                     #
x =                                                                                    #
#set this variable with the univariate time series you want to forecast                #
y =                                                                                    #
#                                                                                      #
#                                                                                      #
maxp = 10; maxq = 10; #you may update these if you would like to search beyond 10      #
results_file_name = "demo_res"; #name of file (string) to save output                  #
#saved file is a list with elements:                                                   #
#mat - matrix of scores, orders, acc (RMSE residuals, RMSE forecast), coefficients     #
#                                                                                      #
#                                                                                      #
#                                                                                      #
#                                                                                      #
########################################################################################


#####################################Bayes ARMA Algorithm###############################
mx = mean(x); x = x - mx; #makes series zero mean                                      
ev = matrix(-Inf,maxp+1,maxq+1)
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
res_barma = list(mat = ev, orders = c(bp,bq), acc = rfitted_acc(x,y,bp,bq, mx), coeffs = lag_terms)
saveRDS(res_barma, results_file_name)

#step 1 initialization
  p = 2; q = 2; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 0; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 1; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  p = 0; q = 1; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  m = which(ev == max(ev), arr.ind = TRUE)
  curr = c(m[1,1]-1, m[1,2]-1)
  curr_val = ev[curr[1]+1, curr[2]+1]
  
  p=starting_vec[starting]; q=starting_vec[starting]; 
  ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))$pdm;
  curr = c(p,q)
  curr_val = ev[p+1,q+1]
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
      if(p>=0 && q>=0 && p<=10 && q<=10){
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
  res_step_barma = list(mat = ev, orders = curr, acc = rfitted_acc(x,y,curr[1],curr[2], mx))
  saverDS(res_step_barma, paste("step_", results_file_name, sep=""))