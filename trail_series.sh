
#!/bin/bash
# Telling how many nodes and processors should be used.
#PBS -l nodes=1:ppn=8
# Naming the file
#PBS -N scalingseries-rt-1-2
# Outputting error
#PBS -j oe
# Not sure what the two next lines do
#PBS -q default
#PBS -S /bin/bash
#PBS -m abe
#PBS -M ss3349@cornell.edu

######## I WONDER IF I NEED THIS LINE ######
cd $PBS_O_WORKDIR

# Telling cluster that you are using R
R --vanilla > ss-rt-2.out <<EOF

# Looking for what machines are available to use.
setwd("/home/fs01/ss3349/bayes_arma")
source("all_code_get_results.R")
#intall.packages("forecast")
library(forecast)
#install.packages("cubature")
library(cubature)
#install.packages("GenSA")
library(GenSA)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("numDeriv")
library(numDeriv)

#library(snowfall)
# library(snow)
#pbsnodefile = Sys.getenv("PBS_NODEFILE")
#machines <- scan(pbsnodefile, what="")
#print(machines)
#nmach = length(machines)
#nmach

# Initializing the nodes
#sfInit(parallel=TRUE,type='SOCK',cpus=nmach,socketHosts=machines)


#################################################################################
#  All of the above 'R --vanilla...' is for the cluster
#  All of the below 'R --vanilla...' is an R file
#  This is the beginning of a 'regular' R file
#################################################################################
#sfSource('lib.R')

######## HERE BEGINS THE SIMULATION ########

#setwd("/home/fs01/ss3349/bayes_arma")
#source("all_code_get_results.R")
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;

#Change this to set the size of the time series
sampls = c(125, 250, 625, 1250, 6250, 12500)
nums = 1
sampls = c(12500)
#v_temp = gen_series(nums, sampls[length(sampls)], 1, rp,rq)
v_temp = readRDS('scaline13_full_sampsize_method_series_21')
vs <- vector("list", nums)
for(i in 1:nums){
  v_samps <- vector("list", length(sampls))
  v_samps[[length(sampls)]] = v_temp[[i]]
  #for(samp in 1:(length(sampls)-1)){
  #      n = sampls[samp]
  #      train = floor(0.8*n)
  #      v_samps[[samp]] = list(series = v_temp[[i]][["series"]][1:train], forc = v_temp[[i]][["forc"]][1:(n-train)], p = v_temp[[i]][["p"]], q = v_temp[[i]][["q"]])
  #}
  vs[[i]] = v_samps
}
# saving vs, do for each rp and rq
saveRDS(vs,file='scaling14_sampsize_method_series_21') # change file name!!!
#sampls = c(125,250,625,1250,6250)


################# Bayes ARMA #############################
results55 <- vector("list", nums)
for(i in 1:nums){
  print(i)
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    print(sampls[j])
    pt <- proc.time()[3]
    ev = matrix(-Inf,maxp+1,maxq+1)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]]
    #run our algorithm
    prev_pdm = -Inf
    best_params = NULL
    for(p in 0:maxp){
      for(q in 0:maxq){
        initialize = rep(0,p+q+1)
        bayes_ar_res = bayes_arima(x,y,p,q, initialize)
        ev[p+1,q+1] = bayes_ar_res[["pdm"]]
        if(ev[p+1, q+1]> prev_pdm){
          prev_pdm = ev[p+1, q+1]
          best_params = bayes_ar_res[["transformed_params"]]
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
saveRDS(results55,file='BARMA14_scaling_size_21') # change file name!!!

# saving results55, do for each rp and rq
#saveRDS(results55,file='scaling_sampsize_method_results_21') # change file name!!!

################# Maximum Likelihood Estimation #########################
results55 <- vector("list", nums)
for(i in 1:nums){
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]]
    mlos = mle_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, mlos[["ords"]], fitted_acc(x,y,mlos[["ords"]][1],mlos[["ords"]][2],rp,rq), mlos[["coeffs"]])
  }
  results55[[i]] = results5
}
saveRDS(results55,file='MLE14_scaling_size_21') # change file name!!!


################# Cross Validation (OFCV) #########################
#based on forecast accuracy iteratively and find the lowest forecast on a held out set 
results55 <- vector("list", nums)
for(i in 1:nums){
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]] 
    cvos = ofcv_arma(x)
    results5[[j]] = c(proc.time()[3] - pt, cvos[["ords"]], fitted_acc(x,y,cvos[["ords"]][1],cvos[["ords"]][2],rp,rq), cvos[["coeffs"]])
  }
  results55[[i]] = results5
}
saveRDS(results55,file='OFCV14_scaling_size_21') # change file name!!!


################# IC based methods #############################
results55 <- vector("list", nums)
for(i in 1:nums){
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]]
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    aicp = a[["arma"]][1]
    aicq = a[["arma"]][2]
    results5[[j]] = c(proc.time()[3] - pt, aicp, aicq, fitted_acc(x,y,aicp,aicq,rp,rq), a[["sigma2"]], a[["coef"]])
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aic14_scaling_size_21') # change file name!!!

results55 <- vector("list", nums)
for(i in 1:nums){
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]]
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aicc", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    aiccp = a[["arma"]][1]
    aiccq = a[["arma"]][2]
    results5[[j]] = c(proc.time()[3] - pt, aiccp, aiccq, fitted_acc(x,y,aiccp,aiccq,rp,rq), a[["sigma2"]], a[["coef"]])
  }
  results55[[i]] = results5
}
saveRDS(results55,file='aicc14_scaling_size_21') # change file name!!!

results55 <- vector("list", nums)
for(i in 1:nums){
  results5 <- vector("list", length(sampls))
  for(j in 1:length(sampls)){
    pt <- proc.time()[3]
    print(i)
    x = vs[[i]][[j]][["series"]]
    y = vs[[i]][[j]][["forc"]]
    a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="bic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
    bicp = a[["arma"]][1]
    bicq = a[["arma"]][2]
    results5[[j]] = c(proc.time()[3] - pt, bicp, bicq, fitted_acc(x,y,bicp,bicq,rp,rq), a[["sigma2"]], a[["coef"]])
  }
  results55[[i]] = results5
}
saveRDS(results55,file='bic14_scaling_size_21') # change file name!!!


EOF
