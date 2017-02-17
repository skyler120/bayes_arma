
#!/bin/bash
# Telling how many nodes and processors should be used.
#PBS -l nodes=1:ppn=8
# Naming the file
#PBS -N step-rt-1-2
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
R --vanilla > step-rt-1.out <<EOF

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



#################################################################################
#  All of the above 'R --vanilla...' is for the cluster
#  All of the below 'R --vanilla...' is an R file
#  This is the beginning of a 'regular' R file
#################################################################################

######## HERE BEGINS THE SIMULATION ########

################# Stepwise implementation of ABARMA #############################
#initialize
rp = 2; rq = 1;  #change these to test different series
maxp = 10; maxq = 10;
num_series = 25
samp_size = 125
nois = 1
starting = 3
starting_vec = c(1,5,9)
vs = gen_series(num_series,samp_size,nois, rp,rq)
#vs = readRDS('stepwise_series_86')
##################Stepwise bayesian ############################
results55 <- vector("list", length(vs))
for(i in 1:length(vs)){
  pt <- proc.time()[3]
  ev = matrix(-Inf,maxp+1,maxq+1)
  x = vs[[i]][["series"]]
  y = vs[[i]][["forc"]]
  
  #step 1 initialization
  p = 2; q = 2; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]];
  p = 0; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]];
  p = 1; q = 0; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]];
  p = 0; q = 1; ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]];
  m = which(ev == max(ev), arr.ind = TRUE)
  curr = c(m[1,1]-1, m[1,2]-1)
  curr_val = ev[curr[1]+1, curr[2]+1]
  # p=starting_vec[starting]; q=starting_vec[starting]; 
  # ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]];
  # curr = c(p,q)
  # curr_val = ev[p+1,q+1]
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
        if(ev[p+1,q+1] == -Inf) {ev[p+1,q+1] = bayes_arima(x,y,p,q,rep(0,p+q+1))[["pdm"]]}
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
  results55[[i]] = c(proc.time()[3] - pt, curr, fitted_acc(x,y,curr[1],curr[2],rp,rq))
}
saveRDS(results55,file='BARMA_stepwise_21') # change file name!!!


EOF
