############### Requirements ###############
require("forecast")
require("numDeriv")
require("GenSA")


############################# Recursive Function for phi and pi ##############################
phiphi = function(k,r){
  phimat = matrix(0,k,k)
  for(j in 1:k){
    for(i in 1:j){
      if(i==j){
        phimat[j,i] =tanh(r[j])
      }
      else{
        phiipre = phimat[i,j-1]
        phiki = phimat[j-i,j-i]
        phimat[j,i] = phiipre - tanh(r[j])*phiki
      }
    }
  }
  return(phimat[k,])
}


############################# Recursive Function for mu_t ################################
muts <- function(phir,pis){
  poq = max(length(phir), length(pis))
  mus = rep(0,poq)
  xapp = c(rep(0,poq), x)
  for(t in (poq+1):(poq+length(x))){
    if(p==0){
      p1 = 0
    }
    else{
      p1 = sum(phir*rev(xapp[(t-p):(t-1)]))
    }
    if(q==0){
      p2 = 0
    }
    else{
      p2 = sum(pis*(rev(xapp[(t-q):(t-1)]) - mus[(t-q):(t-1)]))
    }
    mus = c(mus, p1-p2)
  }
  return(mus)
}


#############################  Function for logp ###############################################

logp <- function(params){
  n = length(x)
  poq = max(p,q)
  phir = 0
  pis = 0
  tsigsq = params[1]
  if(p>0){
    phir = phiphi(p,params[2:(1+p)])
    if(q>0){
      pis = phiphi(q,params[(2+p):length(params)])
    }
  }else if(q>0){
    pis = phiphi(q,params[(2+p):length(params)])
  }else{
    phir = 0
    pis = 0
  }
  mut = muts(phir,pis)
  val = sum((x-mut[(poq+1):length(mut)])^2)
  return(-1*(n+2)/2*tsigsq - 1/2*exp(-tsigsq)*val)
}


############### This is where the Hessians are computed ###############

# Hessian for log p(D,theta | M)
logpTheta.H = function(params){
  return(hessian(logp,params))
}

# diagonal of Hessian for log(J_r)
logJr.H = function(params){
  if (p>0){
    diagH = c(0, 2*tanh(params[2:(1+p)])^2-2, numeric(q))
  }else{
    diagH = numeric(1+p+q)
  }
  return(diagH)
}

# diagonal of Hessian for log(J_s)
logJs.H = function(params){
  if (q>0){
    diagH = c(0, numeric(p), 2*tanh(params[(2+p):(1+p+q)])^2-2)
  }else{
    diagH = numeric(1+p+q)
  }
  return(diagH)
}

# diagonal of Hessian for log(J_phir)
logJphir.H = function(params){
  if (p>0){
    diagH = c(0, (0:(1-p))*(1-tanh(params[2:(1+p)])^2), numeric(q))
  }else{
    diagH = numeric(1+p+q)
  }
  return(diagH)
}

# diagonal of Hessian for log(J_pis)
logJpis.H = function(params){
  if (q>0){
    diagH = c(0, numeric(p), (0:(1-q))*(1-tanh(params[(2+p):(1+p+q)])^2))
  }else{
    diagH = numeric(1+p+q)
  }
  return(diagH)
}



############### This is where we find the argmax of the integrand ###############

# integrand function, note optim function finds argmin i.e. of -f
f = function(params){
  # p(D,theta*|M), off by a constant 
  pfunc = exp(logp(params))
  # J_r
  if (p>0){
    Jr = prod(1-tanh(params[2:(1+p)])^2)
  }else{
    Jr = 1
  }
  # J_s
  if (q>0){
    Js = prod(1-tanh(params[(2+p):(1+p+q)])^2)
  }else{
    Js = 1
  }
  # J_phir
  if (p>=2){
    Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2), 1-tanh(params[1+2*(1:floor(p/2))]))
  }else if (p==1){
    Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2))
  }else{
    Jphir = 1
  }
  # J_pis
  if (q>=2){
    Jpis = prod((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2), 1-tanh(params[1+p+2*(1:floor(q/2))]))
  }else if (q==1){
    Jpis = prod((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2))
  }else{
    Jpis = 1
  }
  # J_sigma
  Jsigma = exp(params[1])
  return(pfunc*Jr*Js*Jphir*Jpis*Jsigma)
}


############### This is where we find the argmax of the log integrand ###############

# negative of log integrand function, note optim function finds argmin
logf_neg = function(params){
  # log p(D,theta*|M), off by a constant 
  pfunc = logp(params)
  # J_r
  if (p>0){
    Jr = sum(log(1-tanh(params[2:(1+p)])^2))
  }else{
    Jr = 1
  }
  # J_s
  if (q>0){
    Js = sum(log(1-tanh(params[(2+p):(1+p+q)])^2))
  }else{
    Js = 1
  }
  # J_phir
  if (p>=2){
    Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2), 1-tanh(params[1+2*(1:floor(p/2))]))
  }else if (p==1){
    Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2))
  }else{
    Jphir = 1
  }
  # J_pis
  if (q>=2){
    Jpis = sum(log((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2)), log(1-tanh(params[1+p+2*(1:floor(q/2))])))
  }else if (q==1){
    Jpis = sum(log((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2)))
  }else{
    Jpis = 1
  }
  # J_sigma
  Jsigma = params[1]
  return(-pfunc-Jr-Js-Jphir-Jpis-Jsigma)
}


# SIMULATION 1

x = arima.sim(n = 100, list(ar = c(0.8897, -0.4858), ma = c(0.4)))
ev2 = matrix(-1000000,5+1,5+1)
for(p in 0:5){
  for(q in 0:5){
    if(p+q>0){
      print(c(p,q))
      init = rnorm(p+q+1)
      res = GenSA(init, logf_neg, rep(-10,1+p+q), rep(10,1+p+q))
      params = res$par
      hess = logpTheta.H(params) + diag(logJr.H(params) + logJs.H(params) + logJphir.H(params) + logJpis.H(params))
      ev2[p+1,q+1] = f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2)
      #ev1 = c(ev1, f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2))
    }
  }
}

set.seed(150)
x1 = arima.sim(n = 100, list(ar = c(0.8897, -0.4858, 0.6, -0.2, 0.1), ma = c(-0.2279, 0.2488, 0.4)))
set.seed(150)
x2 = arima.sim(n = 100, list(ar = c(0.8897, -0.4858), ma = c(0.4)))
set.seed(150)
x3 = arima.sim(n = 100, list(ar = c(0.8897, -0.4858, 0.6, -0.2, 0.1, -0.3, 0.2, -0.2), ma = c(0.1, -0.2, 0.1, -0.2279, 0.2488, 0.4)))

x = x2
for(p in 0:5){
  for(q in 0:5){
    if(p+q>0){
      print(c(p,q))
      init = rnorm(p+q+1)
      res = GenSA(init, logf_neg, rep(-10,1+p+q), rep(10,1+p+q))
      params = res$par
      hess = logpTheta.H(params) + diag(logJr.H(params) + logJs.H(params) + logJphir.H(params) + logJpis.H(params))
      ev3[p+1,q+1] = f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2)
      #ev1 = c(ev1, f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2))
    }
  }
}


x = x3
ev4 = matrix(-1000000,10+1,10+1)
for(p in 0:10){
  for(q in 0:10){
    if(p+q>0){
      print(c(p,q))
      init = rnorm(p+q+1)
      res = GenSA(init, logf_neg, rep(-10,1+p+q), rep(10,1+p+q))
      params = res$par
      hess = logpTheta.H(params) + diag(logJr.H(params) + logJs.H(params) + logJphir.H(params) + logJpis.H(params))
      ev4[p+1,q+1] = f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2)
      #ev1 = c(ev1, f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2))
    }
  }
}
