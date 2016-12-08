############################# Recursive Function for phi and pi ###############################################

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


############################# Recursive Function for mu_t ###############################################

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

############### Compute Accuracy ###############

fitted_acc <- function(x, y, bp, bq, rp, rq){
  arima_x1 = auto.arima(x, d=0, max.p=5, max.q = 5, allowmean = F, approximation = F)
  barima_x1 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
  tarima_x1 = arima(x, order=c(rp,0,rq), include.mean=F, method="ML")
  a = max(arima_x1$residuals)
  b = max(barima_x1$residuals)
  
  if(a>b){
    plot(1:length(x),x - arima_x1$residuals, col="blue", type='l')
    lines(1:length(x), x , col="black")
    lines(1:length(x),x - barima_x1$residuals, col="red")
    lines(1:length(x),x - tarima_x1$residuals, col="green")
    
  }else{
    plot(1:length(x),x - barima_x1$residuals, col="red", type='l')
    lines(1:length(x), x , col="black")
    lines(1:length(x),x - arima_x1$residuals, col="blue")
    lines(1:length(x),x - tarima_x1$residuals, col="green")
  }
  
  rmse_auto_x  = sqrt(sum((arima_x1$residuals)^2  ) / length(x))
  rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))
  rmse_true_x  = sqrt(sum((tarima_x1$residuals)^2  ) / length(x))
  
  fax = forecast(arima_x1, h=length(y))
  fbx = forecast(barima_x1, h=length(y))
  ftx = forecast(tarima_x1, h=length(y))
  
  trmse_auto_x  = sqrt(sum((y-fax$mean)^2  ) / length(x))
  trmse_bayes_x  = sqrt(sum((y-fbx$mean)^2  ) / length(x))
  trmse_true_x  = sqrt(sum((y-ftx$mean)^2  ) / length(x))
  
  
  return(c(arima_x1$arma[1], arima_x1$arma[2], rmse_auto_x, rmse_bayes_x, rmse_true_x, trmse_auto_x, trmse_bayes_x, trmse_true_x))
}

############### Generate Sample Data ###############

gen_arm <- function(pp,qq){
  if((pp+qq)>0){
    phis = runif(pp,-1,1)
    pis = runif(qq,-1,1)
    x = arima.sim(n = 125, list(ar = phis, ma = pis), sd = runif(1,0,1))
    return(list(series=x[1:100],forc = x[101:125], p = pp, q = qq))
  }
}

gen_series <- function(pp,qq){
  res <- vector("list", 3000)
  for(i in 1:3000){
    res[[i]] = try(gen_arm(pp,qq))
    if (typeof(res[[i]])=="character"){
      res[[i]]  = NULL
    }
  }
  valid_series = Filter(Negate(is.null), res)
  return(valid_series)
}

############### Bayesian ARMA Model Order Determination ###############

bayes_arima <- function(x,y, p, q){
  if(p+q>0){
    init = rnorm(p+q+1)
    res = optim(init, logf_neg)
    params = res$par
    hess = logpTheta.H(params) + diag(logJr.H(params) + logJs.H(params) + logJphir.H(params) + logJpis.H(params))
    return(f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2))

  }
  return(-Inf)
}

############### Run Simulations ###############

rp = 8
rq = 6
vs = gen_series(rp,rq)
vs = vs[1:10]
results3 <- vector("list", length(vs))
for(i in 1:length(vs)){
  ev = matrix(-Inf,10+1,10+1)
  print(i)
  x = vs[[i]]$series
  y = vs[[i]]$forc[1:5]
  #run our algorithm
  for(p in 0:10){
    for(q in 0:10){
      ev[p+1,q+1] = bayes_arima(x,y,p,q)
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
  results3[[i]] = c(rp-ap, rp-bp, rq-aq, rq-bq, rp + rq - (ap + aq), rp + rq - (bp+bq), arma_fits[3:length(arma_fits)])
}



