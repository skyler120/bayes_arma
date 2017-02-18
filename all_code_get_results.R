require(forecast)

############################# Recursive Function for phi and pi ###############################################

phiphi = function(k,r){
  phimat = matrix(0,k,k)
  for(j in 1:k){
    for(i in 1:j){
      if(i==j){
        phimat[j,i] =tanh(r[j])
      }
      else{
        phiipre = phimat[j-1,i]
        phiki = phimat[j-1,j-i]
        phimat[j,i] = phiipre - tanh(r[j])*phiki
      }
    }
  }
  return(phimat[k,])
}

############################# Recursive Function for mu_t ###############################################

muts <- function(phir,pis){
  if(p+q == 0){
    return(rep(0,length(x)))
  }
  else{
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
}


muts_bpbq <- function(phir,pis){
  poq = max(length(phir), length(pis))
  mus = rep(0,poq)
  xapp = c(rep(0,poq), x)
  for(t in (poq+1):(poq+length(x))){
    if(bp==0){
      p1 = 0
    }
    else{
      p1 = sum(phir*rev(xapp[(t-bp):(t-1)]))
    }
    if(bq==0){
      p2 = 0
    }
    else{
      p2 = sum(pis*(rev(xapp[(t-bq):(t-1)]) - mus[(t-bq):(t-1)]))
    }
    mus = c(mus, p1-p2)
  }
  return(mus)
}

#############################  Function for logp ###############################################

logp <- function(params){
  #print(x)
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

logp_bpbq <- function(params){
  n = length(x)
  poq = max(bp,bq)
  phir = 0
  pis = 0
  tsigsq = params[1]
  if(bp>0){
    phir = phiphi(bp,params[2:(1+bp)])
    if(bq>0){
      pis = phiphi(bq,params[(2+bp):length(params)])
    }
  }else if(bq>0){
    pis = phiphi(bq,params[(2+bp):length(params)])
  }else{
    phir = 0
    pis = 0
  }
  mut = muts_bpbq(phir,pis)
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
  #}else if (p==1){
    #Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2))
    #Jphir = 1
  }else{
    Jphir = 1
  }
  # J_pis
  if (q>=2){
    Jpis = prod((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2), 1-tanh(params[1+p+2*(1:floor(q/2))]))
  #}else if (q==1){
    #Jpis = prod((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2))
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
    Jr = 0
  }
  # J_s
  if (q>0){
    Js = sum(log(1-tanh(params[(2+p):(1+p+q)])^2))
  }else{
    Js = 0
  }
  # J_phir
  if (p>=2){
    Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2), 1-tanh(params[1+2*(1:floor(p/2))]))
  #}else if (p==1){
    #Jphir = prod((1-tanh(params[2:(1+p)])^2)^floor(((1:p)-1)/2))
  }else{
    Jphir = 0
  }
  # J_pis
  if (q>=2){
    Jpis = sum(log((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2)), log(1-tanh(params[1+p+2*(1:floor(q/2))])))
  #}else if (q==1){
    #Jpis = sum(log((1-tanh(params[(2+p):(1+p+q)])^2)^floor(((1:q)-1)/2)))
  }else{
    Jpis = 0
  }
  # J_sigma
  Jsigma = params[1]
  res = -pfunc-Jr-Js-Jphir-Jpis-Jsigma
  #if(res>10^(308)){return(10^308)
  #}else if(res<-Inf){return(-10^(308))}
  return(res)
}

############### Compute Accuracy ###############
rfitted_acc <- function(x, y, bp, bq, mx){
  y = na.omit(y)
  barima_x1 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
  rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))
  fbx = forecast(barima_x1, h=length(y))
  trmse_bayes_x  = sqrt(sum((y-fbx$mean - mx)^2  ) / length(y))
  return(c(rmse_bayes_x, trmse_bayes_x))
}

fitted_acc <- function(x, y, bp, bq, rp, rq){
  barima_x1 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
  arima_x1 = arima(x, order = c(rp,0,rq), include.mean=F, method="ML")
  
  rmse_auto_x  = sqrt(sum((arima_x1$residuals)^2  ) / length(x))
  rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))
  
  fax = forecast(arima_x1, h=length(y))
  fbx = forecast(barima_x1, h=length(y))

  trmse_auto_x  = sqrt(sum((y-fax$mean)^2  ) / length(y))
  trmse_bayes_x  = sqrt(sum((y-fbx$mean)^2  ) / length(y))
  return(c(rmse_auto_x, rmse_bayes_x, trmse_auto_x, trmse_bayes_x))
}

############### Generate Sample Data ###############

gen_arm <- function(samp, pp, qq, noise_level){
  train = floor(0.8*samp)
  if((pp+qq)>0){
    phis = runif(pp,-1,1)
    pis = runif(qq,-1,1)
    x = arima.sim(n = samp, list(ar = phis, ma = pis), sd = noise_level)
    return(list(ps = phis, qs = pis, series=x[1:train],forc = na.omit(x[train+1:samp]), p = pp, q = qq))
  }
}

gen_coeffs <- function(samp, pp, qq, noise_level){
  train = floor(0.8*samp)
  if((pp+qq)>0){
    phis = runif(pp,-1,1)
    pis = runif(qq,-1,1)
    x = arima.sim(n = samp, list(ar = phis, ma = pis), sd = noise_level)
    return(list(ps = phis, qs = pis))
  }
}

gen_noise <- function(num_series, samp, noises, pp,qq){
  i=1
  train = floor(0.8*samp)
  res <- vector("list", num_series)
  while(i<=num_series){
    temp = try(gen_coeffs(samp, pp, qq, noises[1]))
    if (typeof(temp)=="character"){
      temp  = NULL
    }
    else{
      res5 <- vector("list", length(noises))
      for(j in 1:length(noises)){
        x = arima.sim(n = samp, list(ar = temp$ps, ma = temp$qs), sd = noises[j])
        res5[[j]] = list(series=x[1:train],forc = na.omit(x[train+1:samp]), p = pp, q = qq)
      }
      res[[i]] = res5
      i <- i+1
    }
  }
  return(res)
}

gen_series <- function(num_series, sampl, noise_level, pp,qq){
  i = 1
  res <- vector("list", num_series)
  while(i<=num_series){
    res[[i]] = try(gen_arm(sampl, pp, qq, noise_level))
    if (typeof(res[[i]])=="character"){
      res[[i]]  = NULL
    }
    else{
      i <- i+1
    }
  }
  return(res)
}

# gen_var_samp_series <- function(num_series, samps,pp,qq){
#   i = 1
#   noise_level = 1
#   res <- vector("list", length(samps))
#   for(i in 1:length(samps)){
#     res[[i]] = gen_series(num_series, samps[i], noise_level, pp,qq)
#   }
#   return(res)
# }

# gen_var_noise_series <- function(num_series, noises, pp,qq){
#   i = 1
#   samp = 125
#   res <- vector("list", length(noises))
#   for(i in 1:length(noises)){
#     res[[i]] = try(gen_series(num_series, samp, noises[i], pp, qq))
#   }
#   return(res)
# }

############### Bayesian ARMA Model Order Determination ###############

bayes_arima <- function(x,y, p, q, init){
  #if(p+q>0){
    #init = rep(0,p+q+1)
    res = try(optim(init, logf_neg, method="L-BFGS-B",control = list(maxit=1000), hessian = T))
    if(typeof(res)=="character"){
      res = optim(init, logf_neg, method="SANN", hessian = T)
      print("Done with SANN iter")
    }
    params = res$par
    hess = -res$hessian
    #hess = logpTheta.H(params) + diag(logJr.H(params) + logJs.H(params) + logJphir.H(params) + logJpis.H(params))
    result = f(res$par)*(2*pi)^(dim(hess)[1]/2)*det(-hess)^(-1/2)
    if(is.nan(result)){ result = -Inf}
      return(list(transformed_params = res$par, pdm = result))
  #}
  #return(list(transformed_params = 0, pdm =  -Inf))
}

############### MLE Order Determination ###############

mle_arma <- function(x){
  ev = matrix(-Inf,maxp+1,maxq+1)
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = try(arima(x, order=c(p,0,q), include.mean=F, method="ML"))
      if(typeof(a)=="character"){ ev[p+1, q+1] = -Inf}
      else {ev[p+1,q+1] = a$loglik}
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  a = arima(x, order=c(mlp,0,mlq), include.mean=F, method="ML")
  return(list(ords = c(mlp, mlq), coeffs = a$coef, mat = ev))
}

############### OFCV Order Determination ###############

ofcv_arma <- function(x){
  train = x[1:floor(length(x)*0.8)]
  test = x[floor(length(x)*0.8):length(x)]
  ev = matrix(Inf,maxp+1,maxq+1)
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = try(arima(x, order=c(p,0,q), include.mean=F, method="ML"))
      if(typeof(a)=="character"){
        ev[p+1, q+1] = Inf
      }
      else {
        fax = forecast(a, h=length(test))
        ev[p+1,q+1] = sqrt(sum((test-fax$mean)^2  ) / length(test))
      }

    }
  }
  m = which(ev == min(ev), arr.ind = TRUE)
  cvp = m[1,1] - 1
  cvq = m[1,2] - 1
  a = arima(train, order=c(cvp,0,cvq), include.mean=F, method="ML")
  return(list(ords = c(cvp, cvq), coeffs = a$coef, mat = ev))
}

############### AIC Order Determination ###############
aic_arma <- function(x){
  ev = matrix(-Inf,maxp+1,maxq+1)
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = try(arima(x, order=c(p,0,q), include.mean=F, method="ML"))
      if(typeof(a)=="character"){ ev[p+1, q+1] = -Inf}
      else {ev[p+1,q+1] = 2*length(a$coef) +2 - 2*a$loglik}
    }
  }
  m = which(ev == min(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  a = arima(x, order=c(mlp,0,mlq), include.mean=F, method="ML")
  return(list(ords = c(mlp, mlq), coeffs = a$coef, mat = ev))
}

############### AICc Order Determination ###############
aicc_arma <- function(x){
  ev = matrix(-Inf,maxp+1,maxq+1)
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = try(arima(x, order=c(p,0,q), include.mean=F, method="ML"))
      if(typeof(a)=="character"){ ev[p+1, q+1] = -Inf}
      else {ev[p+1,q+1] = (2*length(a$coef) +2 - 2*a$loglik) + (2*(length(a$coef)+2)*(length(a$coef)+3))/(length(x) - length(a$coef) - 3)}
    }
  }
  m = which(ev == min(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  a = arima(x, order=c(mlp,0,mlq), include.mean=F, method="ML")
  return(list(ords = c(mlp, mlq), coeffs = a$coef, mat = ev))
}

############### BIC Order Determination ###############
bic_arma <- function(x){
  ev = matrix(-Inf,maxp+1,maxq+1)
  for(p in 0:maxp){
    for(q in 0:maxq){
      a = try(arima(x, order=c(p,0,q), include.mean=F, method="ML"))
      if(typeof(a)=="character"){ ev[p+1, q+1] = -Inf}
      else {ev[p+1,q+1] = log(length(x))*(length(a$coef) + 1) - 2*a$loglik}
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  mlp = m[1,1] - 1
  mlq = m[1,2] - 1
  a = arima(x, order=c(mlp,0,mlq), include.mean=F, method="ML")
  return(list(ords = c(mlp, mlq), coeffs = a$coef, mat = ev))
}

############### Find untransformed coeffs ###############
r_to_phi <- function(r){
  k = length(r)
  phi = numeric(k)
  phiprev = numeric(k) 
  for(i in 1:k){
    phi[i] = r[i]
    if(i>1){
      phi[1:(i-1)] = phiprev[1:(i-1)] - r[i]*phiprev[(i-1):1]
    }
    phiprev = phi
  }
  return(phi)
}

get_coeffs<- function(transformed_params, bp, bq){
  params = numeric(length(transformed_params))
  params[1] = exp(transformed_params[1])
  params[2:(bp+1)] = r_to_phi(tanh(transformed_params[2:(bp+1)]))
  params[(bp+2):length(transformed_params)] = r_to_phi(tanh(transformed_params[(bp+2):length(transformed_params)]))
  return(params)
}

############### Find Coeffs of phis and pis (IGNORE) ###############
# joint_data_model_coeffs <- function(params){
#   params[1] = log(params[1])
#   return(exp(logp_bpbq(params)))
# }
# 
# find_coeffs <- function(init_params, bp, bq){
#   res = optim(init_params, joint_data_model_coeffs, method = "L-BFGS-B", lower = c(0, rep(-Inf,bp+bq)), upper = rep(Inf, bp+bq+1))
#   return(res)
# }

