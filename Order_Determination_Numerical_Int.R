############### Order Determination through Numerical Integration ###############
# (still Bayes formulation)

# create function where there is no reparametrization from hypercubes to R

phiphi2 = function(k,r){
  phimat = matrix(0,k,k)
  for(j in 1:k){
    for(i in 1:j){
      if(i==j){
        phimat[j,i] =r[j]
      }
      else{
        phiipre = phimat[i,j-1]
        phiki = phimat[j-i,j-i]
        phimat[j,i] = phiipre - r[j]*phiki
      }
    }
  }
  return(phimat[k,])
}

logp2 <- function(params){
  n = length(x)
  poq = max(p,q)
  phir = 0
  pis = 0
  tsigsq = params[1]
  if(p>0){
    phir = phiphi2(p,params[2:(1+p)])
    if(q>0){
      pis = phiphi2(q,params[(2+p):length(params)])
    }
  }else if(q>0){
    pis = phiphi2(q,params[(2+p):length(params)])
  }else{
    phir = 0
    pis = 0
  }
  mut = muts(phir,pis)
  val = sum((x-mut[(poq+1):length(mut)])^2)
  return(-1*(n+2)/2*tsigsq - 1/2*exp(-tsigsq)*val)
}

# integrand function
f2 = function(params){
  # p(D,theta*|M), off by a constant 
  pfunc = exp(logp2(params))
  # J_phir
  if (p>=2){
    Jphir = prod((1-params[2:(1+p)]^2)^floor(((1:p)-1)/2), 1-params[1+2*(1:floor(p/2))])
  }else if (p==1){
    Jphir = prod((1-params[2:(1+p)]^2)^floor(((1:p)-1)/2))
  }else{
    Jphir = 1
  }
  # J_pis
  if (q>=2){
    Jpis = prod((1-params[(2+p):(1+p+q)]^2)^floor(((1:q)-1)/2), 1-params[1+p+2*(1:floor(q/2))])
  }else if (q==1){
    Jpis = prod((1-params[(2+p):(1+p+q)]^2)^floor(((1:q)-1)/2))
  }else{
    Jpis = 1
  }
  # J_sigma
  Jsigma = exp(params[1])
  return(pfunc*Jphir*Jpis*Jsigma)
}

require(cubature)

hcubature(f=f2, lowerLimit=c(-10,rep(-1,p+q)), upperLimit=c(10,rep(1,p+q)), 
          tol = 1e-05, fDim = 1,
          maxEval = 0, absError = 0, doChecking = FALSE,
          vectorInterface = FALSE, norm = 'L2')