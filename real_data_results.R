dat = read.csv()
data = dat$
 
real_data_acc <- function(data){ 
  train_split = floor(0.75*length(data))
  x = data[1:train_split]
  y = data[train_split+1:length(data)]
  for(p in 0:5){
    for(q in 0:5){
      ev[p+1,q+1] = bayes_arima(x,y,p,q)
      if(is.nan(ev[p+1, q+1])){
        ev[p+1, q+1] = -Inf
      }
    }
  }
  m = which(ev == max(ev), arr.ind = TRUE)
  bp = m[1,1] - 1
  bq = m[1,2] - 1
  arma_fits = fitted_acc(x,y,bp,bq,rp,rq)
  ap = arma_fits[1]
  aq = arma_fits[2]
  acc = c(rp-ap, rp-bp, rq-aq, rq-bq, rp + rq - (ap + aq), rp + rq - (bp+bq), arma_fits[3:length(arma_fits)])
  return(acc)  
}