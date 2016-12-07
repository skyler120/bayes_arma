sp500 = read.csv("quarterly-sp-500-index-19001996.csv")
data = sp500$Quarterly.S.P.500.index..1900.1996
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
acc_sp500 = c(rp-ap, rp-bp, rq-aq, rq-bq, rp + rq - (ap + aq), rp + rq - (bp+bq), arma_fits[3:length(arma_fits)])

fin = read.csv("financial-times-index-leading-eq.csv")
data = fin$Financial.Times.index.leading.equity.prices..quarterly..1960...1971
acc_fin = real_data_acc(data)

eq = read.csv("number-of-earthquakes-per-year-m.csv")
data = eq$Number.of.earthquakes.per.year.magnitude.7.0.or.greater..1900.1998
acc_eq = real_data_acc(data)


