sp500 = read.csv("quarterly-sp-500-index-19001996.csv")
dat = sp500$Quarterly.S.P.500.index..1900.1996
train_split = floor(0.75*length(dat))
x = dat[1:train_split] 
x = x - mean(x)
y = dat[(train_split+1):length(data)]
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
arma_fits = real_fitted_acc(x,y,bp,bq)
ap = arma_fits[1]
aq = arma_fits[2]
acc_sp500 = c(ap,aq,bp,bq,arma_fits[3:length(arma_fits)])

fin = read.csv("financial-times-index-leading-eq.csv")
dat = fin$Financial.Times.index.leading.equity.prices..quarterly..1960...1971
acc_fin = real_data_acc(dat)

eq = read.csv("number-of-earthquakes-per-year-m.csv")
data = eq$Number.of.earthquakes.per.year.magnitude.7.0.or.greater..1900.1998
acc_eq = real_data_acc(data)


############### Real Fitted Accuracy #########################
real_fitted_acc <- function(x, y, p, q){
  arima_x1 = auto.arima(x, d=0, max.p=5, max.q = 5, allowmean = F)
  barima_x1 = arima(x, order=c(p,0,q), include.mean=F)
  a = max(arima_x1$residuals)
  b = max(barima_x1$residuals)
  
  if(a>b){
    plot(1:length(x),x - arima_x1$residuals, col="blue", type='l')
    lines(1:length(x), x , col="black")
    lines(1:length(x),x - barima_x1$residuals, col="red")

  }else{
    plot(1:length(x),x - barima_x1$residuals, col="red", type='l')
    lines(1:length(x), x , col="black")
    lines(1:length(x),x - arima_x1$residuals, col="blue")
  }
  
  rmse_auto_x  = sqrt(sum((arima_x1$residuals)^2  ) / length(x))
  rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))

  fax = forecast(arima_x1, h=length(y)) 
  fbx = forecast(barima_x1, h=length(y)) 

  plot(1:length(y), fbx$mean + mean(x), col="red", type='l')
  lines(1:length(y), y , col="black")
  lines(1:length(y),fax$mean + mean(x), col="blue")
  
  trmse_auto_x  = sqrt(sum((y-(fax$mean+mean(x)))^2  ) / length(x))
  trmse_bayes_x  = sqrt(sum((y-(fbx$mean + mean(x)))^2  ) / length(x))

  
  return(c(arima_x1$arma[1], arima_x1$arma[2], rmse_auto_x, rmse_bayes_x, trmse_auto_x, trmse_bayes_x))
}
