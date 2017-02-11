setwd("~/Desktop/bayes_arma")
source("all_code_get_results.R")
maxp = 10; maxq = 10;
sp500 = read.csv("number-of-earthquakes-per-year-m.csv")
dat = sp500$Number.of.earthquakes.per.year.magnitude.7.0.or.greater..1900.1998
train_split = floor(0.8*length(dat))
x = dat[1:train_split] 
mx = mean(x)
x = x - mx
y = dat[(train_split+1):length(data)]
ev = matrix(-Inf,maxp+1,maxq+1)

######################## bayes arima #######################
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
res_barma = list(mat = ev, res = c(bp, bq, rfitted_acc(x,y,bp,bq, mx), lag_terms))
saveRDS(res_barma, "eq_barma")

######################## MLE, CV arima #######################
mlos = mle_arma(x)
res_mle = list(mat = mlos$mat, res = c(mlos$ords, rfitted_acc(x,y,mlos$ords[1],mlos$ords[2],mx), mlos$coeffs))
saveRDS(res_mle, "eq_mle")
cvos = ofcv_arma(x)
res_cv = list(mat = cvos$mat, res = c(cvos$ords, rfitted_acc(x,y,cvos$ords[1],cvos$ords[2],mx), cvos$coeffs))
saveRDS(res_cv, "eq_cv")

######################## IC auto.arima #######################
a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
aicp = a$arma[1]
aicq = a$arma[2]
res_aic = c(aicp, aicq, rfitted_acc(x,y,aicp,aicq,mx), a$sigma2, a$coef)
saveRDS(res_aic, "eq_aic")

a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="aicc", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
aiccp = a$arma[1]
aiccq = a$arma[2]
res_aicc = c(aiccp, aiccq, rfitted_acc(x,y,aiccp,aiccq,mx), a$sigma2, a$coef)
saveRDS(res_aic, "eq_aicc")

a = auto.arima(x, d=0, max.p=maxp, max.q = maxq, allowmean = F, allowdrift = F, approximation = F, ic="bic", stepwise = F, max.order=maxp + maxq, stationary = T, seasonal = F)
bicp = a$arma[1]
bicq = a$arma[2]
res_bic = c(bicp, bicq, rfitted_acc(x,y,bicp,bicq,mx), a$sigma2, a$coef)
saveRDS(res_bic, "eq_bic")

################ Fitted Acc ############
rfitted_acc <- function(x, y, bp, bq, mx){
  barima_x1 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
  rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))
  fbx = forecast(barima_x1, h=length(y))
    trmse_bayes_x  = sqrt(sum((y-fbx$mean - mx)^2  ) / length(y))
  return(c(rmse_bayes_x, trmse_bayes_x))
}





















############################ IGNORE ################
# arma_fits = fitted_acc(x,y,bp,bq,mx)
# ap = arma_fits[1]
# aq = arma_fits[2]
# acc_sp500 = c(ap,aq,bp,bq,arma_fits[3:length(arma_fits)])
# barima_sp500 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
# plot_resids(barima_sp500$residuals)
# 
# fin = read.csv("financial-times-index-leading-eq.csv")
# dat = fin$Financial.Times.index.leading.equity.prices..quarterly..1960...1971
# train_split = floor(0.75*length(dat))
# x = dat[1:train_split] 
# mx = mean(x)
# x = x - mx
# y = dat[(train_split+1):length(dat)]
# ev = matrix(-Inf,5+1,5+1)
# for(p in 0:5){
#   for(q in 0:5){
#     ev[p+1,q+1] = bayes_arima(x,y,p,q)
#     if(is.nan(ev[p+1, q+1])){
#       ev[p+1, q+1] = -Inf
#     }
#   }
# }
# m = which(ev == max(ev), arr.ind = TRUE)
# bp = m[1,1] - 1
# bq = m[1,2] - 1
# arma_fits = real_fitted_acc(x,y,bp,bq,mx)
# ap = arma_fits[1]
# aq = arma_fits[2]
# acc_fin = c(ap,aq,bp,bq,arma_fits[3:length(arma_fits)])
# barima_fin = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
# plot_resids(barima_fin$residuals)
# 
# eq = read.csv("number-of-earthquakes-per-year-m.csv")
# dat = eq$Number.of.earthquakes.per.year.magnitude.7.0.or.greater..1900.1998
# train_split = floor(0.75*length(dat))
# x = dat[1:train_split] 
# mx = mean(x)
# x = x - mx
# y = dat[(train_split+1):length(dat)]
# ev = matrix(-Inf,5+1,5+1)
# for(p in 0:5){
#   for(q in 0:5){
#     ev[p+1,q+1] = bayes_arima(x,y,p,q)
#     if(is.nan(ev[p+1, q+1])){
#       ev[p+1, q+1] = -Inf
#     }
#   }
# }
# m = which(ev == max(ev), arr.ind = TRUE)
# bp = m[1,1] - 1
# bq = m[1,2] - 1
# arma_fits = real_fitted_acc(x,y,bp,bq,mx)
# ap = arma_fits[1]
# aq = arma_fits[2]
# acc_eq = c(ap,aq,bp,bq,arma_fits[3:length(arma_fits)])
# barima_eq = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
# plot_resids(barima_eq$residuals)

############### Real Fitted Accuracy #########################
# real_fitted_acc <- function(x, y, bp, bq, mx){
#   arima_x1 = auto.arima(x, d=0, max.p=5, max.q = 5, allowmean = F, approximation=T)
#   barima_x1 = arima(x, order=c(bp,0,bq), include.mean=F, method="ML")
#   a = max(arima_x1$residuals)
#   b = max(barima_x1$residuals)
#   y = na.omit(y)
#   if(a>b){
#     lines(1:length(x),x - arima_x1$residuals, col="blue")
#     plot(1:length(x), x , col="black", type='l')
#     lines(1:length(x),x - barima_x1$residuals, col="red")
# 
#   }else{
#     lines(1:length(x),x - barima_x1$residuals, col="red")
#     plot(1:length(x), x , col="black", type='l')
#     lines(1:length(x),x - arima_x1$residuals, col="blue")
#   }
#   
#   rmse_auto_x  = sqrt(sum((arima_x1$residuals)^2  ) / length(x))
#   rmse_bayes_x  = sqrt(sum((barima_x1$residuals)^2  ) / length(x))
# 
#   fax = forecast(arima_x1, h=length(y)) 
#   fbx = forecast(barima_x1, h=length(y)) 
# 
#   plot(1:length(y), (fbx$mean+mx), col="red", type='l')
#   lines(1:length(y), y , col="black")
#   lines(1:length(y),(fax$mean+mx), col="blue")
#   
#   trmse_auto_x  = sqrt(sum((y-(fax$mean + mx))^2  ) / length(x))
#   trmse_bayes_x  = sqrt(sum((y-(fbx$mean + mx))^2  ) / length(x))
# 
#   
#   return(c(arima_x1$arma[1], arima_x1$arma[2], rmse_auto_x, rmse_bayes_x, trmse_auto_x, trmse_bayes_x))
# }
