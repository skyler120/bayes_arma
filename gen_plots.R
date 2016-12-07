ARA = c()
ARB = c()
MAA = c()
MAB = c()
ALLA = c()
ALLB = c()
rmsea = c()
rmseb = c()
for(i in 1:length(results3)){
  ARA = c(ARA,results3[[i]][1])
  ARB = c(ARB,results3[[i]][2])
  MAA = c(MAA,results3[[i]][3])
  MAB = c(MAB,results3[[i]][4])
  ALLA = c(ALLA,results3[[i]][5])
  ALLB = c(ALLB,results3[[i]][6])
  rmsea = c(rmsea,results3[[i]][7])
  rmseb = c(rmseb,results3[[i]][8])
}
rmsea
mean(rmsea - rmseb) #-0.00135 0.1199275
mean(c(rmsea[1:3], rmsea[5:length(rmsea)]) - c(rmseb[1:3], rmseb[5:length(rmseb)])) #0.15, 0.01413537
sum((rmsea-rmseb)>0) #13 6 11

sum((ARA>0)) #24 13 25
sum((ARA<0)) #1 2 0 
sum((MAA>0)) #17 10 24
sum((MAA<0)) #3 3 0

sum((ARB>0)) #20 9 25
sum((ARB<0)) #0 3 0 
sum((MAB>0)) #16 11 25
sum((MAB<0)) #4 2 0



#hist(rmsea - rmseb, breaks = seq(-2,5,by=1),xlab = "RMSE auto.arima - RMSE Bayesian", main = "RMSE of auto.arima - RMSE of Bayesian ARMA for ARMA(5,3)")
df = data.frame(rmsea-rmseb)
ggplot(data = df, aes(df$rmsea...rmseb), col=I("blue")) + geom_histogram()
qplot(df$rmsea...rmseb, binwidth=0.25, main = "Histogram of RMSE of auto.arima - RMSE of ABARMA", xlab = "Difference in RMSE", fill = I("blue"),col=I("black"))

methods = c(rep("auto.arima", length(ARA)), rep("ABARMA", length(ARB)))
differences = c(ARA, ARB)
d = data.frame(methods, differences)

p = ggplot(d, aes(differences, fill = methods)) + geom_bar(position = 'identity', alpha = .3) 
p + labs(title = "True AR Lag Order - Approximated AR Lag Order")

methods = c(rep("auto.arima", length(MAA)), rep("ABARMA", length(MAB)))
differences = c(MAA, MAB)
d = data.frame(methods, differences)

p = ggplot(d, aes(differences, fill = methods)) + geom_bar(position = 'identity', alpha = .3) 
p + labs(title = "True MA Lag Order - Approximated MA Lag Order")

