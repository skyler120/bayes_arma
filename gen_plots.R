ARA = c()
ARB = c()
MAA = c()
MAB = c()
ALLA = c()
ALLB = c()
rmsea = c()
rmseb = c()
rmset = c()
trmsea = c()
trmseb = c()
trmset = c()
sum_stats = rep(0,16)


for(i in 1:length(results)){
  ARA = c(ARA,results[[i]][1])
  ARB = c(ARB,results[[i]][2])
  MAA = c(MAA,results[[i]][3])
  MAB = c(MAB,results[[i]][4])
  ALLA = c(ALLA,results[[i]][5])
  ALLB = c(ALLB,results[[i]][6])
  rmsea = c(rmsea,results[[i]][7])
  rmseb = c(rmseb,results[[i]][8])
  rmset = c(rmset,results[[i]][9])
  trmsea = c(trmsea,results[[i]][10])
  trmseb = c(trmseb,results[[i]][11])
  trmset = c(trmset,results[[i]][12]) 
}


sum_stats[1] = sum((ARA>0)) #under ar
sum_stats[2] = sum((ARA<0)) #over ar
sum_stats[3] =sum((ARA==0)) #correct ar
sum_stats[4] =sum((MAA>0)) #under ma
sum_stats[5] =sum((MAA<0)) #over ma
sum_stats[6] =sum((MAA==0)) #correct ma 

sum_stats[7] = sum((ARB>0)) #under ar
sum_stats[8] = sum((ARB<0)) # over ar
sum_stats[9] = sum((ARB==0))#correct ar
sum_stats[10] = sum((MAB>0)) #under ma
sum_stats[11] = sum((MAB<0)) #over ma
sum_stats[12] =sum((MAB==0))#correct ma

sum_stats[13] = mean(rmset - rmsea)
sum_stats[14] = mean(rmset - rmseb)
sum_stats[15] = mean(trmset - trmsea)
sum_stats[16] = mean(trmset - trmseb)


#df = data.frame(rmset-rmsea)
#ggplot(data = df, aes(df$rmset...rmsea), col=I("blue")) + geom_histogram()
#qplot(df$rmset...rmsea, binwidth=0.1, main = "Histogram of Model Fit Error", xlab = "RMSE of True Model - RMSE of auto.arima", fill = I("blue"),col=I("black"))

#df = data.frame(rmset-rmseb)
#ggplot(data = df, aes(df$rmset...rmseb), col=I("blue")) + geom_histogram()
#qplot(df$rmset...rmseb, binwidth=0.1, main = "Histogram of Model Fit Error", xlab = "RMSE of True Model - RMSE of ABARMA Model", fill = I("blue"),col=I("black"))

#df = data.frame(trmset-trmsea)
#ggplot(data = df, aes(df$trmset...trmsea), col=I("blue")) + geom_histogram()
#qplot(df$trmset...trmsea, binwidth=0.1, main = "Histogram of Forecast Error", xlab = "RMSE of True Model - RMSE of auto.arima", fill = I("blue"),col=I("black"))

#df = data.frame(trmset-trmseb)
#ggplot(data = df, aes(df$trmset...trmseb), col=I("blue")) + geom_histogram()
#qplot(df$trmset...trmseb, binwidth=0.1, main = "Histogram of Forecast Error", xlab = "RMSE of True Model - RMSE of ABARMA Model", fill = I("blue"),col=I("black"))


methods = c(rep("Model Fit Error - auto.arima", length(ARA)), rep("Model Fit Error - ABARMA", length(ARB)))
differences = c(rmset - rmsea, rmset - rmseb)
d = data.frame(methods, differences)

p = ggplot(d, aes(differences, fill = methods)) + geom_histogram(binwidth=0.01,position = 'identity', alpha = .3) 
p + labs(title = "Model Fit Error of ABARMA and auto.arima")

methods = c(rep("Forecast Error - auto.arima", length(ARA)), rep("Forecast Error - ABARMA", length(ARB)))
differences = c(trmset - trmsea, trmset - trmseb)
d = data.frame(methods, differences)

p = ggplot(d, aes(differences, fill = methods)) + geom_histogram(binwidth=0.01,position = 'identity', alpha = .3) 
p + labs(title = "Forecast Error of ABARMA and auto.arima")

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

