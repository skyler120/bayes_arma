require(ggplot2)
require(ggthemes)
require(reshape2)
require(gridExtra)
require(scales)

setwd('/home/wenyu/Desktop/Cornell/Class/ORIE 6741/Results')

### using results from one time series
# plot grid-search heat map for each of six methods

# sample dataset
ev = matrix(1:121, 11)
rownames(ev) = 1:11
colnames(ev) = 1:11
ev.m = melt(ev, varnames = c('p','q'))
# heatmap where darker color = higher value
ggplot(ev.m, aes(x=p,y=q, fill=value)) +
  geom_tile(color='white', size=0.1) +
  scale_fill_gradient(low='white', high='steelblue') +
  coord_equal() +
  theme_tufte() +
  theme(text=element_text(size=20), 
        legend.position='none',
        plot.title=element_text(hjust=0.5)) +
  ggtitle('Method 1')

#setwd('/home/wenyu/Desktop/Cornell/Class/ORIE 6741/Results/ev5')
# 
ev_abarma_21 = readRDS('ev4_BARMA_approach_53')[[1]]$mat
ev_aic_21 = -readRDS('ev4_aic_approach_53')[[1]]$mat
ev_aicc_21 = -readRDS('ev4_aicc_approach_53')[[1]]$mat
ev_bic_21 = -readRDS('ev4_bic_approach_53')[[1]]$mat
ev_mle_21 = readRDS('ev4_MLE_approach_53')[[1]]$mat
ev_cv_21 = -readRDS('ev4_OFCV_approach_53')[[1]]$mat
plotpq = function(ev, methodname){
  rownames(ev) = 0:10
  colnames(ev) = 0:10
  ev.m = melt(ev, varnames = c('p','q'))
  # heatmap where darker color = higher value
  ggplot(ev.m, aes(x=p,y=q, fill=value)) +
    geom_tile(color='white', size=0.1) +
    scale_fill_gradient(low='white', high='steelblue') +
    coord_equal() +
    theme_tufte() +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_continuous(breaks=pretty_breaks()) +
    theme(text=element_text(size=15), 
          legend.position='none',
          plot.title=element_text(hjust=0.5)) +
    ggtitle(methodname)
}

# plot heatmap for all methods
plot_ev_abarma_21 = plotpq(ev_abarma_21, 'ABARMA')
plot_ev_cv_21 = plotpq(ev_cv_21, 'Cross-validation')
plot_ev_mle_21 = plotpq(ev_mle_21, 'MLE')
plot_ev_aic_21 = plotpq(ev_aic_21, 'AIC')
plot_ev_aicc_21 = plotpq(ev_aicc_21, 'AICc')
plot_ev_bic_21 = plotpq(ev_bic_21, 'BIC')
grid.arrange(plot_ev_abarma_21, plot_ev_cv_21, plot_ev_mle_21, 
             plot_ev_aic_21, plot_ev_aicc_21, plot_ev_bic_21, ncol=3)

### using results from multiple ts with different initializations
# plot shows robustness to initializations

init0 = readRDS('iter_BARMA_0init_21')
changeinit = readRDS('change_init_pm3_results_21')

plotinit = function(v, seriesno){
  mat = matrix(0, 4, 4)
  rownames(mat) = 0:3
  colnames(mat) = 0:3
  for (i in 1:10){
    mat[v[[i]]$orders[1]+1,v[[i]]$orders[2]+1] = mat[v[[i]]$orders[1]+1,v[[i]]$orders[2]+1]+1
  }
  mat.m = melt(mat, varnames = c('p','q'))
  # heatmap where darker color = higher value
  ggplot(mat.m, aes(x=p,y=q, fill=value)) +
    geom_tile(color='white', size=0.1) +
    scale_fill_gradient(low='white', high='green') +
    coord_equal() +
    theme_tufte() +
    scale_x_continuous(breaks=pretty_breaks()) +
    scale_y_continuous(breaks=pretty_breaks()) +
    theme(text=element_text(size=15), 
          legend.position='none',
          plot.title=element_text(hjust=0.5)) +
    ggtitle(paste('Series',seriesno))
}

s1 = plotinit(changeinit[[1]],1)
s2 = plotinit(changeinit[[2]],2)
s3 = plotinit(changeinit[[3]],3)
s4 = plotinit(changeinit[[4]],4)
s5 = plotinit(changeinit[[5]],5)
grid.arrange(s1,s2,s3,s4,s5,ncol=5)

# sample dataset
# order estimates using initialization = 0
p0 = 2; q0 = 1
# order estimates using random initializations
pq = matrix(floor(runif(50,1,5)), 25, 2)
colnames(pq) = c('p','q')
pq = data.frame(pq)
# counts for distance from p0
diff_p = data.frame(abs(pq$p-p0))
names(diff_p) = 'diff'
ggplot(data=diff_p, aes(diff)) +
  geom_bar() +
  theme(text=element_text(size=20),
        plot.title=element_text(hjust=0.5)) +
  xlab('differences') +
  ggtitle('|p0 - Approximated AR Lag Order with Random Initializations|')

### using results from multiple ts with different initializations
# plot shows initializations at 0 find higher modes

### using results from multiple ts
# plots true lag - approximated lag for each method

# sample dataset similar to previous
diff_p = data.frame(p0 - pq$p)
names(diff_p) = 'diff'
ggplot(data=diff_p, aes(diff)) +
  geom_bar() +
  xlim(-2,2) + # set xlim for comparison btw histograms
  theme(text=element_text(size=20),
        plot.title=element_text(hjust=0.5)) +
  xlab('differences') +
  ggtitle('True AR Lag Order - Approximated AR Lag Order')

# 
# time, p, q, true fit, model fit, true forecast, model forecast, coefficients
abarma21 = readRDS('BARMA_approach_86')
cv21 = readRDS('OFCV_approach_86')
mle21 = readRDS('MLE_approach_86')
aic21 = readRDS('aic_approach_86')
aicc21 = readRDS('aicc_approach_86')
bic21 = readRDS('bic_approach_86')
pp = 8
qq = 6

plotlag = function(res, methodname){
  diffp = pp - sapply(res, function(x){return(x[2])}) # change true p
  diffq = qq - sapply(res, function(x){return(x[3])}) # change true q
  df = data.frame(diff = c(diffp,diffq), ord = c(rep('p',25),rep('q',25)))
  ggplot(data=df, aes(diff,fill=factor(ord))) +
    geom_bar(position='dodge') +
    ylim(0,20) +
    theme(text=element_text(size=15),
          plot.title=element_text(hjust=0.5),
          legend.position = 'none') +
    scale_x_continuous(breaks=pretty_breaks(),
                       limits=c(-5,10)) + ### change limits
    xlab('differences') +
    ggtitle(methodname)
}

# plots for all methods
lag_abarma21 = plotlag(abarma21, 'ABARMA')
lag_cv21 = plotlag(cv21, 'Cross-validation')
lag_mle21 = plotlag(mle21, 'MLE')
lag_aic21 = plotlag(aic21, 'AIC')
lag_aicc21 = plotlag(aicc21, 'AICc')
lag_bic21 = plotlag(bic21, 'BIC')

# red p, blue q
grid.arrange(lag_abarma21, lag_cv21, lag_mle21, 
             lag_aic21, lag_aicc21, lag_bic21, ncol=3)

### using results from multiple ts
# plots true model fit error - model fit error for each method

# sample dataset similar to previous
errtrue = 0.5
err = runif(25,0,2)
diff_err = data.frame(errtrue - err)
names(diff_err) = 'diff'
ggplot(data=diff_err, aes(diff)) +
  geom_histogram(binwidth = 0.11) +
  xlim(-2,2) + # set xlim for comparison btw histograms
  theme(text=element_text(size=20),
        plot.title=element_text(hjust=0.5)) +
  xlab('differences') +
  ggtitle('True Model Fit Error - Method 1 Fit Error')

# fit
plotfit = function(res, methodname){
  diff = sapply(res, function(x){return(x[4])}) - sapply(res, function(x){return(x[5])})
  df = data.frame(diff = diff)
  ggplot(data=df, aes(diff)) +
    geom_histogram(binwidth = 0.05) +
    xlim(-2,0.2) +
    ylim(0,15) +
    theme(text=element_text(size=15),
          plot.title=element_text(hjust=0.5),
          legend.position = 'none') +
    xlab('differences') +
    ggtitle(methodname)
}

# plots for all methods
fit_abarma21 = plotfit(abarma21, 'ABARMA')
fit_cv21 = plotfit(cv21, 'Cross-validation')
fit_mle21 = plotfit(mle21, 'MLE')
fit_aic21 = plotfit(aic21, 'AIC')
fit_aicc21 = plotfit(aicc21, 'AICc')
fit_bic21 = plotfit(bic21, 'BIC')

grid.arrange(fit_abarma21, fit_cv21, fit_mle21, 
             fit_aic21, fit_aicc21, fit_bic21, ncol=3)

# forecast
plotforecast = function(res, methodname){
  diff = sapply(res, function(x){return(x[6])}) - sapply(res, function(x){return(x[7])})
  df = data.frame(diff = diff)
  ggplot(data=df, aes(diff)) +
    geom_histogram(binwidth = 0.05) +
    xlim(-2,1) +
    ylim(0,15) +
    theme(text=element_text(size=15),
          plot.title=element_text(hjust=0.5),
          legend.position = 'none') +
    xlab('differences') +
    ggtitle(methodname)
}

# plots for all methods
fore_abarma21 = plotforecast(abarma21, 'ABARMA')
fore_cv21 = plotforecast(cv21, 'Cross-validation')
fore_mle21 = plotforecast(mle21, 'MLE')
fore_aic21 = plotforecast(aic21, 'AIC')
fore_aicc21 = plotforecast(aicc21, 'AICc')
fore_bic21 = plotforecast(bic21, 'BIC')

grid.arrange(fore_abarma21, foer_cv21, fore_mle21, 
             fore_aic21, fore_aicc21, fore_bic21, ncol=3)

### using results from multiple ts
# plots approximation results with increasing noise

#setwd('/home/wenyu/Desktop/Cornell/Class/ORIE 6741/Results/scaling_res')

abarma_var = readRDS('BARMA2_scaling_variance_results_21')
cv_var = readRDS('OFCV2_scaling_var_21')
mle_var = readRDS('MLE2_scaling_var_21')
aic_var = readRDS('aic2_scaling_var_21')
aicc_var = readRDS('aicc2_scaling_var_21')
bic_var = readRDS('bic2_scaling_var_21')

vardf = data.frame(Method=c('ABARMA','Cross-validation','MLE','AIC','AICc','BIC'),
                   Low=numeric(6), Medium=numeric(6), High=numeric(6))
vardf[1,2:4]=rowSums(sapply(abarma_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf[2,2:4]=rowSums(sapply(cv_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf[3,2:4]=rowSums(sapply(mle_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf[4,2:4]=rowSums(sapply(aic_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf[5,2:4]=rowSums(sapply(aicc_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf[6,2:4]=rowSums(sapply(bic_var, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
vardf.m = melt(vardf, id.vars='Method')
colnames(vardf.m) = c('Method','Noise','Accuracy')
vardf.m$Method = factor(vardf.m$Method, levels=c('ABARMA','Cross-validation','MLE','AIC','AICc','BIC'))
ggplot(vardf.m, aes(Method, Accuracy)) +   
  theme(text=element_text(size=15)) +
  geom_bar(aes(fill = Noise), position = "dodge", stat="identity")

# plots approximation results with increasing length

abarma_len = readRDS('BARMA2_scaling_size_21')
cv_len = readRDS('OFCV2_scaling_size_21')
mle_len = readRDS('MLE2_scaling_size_21')
aic_len = readRDS('aic2_scaling_size_21')
aicc_len = readRDS('aicc2_scaling_size_21')
bic_len = readRDS('bic2_scaling_size_21')

lendf = data.frame(Method=c('ABARMA','Cross-validation','MLE','AIC','AICc','BIC'),
                   'L100'=numeric(6), 'L200'=numeric(6))
lendf[1,2:3]=rowSums(sapply(abarma_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf[2,2:3]=rowSums(sapply(cv_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf[3,2:3]=rowSums(sapply(mle_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf[4,2:3]=rowSums(sapply(aic_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf[5,2:3]=rowSums(sapply(aicc_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf[6,2:3]=rowSums(sapply(bic_len, function(i){return((sapply(i, function(x){return(x[2]==2 & x[3]==1)})))})/10)
lendf.m = melt(lendf, id.vars='Method')
colnames(lendf.m) = c('Method','Length','Accuracy')
lendf.m$Method = factor(lendf.m$Method, levels=c('ABARMA','Cross-validation','MLE','AIC','AICc','BIC'))
ggplot(lendf.m, aes(Method, Accuracy)) +   
  theme(text=element_text(size=15)) +
  geom_bar(aes(fill = Length), position = "dodge", stat="identity")

# plots approximated coefficients
a = readRDS('BARMA_test_coeff_21_v2')
t = readRDS('test_coeffs_11_v2')
phis = numeric(length(a))
pis = numeric(length(a))
sig = numeric(length(a))
tphis = numeric(length(a))
tpis = numeric(length(a))
for(i in 1:length(a)){
  if(a[[i]]$res[2]==1 &&  a[[i]]$res[3]==1){
    pis[i] = a[[i]]$res[length(a[[i]]$res)]
    phis[i] = a[[i]]$res[length(a[[i]]$res) - 1]
    sig[i] = a[[i]]$res[length(a[[i]]$res) - 2]
    tphis[i] = t[[i]]$ps
    tpis[i] = t[[i]]$qs
  }
}

plot(tPhis,Phis,xlab=expression(phi),ylab=expression(paste('Estimated ',phi)))
lines(-2:2,-2:2,lty=3)
plot(tPis,Pis*-1,xlab=expression(pi),ylab=expression(paste('Estimated ',pi)))
lines(-2:2,-2:2,lty=3)
Sigdf = data.frame(Sig)
ggplot(Sigdf, aes(Sig)) + 
  geom_histogram(binwidth = 0.05) +
  theme(text=element_text(size=15)) +
  xlab(expression(paste('Estimated ',sigma^2)))
hist(Sig,main='',xlab=expression(paste('Estimated ',sigma^2)))

par(mfrow=c(1,3))
plot(tPhis,Phis,xlab=expression(phi),ylab=expression(paste('Estimated ',phi)))
lines(-2:2,-2:2,lty=3)
plot(tPis,Pis*-1,xlab=expression(pi),ylab=expression(paste('Estimated ',pi)))
lines(-2:2,-2:2,lty=3)
hist(Sig,main='',xlab=expression(paste('Estimated ',sigma^2)))
abline(v=1,col='red')
par(mfrow=c(1,1))
