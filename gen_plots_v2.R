require(ggplot2)
require(reshape2)

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

par(mfrow=c(3,2))
# plot heatmap for all methods
par(mfrow=c(1,1))


### using results from multiple ts with different initializations
# plot shows robustness to initializations

int.hist = function(x,ylab="Frequency",...) {
  barplot(table(factor(x,levels=min(x):max(x))),space=0,xaxt="n",ylab=ylab,...);axis(1)
}

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
  ggtitle('|True AR Lag Order - Approximated AR Lag Order|')

par(mfrow=c(3,2))
# plots for all methods
par(mfrow=c(1,1))
