require(plot3D)

setwd('/home/wenyu/Desktop/Cornell/Class/ORIE 6741/bayes_arma')
source('all_code_get_results.R')

p = 1
q = 0
# AR(1) time series, 2 parameters tilde_sigma^2 and tilde_r1 in our model
x = arima.sim(n = 100, list(ar = c(0.8)))

xaxis <- seq(-2, 2, length.out=100)
yaxis <- seq(-2, 2, length.out=100)
M <- mesh(xaxis, yaxis)

alpha <- M$x
beta <- M$y

ff = matrix(NA,100,100)
for (i in 1:100){
  for (j in 1:100){
    ff[i,j] = f(c(alpha[i,j],beta[i,j]))
  }
}

# xaxis is sigma^2, yaxis is tilde_r1, zaxis is value of integrand
# shape suggests that Laplace approximation is reasonable
surf3D(x=alpha, y=beta, z=ff, phi=10, theta=150, colkey=FALSE, bty="b2",
       main=expression(paste(f,'(',theta,"*",')',' for',' AR(1)')))