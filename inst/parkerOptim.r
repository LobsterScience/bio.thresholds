
theta1<-c(Tm=15,Tu=27,sh=12)
x=rep(1:25,each=10)
pk = function(pars,x) ((x/pars[1])*((pars[2]-x)/(pars[2]-pars[1]))^((pars[2]-pars[1])/pars[1]))^pars[3]
y = pk(theta1,x)

ym = min(y)
y = y + runif(length(y),-.15,.15)
if(any(y<0)) y[which(y<0)]<- ym

parks<- function(params=theta1,predictors=x,data=y, dis = 'gamma') {
Tu = params[2]
Tm = params[1]
sh = params[3]
	U=(Tu-Tm)/Tm
	Z=(Tu-predictors)/(Tu-Tm)
	if(any(Z<0)) Z[which(Z<=0)]<-0.00001
	y.hat=(predictors/Tm*(Z^U))^sh
	#y.hat = pk(params,predictors)
	shape = y.hat^2 / var(data)
  	scale = var(data) / y.hat
	if(any(c(shape,scale)<=0)) return(Inf)
	#residuals = data - y.hat
	#sigma2 = sum(residuals^2) / length(data)
  	if(dis == 'gamma')	nloglike = -1 *sum((shape - 1) * log(data) - data / scale - lgamma(shape) - shape * log(scale)) #gamma with indentity link
  	#if(dis == 'norm')	nloglike = -0.5 * length(data) * log(2 * pi * sigma2) + sum(residuals^2) / (2 * sigma2)
  	return(nloglike )

}


fit = optim(par = theta1,fn = parks, data=y,predictors=x, method='L-BFGS-B',hessian=T,lower=c(5,10,6))


plot(x,y)
lines(x,pk(fit$par,x))


