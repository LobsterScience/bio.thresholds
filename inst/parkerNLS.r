require(ggplot2)

theta1<-c(Tm=15,Tu=24,sh=2)
x=rep(1:25,each=2)
pk = function(pars,x) ((x/pars[1])*((pars[2]-x)/(pars[2]-pars[1]))^((pars[2]-pars[1])/pars[1]))^pars[3]
y = pk(theta1,x)
ym = min(y)
y = y + runif(length(y),-.15,.15)
if(any(y<0)) y[which(y<0)]<- ym

pkn ='y~((x/a*(b-x)/(b-a))^((b-a)/a))^c'

da = data.frame(x=x,y=y)
da = na.omit(da)
nls_fit = nls(pkn,data=da, start=list(a=15,b=25,c=6))

da$y_pred <- predict(nls_fit, newdata = da)

# Plot data and fitted line
ggplot(da, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = y_pred), color = "blue") +
  labs(title = "Nonlinear Least Squares Fit",
       x = "x", y = "y")


area.under.curve<-function(fit=nls_fit){
  A=coef(fit)[1]
  B=coef(fit)[2]
  L=coef(fit)[3]
  require(MASS)
  X<-seq(0,B,by=.1)
  pk.fn<-function(x=X,a=A,b=B,c=L){
    ((x/a*(b-x)/(b-a))^((b-a)/a))^c
  }
  p<-pk.fn()
  w<-which(p>0) #0 can change this to remove tails
  x1<-X[min(w)]
  x2<-X[max(w)]
  a = data.frame(x=X,y=p)
  a2 = subset(a, x<A)
  a3 = subset(a, x>A)

  l75=a2[which.min(abs(a2$y-0.75)),1]
  u75=a3[which.min(abs(a3$y-0.75)),1]
  a3<-a[a$x>=l75 & a$x<=u75,]
  ot<-a3[,1]
  og<-a3[,2]

  ar<-area(pk.fn, x1, x2) #full area
  ar1<-area(pk.fn, A,B) #between Tm and Tu
  ar2<-area(pk.fn,l75,u75)
  return(c(lower75=l75,upper75=u75,Topt=A,Tmax=B,shape=L,TotalAUC=ar,ToptTmaxAUC=ar1,OptArea = ar2))
}
areas<-area.under.curve(fit)

