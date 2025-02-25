#ratowski Ratkowsky, D.A., Lowry, R.K., McMeekin, T.A., Stokes, A.N., & Chandler, R.A. (1983). Model for bacterial culture growth rate through out the entire biokinetic temperature range. Journal of Bacteriology 154: 1222-1226.

theta <- c(d=0.1,Tl=2,g=0.2,Tu=20)
x=rep(2:20,each=2)
rat = function(pars,x) pars[1]*(x-pars[2])*(1-exp(pars[3]*(x-pars[4]))) #these function needs to be constrainted; anything beyond Tmin and Tmax is negative

y = rat(theta,x)
ym = min(y)
y = y + runif(length(y),-.15,.15)
if(any(y<0)) y[which(y<0)]<- ym

ratn ='y~a*(x-b)*(1-exp(d*(x-f)))'

da = data.frame(x=x,y=y)
da = na.omit(da)
nls_fit = nls(ratn,data=da, start=list(a=.1,b=2,d=.2,f=20))

da$y_pred <- predict(nls_fit, newdata = da)

# Plot data and fitted line
ggplot(da, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = y_pred), color = "blue") +
  labs(title = "Nonlinear Least Squares Fit",
       x = "x", y = "y")




area.under.curve<-function(fit=nls_fit){
  d1=coef(fit)[1]
  TL1=coef(fit)[2]
  g1=coef(fit)[3]
  Tu1=coef(fit)[4]
  require(MASS)
  X<-seq(TL1,Tu1,by=.1)
  rt.fn = function(x=X,dd=d1,TL=TL1,g=g1,Tu=Tu1){
    dd*(x-TL)*(1-exp(g*(x-Tu))) #these function needs to be constrainted; anything beyond Tmin and Tmax is negative
  }
  p<-rt.fn()
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

  #finding Tm
  Tm=seq(TL1,Tu1,0.01)
  X1<-log(1+g1*(Tm-TL1))+g1*(Tm-Tu1) #equation 5
  Topt<-mean(Tm[which(round(X1*100)/100==0)])

  ar<-area(rt.fn, x1, x2,fa=rt.fn(x=x1),fb =rt.fn(x=x2) ) #full area; for this to work; the first variable your your function call (rt.fn) needs to be xvariable
  ar1<-area(rt.fn, Topt,Tu1) #between Tm and Tu
  ar2<-area(rt.fn,l75,u75)
  return(c(lower75=l75,upper75=u75,Topt=Topt,Tmax=Tu1,shape=c(g1,d1),TotalAUC=ar,ToptTmaxAUC=ar1,OptArea = ar2))
}
areas<-area.under.curve(fit)


