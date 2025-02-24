#simulate data using ratkowsky
theta <- c(d=0.1,Tl=2,g=0.2,Tu=20)
x=rep(2:20,each=2)
rat = function(pars,x) pars[1]*(x-pars[2])*(1-exp(pars[3]*(x-pars[4]))) #these function needs to be constrainted; anything beyond Tmin and Tmax is negative

y = rat(theta,x)
ym = min(y)
y = y + runif(length(y),-.15,.15)
if(any(y<0)) y[which(y<0)]<- ym


da = data.frame(x=x,y=y)
da = na.omit(da)

elliott_func = function(x,a,b,c){
  ifelse(x<b,((x-a)/(b-a)),((x-c)/(b-c)))
}
nls_fit = nls(y~ifelse(x<=b,((x-a)/(b-a)),((x-c)/(b-c))),data=da, start=list(a=2,b=15,c=20))

da$y_pred <- predict(nls_fit, newdata = da)

ggplot(da, aes(x = x, y = y)) +
  geom_point() +
  geom_line(aes(y = y_pred), color = "blue") +
  labs(title = "Nonlinear Least Squares Fit",
       x = "x", y = "y")


area.under.curve<-function(fit=nls_fit) {
  A=coef(fit)[1]
  B=coef(fit)[2]
  L=coef(fit)[3]
  require(MASS)

  X<-seq(0,L,by=.1)
  ells1<-function(o=X,Tm=B,Tl=A,Tu=L) {
    y=rep(NA,length=length(o))
    for(i in 1:length(o)) {
      if(o[i]<=Tm) {
        y[i]<-	((o[i]-Tl)/(Tm-Tl))
      }
      if(o[i]>Tm) {
        y[i]<-	((o[i]-Tu)/(Tm-Tu))
      }
    }
    y
  }
  p<-ells1()
  w<-which(p>0) #0 can change this to remove tails
  x1<-X[min(w)]
  x2<-X[max(w)]
  a = data.frame(x=X,y=p)
  a2 = subset(a, x<B)
  a3 = subset(a, x>B)

  l75=a2[which.min(abs(a2$y-0.75)),1]
  u75=a3[which.min(abs(a3$y-0.75)),1]
  a3<-a[a$x>=l75 & a$x<=u75,]
  ot<-a3[,1]
  og<-a3[,2]

  ar<-area(ells1,A,L )
  ar1<-area(ells1, B,L)
  ar2<-area(ells1,l75,u75)
  return(c(lower75=l75,upper75=u75,Topt=B,Tmax=L,TL=A,TotalAUC=ar,ToptTmaxAUC=ar1,OptArea = ar2))
}

area.under.curve(nls_fit)
