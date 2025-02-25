setwd('C:/Users/cooka/Documents/bio_data/bio.thresholds/BQ_Compiled/')
require(ggplot2)
require(dplyr)
require(tidyr)

x = readxl::read_xlsx('OATholds_combinedJan2025.xlsx',sheet=1)
###############################################################################################################################################
# notes
#Harrington and Hamlin Adult Heart rate
#n=18 lobster each one had HR measured as temp increased, ie temp-hr tracked through warming on each individual (lacking independence)
###############################################################################################################################################
###load this in
AUC_elliott<-function(fit=nls_fit,data=data) {
  #this ignores scalar parameter D, which is fine
  A=coef(fit)[1]
  B=coef(fit)[2]
  L=coef(fit)[3]
  require(MASS)

  X<-seq(0,L,by=.01)
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
  l75_src = 'model'
  u75_src = 'model'

  l75=a2[which.min(abs(a2$y-0.75)),1]
  u75=a3[which.min(abs(a3$y-0.75)),1]
  if(min(data$x)>l75){ l75= min(data$x); print('L75 from data'); l75_src='data'}
  if(max(data$x)<u75){ u75= max(data$x); print('U75 from data'); u75_src='data'}

  a3<-a[a$x>=l75 & a$x<=u75,]
  ot<-a3[,1]
  og<-a3[,2]

  ar<-area(ells1,A,L )
  ar1<-area(ells1, B,L)
  ar2<-area(ells1,l75,u75)
  return(c(lower75=l75,lower75_source=l75_src, upper75=u75,upper75_source=u75_src, Topt=B,Tmax=L,TL=A,TotalAUC=ar,ToptTmaxAUC=ar1,OptArea = ar2))
}
#################################################################################################

# Get Data

###################################Acid Treatment
a = subset(x,Grouping==1,select=c(Temp, Treat_Mean, Treat_Se))
a$sd = as.numeric(a$Treat_Se)*sqrt(18)
a$y = as.numeric(a$Treat_Mean)
a$x= as.numeric(a$Temp)
ggplot(a, aes(x = x, y =y,ymin=y-sd,ymax=y+sd)) +
  geom_point() +
  geom_errorbar(width=0)+
  theme_test()+
  labs(x = "Temperature", y = "Heart Rate")

gen_rep_meas <- function(mu, sigma, rho, i) {
  rms <- numeric(i)
  rms[1] <- rnorm(1, mean = mu[1], sd = sigma[1])
  for (t in 2:i) {
    rms[t] <- mu[t] + rho * (rms[t-1] - mu[t-1]) + rnorm(1, mean = 0, sd = sigma[t] * sqrt(1 - rho^2))
  }
  return(rms)
}


niter=5000

out = list()
for(k in 1:niter){
        sim_data_ac <- list()
  rho1 = runif(1, 0.01,.5)
  for(j in 1:18){
        re = gen_rep_meas(mu = a$y,sigma = a$sd, rho=rho1,i=nrow(a))
        sim_data_ac[[j]] = data.frame(ind=j,Temp = a$x,HR = re)
      }

sim_df_ac <- do.call(rbind, sim_data_ac)

##elliott on data with individual autocorrelation
da = subset(sim_df_ac,select=c(Temp,HR))
names(da) = c('x','y')
try({
nls_e_ac = nls(y~ ifelse(x<=b,d*((x-a)/(b-a)),d*((x-c)/(b-c))),data=da, start=list(a=2,b=25,c=35,d=40))

a$e_pred <- predict(nls_e_ac, newdata = a)

out[[k]] = AUC_elliott(nls_e_ac,data=a)
}, silent=T)
}

outAC = as.data.frame(do.call(rbind,out))
outAC = bio.utilities::toNums(outAC,c(1,3,5:10))
###########################################################################################################################################
#######Control

# Get Data

###################################Acid Treatment
aa = subset(x,Grouping==1,select=c(Temp, Control_Mean, Control_SE))
aa$sd = as.numeric(aa$Control_SE)*sqrt(18)
aa$Control_Mean = as.numeric(aa$Control_Mean)
aa$x = aa$Temp
aa$y = aa$Control_Mean
out=list()
for(k in 1:niter){
rho1 = runif(1,.001,.5)
# Generate time series for each individual
sim_data_co <- list()
for(j in 1:18){
  re = gen_rep_meas(mu = aa$y,sigma = aa$sd, rho=rho1,i=nrow(aa))
  sim_data_co[[j]] = data.frame(ind=j,Temp = aa$x,HR = re)
}

sim_df_co <- do.call(rbind, sim_data_co)

da = subset(sim_df_co,select=c(Temp,HR))
names(da) = c('x','y')
try({
nls_e_co = nls(y~ ifelse(x<=b,d*((x-a)/(b-a)),d*((x-c)/(b-c))),data=da, start=list(a=-1,b=26,c=60,d=100))
aa$x = aa$Temp
aa$e_pred <- predict(nls_e_co, newdata = aa)

out[[k]] = AUC_elliott(nls_e_co,data=aa)},
silent = T)
}
outCo = as.data.frame(do.call(rbind,out))
outCo = bio.utilities::toNums(outCo,c(1,3,5:10))



outAC$Group = 'Acid'
outCo$Group = 'Control'

ab = rbind(outCo,outAC)

a1 = ggplot(ab, aes(x =lower75,group=Group,colour=Group)) +
  geom_density(position='identity') +
  theme_test(base_size = 14)+labs(x='75% Max HR')


b1 = ggplot(ab, aes(x =Topt.b,group=Group,colour=Group)) +
  geom_density(position='identity') +
  theme_test(base_size = 14)+labs(x='Max HR')

gridExtra::grid.arrange(a1,b1)
