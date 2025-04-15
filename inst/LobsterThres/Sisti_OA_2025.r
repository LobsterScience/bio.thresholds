setwd('C:/Users/cooka/Documents/bio_data/bio.thresholds/BQ_Compiled/Sisti_OA_2025/')
require(ggplot2)
require(dplyr)
require(tidyr)

###############################################################################################################################################
# notes
#Sisti Multiple Variables; 5 females, eggs excised from females, measurements were taken on 3 embryos from each brood at each time point and treatment and measurement

###############################################################################################################################################
###load this in
AUC_elliott<-function(fit=nls_fit,data=data) {
  #this ignores scalar parameter D, which is fine
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
# O2
#################################################################################################

a = read.csv('Sisti_OA_O2.csv')

a$y = as.numeric(a$Oxygen_Consumption_Rate_OCR )
a$x= as.numeric(a$pH_range_midpoint)
a$pei = a$PEI_microns
a$Hion = 10^-a$x
breaks <- c(0, 100, 200, 300, Inf)
labels <- c("0-100", "101-200", "201-300", ">300")
# Apply the cut function
a$g_pei <- cut(a$pei, breaks = breaks, labels = labels, right = TRUE)


#o2 consump by PEI and pH
ggplot(a, aes(x = x, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "pH", y = "O2 Comp")

#by pH
ggplot(a, aes(x = pei, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "PEI", y = "O2 Comp")+
  facet_wrap(~x)+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)

#by pei group
ggplot(a, aes(x = Hion, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "H+", y = "O2 Comp")+
  facet_wrap(~g_pei)+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)
##############################################################
###conclusion
#nothing to analyse that would be meaningful for o2 Comp
###############################################################

######################################################################################################################

#################################################################################################
# FRAP
#################################################################################################

a = read.csv('Sisti_OA_FRAP.csv')

a$y = as.numeric(a$Ferric.reducing_antioxidant_power_FRAP )
a$x= as.numeric(a$pH_range_midpoint)
a$pei = a$PEI_microns
a$Hion = 10^-a$x
breaks <- c(0, 100, 200, 300, Inf)
labels <- c("0-100", "101-200", "201-300", ">300")
# Apply the cut function
a$g_pei <- cut(a$pei, breaks = breaks, labels = labels, right = TRUE)


#FRAP consump by PEI and pH
ggplot(a, aes(x = x, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "pH", y = "FRAP")

#by pH
ggplot(a, aes(x = pei, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "PEI", y = "FRAP")+
  facet_wrap(~x)+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)

#by pei group
ggplot(a, aes(x = Hion, y =y)) +
  geom_point(aes(colour=pei)) +
  theme_test()+
  labs(x = "H+", y = "FRAP")+
  facet_wrap(~g_pei)+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)
##############################################################
###conclusion
#nothing to analyse that would be meaningful for FRAP
###############################################################
#################################################################################################
# NAK
#################################################################################################

a = read.csv('Sisti_OA_NAK.csv')

a$y = as.numeric(a$Na_K_ATPase_activity )
a$x= as.numeric(a$pH_fr_plot)
a$Hion = 10^-a$x

#NAK and pH
ggplot(subset(a), aes(x = x, y =y)) +
  geom_point() +
  theme_test()+
  labs(x = "pH", y = "Na K ATPase")+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)+
  geom_smooth(data=subset(a,x<8),method='gam',formula=y~s(x,k=3),se=FALSE)

##############################################################
###conclusion
#some differences from control but essentially flat after pH drops below 8
###############################################################

#################################################################################################
# HR
#################################################################################################

  a = read.csv('Sisti_OA_Pro.csv')

a$y = as.numeric(a$Protein_carbonyl_concentration )
a$x= as.numeric(a$pH_fr_plot)
a$Hion = 10^-a$x

#NAK and pH
ggplot(subset(a), aes(x = x, y =y)) +
  geom_point() +
  theme_test()+
  labs(x = "pH", y = "Protein Carbonyl")+
  geom_smooth(method='gam',formula=y~s(x,k=5),se=FALSE)
##############################################################
###conclusion
#some differences from control but essentially flat after pH drops below 8
###############################################################

############################
#Essentially all the metrics in this paper do not show a 'threshold' per se they are all linear relationships and thus represent a gradient of response (i.e. all pHs are worse than the control and the further you get from control the worse the activity is). I would suggest for responses like this you use the linear relationship and define quantile changes from the control. The threshold calculations I was focusing on were more along the lines of we see no effect between pH 8-7.6, then it gets worse and identifying where that 'break point' occurs.

