############################################
############################################
#VERSION 3

library(ggplot2)
library(devtools)
library(rjags)
library(portalr)
library(dplyr)
library(remotes)
library(imputeTS)
library(ggpubr)

df.a = abundance(level="plot",
                 time="newmoon",
                 clean=FALSE,
                 plots = "Longterm") %>% 
  group_by(newmoonnumber, treatment) %>% 
  filter(treatment=="control")%>%
  summarize(DM = sum(DM), PP = sum(PP)) %>% 
  data.frame()

temp=weather(level="newmoon", fill=TRUE)
temp=subset(temp, !newmoonnumber=="NA")
temp=subset(temp, !newmoonnumber>519)#max moon for rodents
df.a=subset(df.a, !newmoonnumber<200)#truncate to ignore early 0's
temp=subset(temp, !newmoonnumber<200)#truncate to ignore early 0's

str(df.a)
str(temp)

df<-merge(df.a, temp, by="newmoonnumber", all=TRUE)

head(df)

df2<-subset(df, select=c("newmoonnumber","PP","DM","mintemp"))
df2$DM<-round(na_interpolation(df2$DM), digits=0)
head(df2)

##FILL IN EMPTY PERIODS##

df3<-data.frame(newmoonnumber=min(df$newmoonnumber):max(df$newmoonnumber))
str(df3)

df4<-merge(df2, df3, by="newmoonnumber", all=TRUE) 
str(df4)

#dev.off()
ggplot(df4, aes(x=newmoonnumber, y=PP))+geom_point()+geom_line()

tp<-ggplot(data=df4, aes(x=newmoonnumber, y=mintemp))+ geom_line()+xlab("")+ylab("Min. temp. (C)")
pp<-ggplot(data=df4, aes(x=newmoonnumber, y=PP))+geom_line()+xlab("")+ylab("C. penicillatus")
dm<-ggplot(data=df4, aes(x=newmoonnumber, y=DM))+geom_line()+xlab("Time")+ylab("D. merriami")

ggarrange(tp, pp, dm + rremove("x.text"), 
          nrow = 3)

y = df4$PP
y2=df4$PP
newmoonnumber=df4$newmoonnumber
y[floor(length(y)-0.25*(length(y))):length(y)]<-NA
y
temp=df4$mintemp
dm=df4$DM

RandomWalk = "
model{
  
#### Process Model
  for(t in 2:n){
    mu[t] <- max(0,x[t-1] + betaIntercept + betaDTmin*(temp[t]-temp[t-1]) + betaDM*dm[t-1])
    x[t]~dnorm(mu[t],tau_add) T(0,)
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_add ~ dgamma(a_add,r_add)
  betaIntercept~dnorm(0,0.001)
  betaDTmin~dnorm(0,0.001)
  betaDM~dnorm(0,0.001)
}
"

data <- list(y=y,n=length(y),temp=temp, x_ic=2 ,tau_ic=1/10000,a_add=1,r_add=1, dm=dm) 

nchain = 3
init <- list()
for(i in 1:nchain){
  y.samp = sample(y,length(y),replace=TRUE)
  init[[i]] <- list(tau_add=1/var(diff(y.samp)))
}

j.model   <- jags.model (file = textConnection(RandomWalk),
                         data = data,
                         inits = init,
                         n.chains = 3)

## burn-in
jags.out   <- coda.samples (model = j.model,
                            variable.names = c("tau_add"),
                            n.iter = 10000)
plot(jags.out)

jags.out   <- coda.samples (model = j.model,
                            variable.names = c("x","tau_add", "betaIntercept", "betaDTmin", "betaDM"),
                            n.iter = 10000)


summary(jags.out)

par(mar=c(1,1,1,1))
plot(jags.out[,c("tau_add", "betaIntercept", "betaDTmin","betaDM")])

newmoon.rng = c(1,length(newmoonnumber)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 
str(ci)

plot(newmoonnumber,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Abundance",xlim=newmoonnumber[newmoon.rng])

## adjust x-axis label to be monthly if zoomed

ecoforecastR::ciEnvelope(newmoonnumber,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(newmoonnumber,y,pch="+",cex=1.5)
points(newmoonnumber,y2,pch="+",cex=0.5)
lines(newmoonnumber, ci[2,])

hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])

plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])