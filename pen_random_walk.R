############################################
############################################
#VERSION 3

install.packages('portalr')
install.packages('remotes')
remotes::install_github('ha0ye/rEDM')
install.packages("rEDM")
install.packages("devtools")
install.packages("rjags")
install.packages("Rcpp")
install.packages("ggplot2")
install.packages("usethis")
library(ggplot2)
library(Rcpp)
library(devtools)
library(rjags)
library(rEDM)
library(portalr)
library(dplyr)
library(remotes)

df.a=abundance(time="newmoon", clean=FALSE)
temp=weather(level="newmoon", fill=TRUE)
temp=subset(temp, !newmoonnumber=="NA")
temp=subset(temp, !newmoonnumber>519)
df.a=subset(df.a, !newmoonnumber<32)

str(df.a)
str(temp)

max(df.a$newmoonnumber)
max(temp$newmoonnumber)
min(df.a$newmoonnumber)
min(temp$newmoonnumber)
519-32

df<-merge(df.a, temp, by="newmoonnumber", all=TRUE)
str(df)
head(df)
df2<-subset(df, select=c("newmoonnumber","PP","mintemp"))
str(df2)

##FILL IN EMPTY PERIODS##

df3<-data.frame(newmoonnumber=min(df$newmoonnumber):max(df$newmoonnumber))
str(df3)

df4<-merge(df2, df3, by="newmoonnumber", all=TRUE) 
str(df4)

dev.off()
ggplot(df4, aes(x=newmoonnumber, y=PP))+geom_point()+geom_line()

y = df4$PP
y2=df4$PP
newmoonnumber=df4$newmoonnumber
y[floor(length(y)-0.25*(length(y))):length(y)]<-NA
y
temp=df4$mintemp


RandomWalk = "
model{
  
#### Process Model
  for(t in 2:n){
    mu[t] <- x[t-1] + betaIntercept + betaTmin*temp[t-1]
    x[t]~dnorm(mu[t],tau_add)
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_add ~ dgamma(a_add,r_add)
  betaIntercept~dnorm(0,0.001)
  betaTmin~dnorm(0,0.001)  
}
"

data <- list(y=y,n=length(y),temp=temp, x_ic=100 ,tau_ic=100,a_add=1,r_add=1)

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
                            variable.names = c("x","tau_add"),
                            n.iter = 10000)

newmoon.rng = c(1,length(newmoonnumber)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 
str(ci)

par(mar=c(1,1,1,1))
plot(newmoonnumber,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Abundance",xlim=newmoonnumber[newmoon.rng])

## adjust x-axis label to be monthly if zoomed

ecoforecastR::ciEnvelope(newmoonnumber,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(newmoonnumber,y,pch="+",cex=1.5)
points(newmoonnumber,y2,pch="+",cex=0.5)

hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])

plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])