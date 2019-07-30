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

df=abundance()
str(df)

##FILL IN EMPTY PERIODS##
print(df$period)
475-27#max period-min period
448-416#number of missing periods
df2<-data.frame(period=min(df$period):max(df$period))
str(df2)

df3<-merge(df, df2, by="period", all=TRUE)
str(df3)

df4<-subset(df3, select=c("period","PP"))
str(df4)

ggplot(df4, aes(x=period, y=PP))+geom_point()+geom_line()

y = df4$PP
y2=df4$PP
period=df4$period
y[floor(length(y)-0.25*(length(y))):length(y)]<-NA
y


RandomWalk = "
model{
  
#### Process Model
  for(t in 2:n){
    x[t]~dnorm(x[t-1],tau_add)
  }
  
  #### Data Model
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }

  #### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_add ~ dgamma(a_add,r_add)
}
"

data <- list(y=y,n=length(y),x_ic=100 ,tau_ic=100,a_add=1,r_add=1)

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

period.rng = c(1,length(period)) ## adjust to zoom in and out
out <- as.matrix(jags.out)
x.cols <- grep("^x",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 
str(ci)

plot(period,ci[2,],type='n',ylim=range(y,na.rm=TRUE),ylab="Abundance",xlim=period[period.rng])

## adjust x-axis label to be monthly if zoomed

ecoforecastR::ciEnvelope(period,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(period,y,pch="+",cex=1.5)
points(period,y2,pch="+",cex=0.5)

hist(1/sqrt(out[,1]),main=colnames(out)[1])
hist(1/sqrt(out[,2]),main=colnames(out)[2])

plot(out[,1],out[,2],pch=".",xlab=colnames(out)[1],ylab=colnames(out)[2])
cor(out[,1:2])