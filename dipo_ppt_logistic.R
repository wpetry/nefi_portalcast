#################################################-
## Logistic growth with precipitation-dependent K model for Dipodymus ----
#################################################-
## Preliminaries ----
#################################################-
library(portalr)
library(ecoforecastR)
library(tidyverse)
library(rjags)
library(coda)

theme_set(theme_bw())

##` @param IC    Initial Conditions
##` @param r     Intrinsic growth rate
##` @param Kg    Across-site ('global') mean carrying capacity
##` @param alpha Site random effect
##` @param beta  Slope of precipitation effect on K
##` @param ppt   Precipitation forecast
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
forecastN <- function(IC,r,Kg,alpha,beta,ppt,Q=0,n=Nmc){
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    K = pmax(1,Kg + alpha + beta*log(ppt[,t]/23.70712))  ## calculate carrying capacity
    mu = log(pmax(1,Nprev + r*Nprev*(1-Nprev/K)))   ## calculate mean
    N[,t] <- rlnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

plot.run <- function(){
  sel = seq(s,ncol(ci),by=NS)
  plot(time,time,type='n',ylim=ylim,ylab="N")
  ecoforecastR::ciEnvelope(time1,ci[1,sel],ci[3,sel],col=col.alpha("lightBlue",0.6))
  lines(time1,ci[2,sel],col="blue")
  points(time1,No[s,])
}

#################################################-
## Get & subset data to Dipodymus ----
#################################################-
data <- abundance(level = "plot",
                  time = "newmoon",
                  clean = FALSE,
                  plots = "Longterm") %>% 
  group_by(newmoonnumber, treatment) %>% 
  summarize(DM = sum(DM), PP = sum(PP)) %>% 
  data.frame() %>%
  filter(treatment == "control") %>%
  as_tibble()

ggplot(data, aes(x = newmoonnumber, y = DM))+
  geom_line()

dipo_dat <- data %>%
  select(newmoonnumber, DM) %>%
  right_join(y = tibble(newmoonnumber = min(data$newmoonnumber):max(data$newmoonnumber))) %>%
  arrange(newmoonnumber) %>%
  as.data.frame()

# set aside a training subset
dipo_dat_train <- dipo_dat
dipo_dat_train$DM[floor(3/4 * nrow(dipo_dat)):nrow(dipo_dat)] <- NA

#################################################-
## Get weather data for Portal ----
#################################################-
wx <- weather(level = "newmoon", fill = TRUE)

ggplot(wx, aes(y = precipitation, x = newmoonnumber))+
  geom_line()

#################################################-
## Merge abundance and weather data, prepare for JAGS ----
#################################################-
merged_dat <- left_join(dipo_dat_train, wx,
                        by = "newmoonnumber") %>%
  filter(newmoonnumber >= 32)

#################################################-
## Fit ppt-dependent K logistic growth model ----
#################################################-
logisticRE <- "
model{

  ## priors
  r_global ~ dnorm(0, 0.1)     ## across-site mean growth rate
  K_global ~ dlnorm(6, 0.01)   ## across-site mean carrying capacity
  beta_1 ~ dnorm(0, 0.000001)    ## slope of K response to precip - current mo
  beta_2 ~ dnorm(0, 0.000001)    ## slope of K response to precip - previous mo
  R ~ dgamma(0.01, 0.00000001) ## Observation error precision
  Q ~ dgamma(0.01, 0.00000001) ## Process errror precision 

  ## initial conditions, s = site
  for(s in 1:NS){
    lN[s, 1] ~ dnorm(6, 0.001)           ## prior on IC, log scale
    N[s, 1] <- exp(lN[s, 1])             ## IC, linear scale
  }

  ## process model, t = time, s = site
  for(t in 2:NT){
    for(s in 1:NS){

      ## K is a linear model with a site random effect and fixed effect on log(precip), centered on the long-run mean precipitation
      K[s, t]  <- max(1, K_global + beta_1 * log(precip[t] / 23.70712) + beta_2 * log(precip[t - 1] / 23.70712))  

      ## standard logistic growth process model, logged     
      mu[s, t] <- log(max(1, N[s, t - 1] + r_global * N[s, t - 1] * (1 - N[s, t - 1] / K[s, t])))

      ## process error
      lN[s, t] ~ dnorm(mu[s, t], Q)
      N[s, t] <- exp(lN[s, t])
    }
  }
  ## observation model
  for(t in 1:NT){
    for(s in 1:NS){
      No[s, t] ~ dlnorm(lN[s, t], R)
    }
  }
}
"

jags_dat <- merged_dat %>%
  select(No = DM, precip = precipitation) %>%
  mutate(alpha_site = 1) %>%
  as.list()
jags_dat$No <- t(as.matrix(jags_dat$No))
jags_dat$NT <- length(jags_dat[["No"]])
jags_dat$NS <- 1

# compile the JAGS model
jags_model <- jags.model(textConnection(logisticRE),
                         data = jags_dat,
                         n.chains = 3)

# sample that shit!!!
jags_samps <- coda.samples(jags_model,
                           variable.names = c("No","N", "K_global",
                                              "r_global",
                                              "beta_1", "beta_2"),
                           n.iter = 10000)


out <- list(params = NULL, predict = NULL, model = jags_model, 
            data = jags_samps)
mfit <- as.matrix(jags_samps, chains = TRUE)
pred.cols <- union(grep("N[", colnames(mfit), fixed = TRUE), 
                   grep("mu[", colnames(mfit), fixed = TRUE))

chain.col <- which(colnames(mfit) == "CHAIN")

out$predict <- mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
out$params <- mat2mcmc.list(mfit[, -pred.cols])

Nmc <- 1000

plot(jags_samps[,c("K_global","beta_1","beta_2")])

#################################################-
## Make credible intervals & prediction intervals ----
#################################################-
# try again
newmoon.rng = c(1,length(merged_dat$newmoonnumber)) ## adjust to zoom in and out
out <- as.matrix(jags_samps)
x.cols <- grep("^No",colnames(out)) ## grab all columns that start with the letter x
ci <- apply(out[,x.cols],2,quantile,c(0.025,0.5,0.975)) 
str(ci)

plot(merged_dat$newmoonnumber,ci[2,],type='n',ylim=c(min(ci[1,]), max(ci[3,])),ylab="Abundance",xlim=c(200,500), log="y", xlab = "Time (months)")

ecoforecastR::ciEnvelope(merged_dat$newmoonnumber,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(merged_dat$newmoonnumber,merged_dat$DM,pch="+",cex=1.5)
points(dipo_dat$newmoonnumber,dipo_dat$DM,pch="+",cex=0.5)
lines(merged_dat$newmoonnumber, ci[2,])

#################################################-
## Paritition errors [NOT RUN] ----
#################################################-
## Initial conditions
## sample parameter rows from previous analysis
# prow <- sample.int(nrow(as.matrix(out$params)), Nmc, replace = TRUE)
# 
# ## initial conditions
# IC <- as.matrix(out$predict)
# 
# N_IC <- forecastN(IC = IC[prow,"x[364]"],
#                   r=params[prow,"r_global"],  ## sample parameters
#                   Kg=params[prow,"K_global"],
#                   alpha=0,
#                   beta=params[prow,"beta"],
#                   ppt=precip,
#                   Q=0,  ## process error off
#                   n=Nmc)
# 
# ## Plot run
# plot.run()
# lines(time2,N.det,col="purple",lwd=3)
