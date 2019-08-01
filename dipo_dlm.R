#################################################-
## Dynamic Linear Model for Dipodymus (abiotic) ----
#################################################-
## Preliminaries ----
#################################################-
library(portalr)
library(ecoforecastR)
library(tidyverse)
library(rjags)
library(coda)

##` @param IC    Initial Conditions
##` @param r     Intrinsic growth rate
##` @param Kg    Across-site ('global') mean carrying capacity
##` @param alpha Site random effect
##` @param beta  Slope of precipitation effect on K
##` @param ppt   Precipitation forecast
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
##` @param NT    number of forecast timesteps
forecastN <- function(IC,r,Kg,alpha,beta,ppt,Q=0,n=Nmc,NT){
  N <- matrix(NA,n,NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    mu <- Nprev + betaIntercept + betaTmin * mintemp[,t]
    
    # if(Q=0), mu
    # if(Q!=0), rpois(n, mu)
    
    x[t]~dnorm(mu[t],tau_add)
    
    
    K = pmax(1,Kg + alpha + beta*log(ppt[,t]/800))  ## calculate carrying capacity
    mu = log(pmax(1,Nprev + r*Nprev*(1-Nprev/K)))   ## calculate mean
    N[,t] <- rlnorm(n,mu,Q)                         ## predict next step
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
}

#################################################-
## Get & subset data to Dipodymus ----
#################################################-
data <- abundance(getwd(), time = "newmoon",
                  clean = FALSE) %>%
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

ggplot(wx, aes(y = mintemp, x = newmoonnumber))+
  geom_line()

#################################################-
## Merge abundance and weather data, prepare for JAGS ----
#################################################-
merged_dat <- left_join(dipo_dat_train, wx,
                        by = "newmoonnumber") %>%
  filter(newmoonnumber >= 32)

#################################################-
## Fit DLM with mintemp as the predictor ----
#################################################-
# The strategy here is to use MD's automatic DLM code,
# then tailor it to our problem. We can get this autocode
# as follows:
# ef.out <- fit_dlm(model = list(obs = "DM",
#                                fixed = "~ 1 + X + mintemp"),
#                   merged_dat)
# write.table(ef.out$model,file = "test")

dlm_mintemp <-
"model{
  
  #### Process Model
  for(t in 2:n){
    mu[t] <- x[t-1] + betaIntercept + betaTmin*mintemp[t-1]
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
# make data for JAGS & add prior hyperparameters
jags_dat <- merged_dat %>%
  select(y = DM, mintemp) %>%
  as.list()
jags_dat$x_ic <- 30
jags_dat$tau_ic <- 100
jags_dat$a_add <- 1
jags_dat$r_add <- 1
jags_dat$n <- length(jags_dat[["y"]])

#################################################-
## Run the JAGS model ----
#################################################-
jags_model <- jags.model(textConnection(dlm_mintemp),
                         data = jags_dat,
                         n.chains = 3)

dlm_mintemp_samps <- coda.samples(jags_model,
                                  variable.names = c("x", "tau_add",
                                                     "betaIntercept",
                                                     "betaTmin"),
                                  n.iter = 5000)
#################################################-
## DLM fit diagnostics ----
#################################################-
effectiveSize(dlm_mintemp_samps)
gelman.plot(dlm_mintemp_samps[, c("tau_add",
                                 "betaIntercept",
                                 "betaTmin")])
plot(dlm_mintemp_samps[, c("tau_add",
                           "betaIntercept",
                           "betaTmin")])

# remove burn in of 1000 samples
dlm_mintemp_samps <- window(dlm_mintemp_samps, start = 1000)

# add predictions
out <- list(params = NULL, predict = NULL, model = jags_model, 
           data = jags_dat)
mfit <- as.matrix(dlm_mintemp_samps, chains = TRUE)
pred.cols <- union(grep("x[", colnames(mfit), fixed = TRUE), 
                  grep("mu[", colnames(mfit), fixed = TRUE))
chain.col <- which(colnames(mfit) == "CHAIN")
out$predict <- mat2mcmc.list(mfit[, c(chain.col, pred.cols)])
out$params <- mat2mcmc.list(mfit[, -pred.cols])

#################################################-
## Make a deterministic forecast ----
#################################################-
# calculate the mean of the driver (tmin) and fitted model parameters
tmin_mean <- mean(jags_dat$mintemp, na.rm = TRUE)
## parameters
params <- as.matrix(out$params)
param.mean <- apply(params, 2, mean)
## initial conditions
IC <- as.matrix(out$predict)

N.det <- forecastN(IC = mean(IC[,"X[364]"]),
                   betaIntercept = param.mean["betaIntercept"],
                   betaTmin = param.mean["betaTmin"],
                   tau_add = param.mean["betaIntercept"],
                   mintemp = tmin_mean,
                   Q=0,  ## process error off
                   n=1)






## calculate mean of all inputs
ppt.mean <- matrix(apply(ppt_ensemble,2,mean),1,NT) ## driver
## parameters
params <- as.matrix(out$params)
param.mean <- apply(params,2,mean)
## initial conditions
IC <- as.matrix(out$predict)

N.det <- forecastN(IC=mean(IC[,"N[6,30]"]),
                   r=param.mean["r_global"],
                   Kg=param.mean["K_global"],
                   alpha=param.mean["alpha_site[6]"],
                   beta=param.mean["beta"],
                   ppt=ppt.mean,
                   Q=0,  ## process error off
                   n=1)

## Plot run
plot.run()
lines(time2,N.det,col="purple",lwd=3)












