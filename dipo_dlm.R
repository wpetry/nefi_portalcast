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

#################################################-
## Merge abundance and weather data, prepare for JAGS ----
#################################################-
merged_dat <- left_join(dipo_dat_train, wx,
                        by = "newmoonnumber")

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

dlm_mintemp <- "
model{
#### Priors
  x[1] ~ dnorm(x_ic,tau_ic)
  tau_add ~ dgamma(a_add,r_add)
  
  #### Fixed Effects
  betaX ~dnorm(0,0.001)
  betaIntercept~dnorm(0,0.001)
  betamintemp~dnorm(0,0.001)
  for(j in 1: 2 ){
    muXf[j] ~ dnorm(0,0.001)
    tauXf[j] ~ dgamma(0.01,0.01)
  }
  
  #### Data Model
  for(t in 1:n){
    OBS[t] ~ dpois(x[t])
    Xf[t,1] ~ dnorm(muXf[1],tauXf[1])
    Xf[t,2] ~ dnorm(muXf[2],tauXf[2])
  }
  
  #### Process Model
  for(t in 2:n){
    mu[t] <- x[t-1]  + betaX*x[t-1] + betaIntercept*Xf[t,1] + betamintemp*Xf[t,2]
    x[t]~dnorm(mu[t],tau_add)
  }
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


#################################################-
## DLM fit diagnostics ----
#################################################-
params <- window(ef.out$params,start=1000) ## remove burn-in
plot(params)

