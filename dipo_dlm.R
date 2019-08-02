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

theme_set(theme_bw())

##` @param IC    Initial Conditions
##` @param r     Intrinsic growth rate
##` @param Kg    Across-site ('global') mean carrying capacity
##` @param alpha Site random effect
##` @param beta  Slope of precipitation effect on K
##` @param ppt   Precipitation forecast
##` @param Q     Process error (default = 0 for deterministic runs)
##` @param n     Size of Monte Carlo ensemble
##` @param NT    number of forecast timesteps
forecastN <- function(IC, betaIntercept, betaTmin, mintemp,
                      Q = 0, n = Nmc, NT){
  N <- matrix(NA, n, NT)  ## storage
  Nprev <- IC           ## initialize
  for(t in 1:NT){
    mu <- Nprev + betaIntercept + betaTmin * mintemp[,t]
    if (Q==0) {
      N[, t] <- mu
    } else {
      N[, t] <- rpois(n, mu)
    }
    Nprev <- N[,t]                                  ## update IC
  }
  return(N)
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
NT <- 128

# calculate the mean of the driver (tmin) and fitted model parameters
tmin_mean <- t(as.matrix(rep(mean(jags_dat$mintemp, na.rm = TRUE), NT)))

## parameters
params <- as.matrix(out$params)
param_mean <- apply(params, 2, mean, na.rm = TRUE)
## initial conditions
IC <- as.matrix(out$predict)

N_det_mean <- forecastN(IC = mean(IC[,"x[364]"], na.rm = TRUE),
                        betaIntercept = param_mean["betaIntercept"],
                        betaTmin = param_mean["betaTmin"],
                        mintemp = tmin_mean,
                        Q = 0,  ## process error off
                        n = 1,
                        NT = NT)

## Plot run
plot(t(N_det_mean))

#################################################-
## Propogate errors ----
#################################################-
Nmc <- 1000

## Initial conditions
## sample parameter rows from previous analysis
prow <- sample.int(nrow(params), Nmc, replace = TRUE)

N_IC <- forecastN(IC = IC[prow,"x[364]"],
                 betaIntercept = param_mean["betaIntercept"],
                 betaTmin = param_mean["betaTmin"],
                 mintemp = tmin_mean,
                 Q = 0,  ## process error off
                 n = Nmc,
                 NT = NT)

## Plot run
N_IC_ci <- apply(N_IC, 2, quantile, c(0.025, 0.5, 0.975)) %>%
  t() %>%
  as_tibble() %>%
  rename(q2.5 = `2.5%`, q50 = `50%`, q97.5 = `97.5%`) %>%
  rownames_to_column("time") %>%
  mutate(newmoonnumber = as.numeric(time) + 364)

dipo_dat %>%
  full_join(N_IC_ci, by = "newmoonnumber") %>%
  ggplot(aes(x = newmoonnumber, y = DM)) +
  geom_ribbon(aes(ymin = q2.5, ymax = q97.5,
                                  x = newmoonnumber), fill = "blue") +
  geom_line(size = 1.5)





#################################################-
## FC: Known, variable mintemp ----
#################################################-
tmin_obs <- t(as.matrix(jags_dat$mintemp[365:length(jags_dat$mintemp)]))
params <- as.matrix(out$params)
param_mean <- apply(params, 2, mean, na.rm = TRUE)
IC <- as.matrix(out$predict)

N_det_obstmin <- forecastN(IC = mean(IC[,"x[364]"], na.rm = TRUE),
                        betaIntercept = param_mean["betaIntercept"],
                        betaTmin = param_mean["betaTmin"],
                        mintemp = tmin_obs,
                        Q = 0,  ## process error off
                        n = 1,
                        NT = NT)

## Plot run
dipo_dat %>%
  bind_cols(tibble(forecast = c(rep(NA, 364), N_det_obstmin))) %>%
  ggplot(aes(x = newmoonnumber, y = DM)) +
  geom_line(size = 1.5) +
  geom_line(aes(y = forecast), color = 'blue', size = 1.5) 











