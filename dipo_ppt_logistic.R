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
  beta ~ dnorm(0, 0.000001)    ## slope of K response to precip
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
      K[s, t]  <- max(1, K_global + beta * log(precip[t] / 23.70712))  

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
                           variable.names = c("No", "K_global",
                                              "r_global",
                                              "beta"),
                           n.iter = 5000)







