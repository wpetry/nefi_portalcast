#################################################-
## Random walk forecast for Dipodymus ----
#################################################-
## Preliminaries ----
#################################################-
library(portalr)
library(tidyverse)
library(rjags)
library(coda)
#################################################-
## Get & subset data to Dipodymus ----
#################################################-
data <- abundance(getwd()) %>%
  as_tibble()

dipo_dat <- data %>%
  select(period, DM) %>%
  right_join(y = tibble(period = min(data$period):max(data$period))) %>%
  arrange(period) %>%
  as.data.frame()

# set aside a training subset
dipo_dat_train <- dipo_dat
dipo_dat_train$DM[floor(3/4 * nrow(dipo_dat)):nrow(dipo_dat)] <- NA

# set up data for JAGS
y <- dipo_dat_train$DM
n <- nrow(dipo_dat_train)
x_ic <- 30
tau_ic <- 100
a_add <- 1
r_add <- 1
period <- dipo_dat$period

dipo_dat_jags <- list(y = y, n = n, x_ic = x_ic,
                      tau_ic = tau_ic, a_add = a_add,
                      r_add = r_add)

#################################################-
## Define model ----
#################################################-
dipo_mod_rw <- "
model{
  #### Data Model
  for(t in 1:n){
    y[t] ~ dpois(x[t])
  }
  
  #### Process Model
  for(t in 2:n){
    x[t] ~ dnorm(x[t-1], tau_add)
  }
  
  #### Priors
  x[1] ~ dnorm(x_ic, tau_ic)
  tau_add ~ dgamma(a_add, r_add)
}
"

#################################################-
## Compile the model ----
#################################################-
# set options
nchains <- 3
# init <- list()  # comment out to just do these from the priors
# for(i in 1:nchain){
#   DM.samp <- sample(dipo_dat_train$DM,
#                     length(dipo_dat_train$DM), replace = TRUE)
#   init[[i]] <- list(tau_add = 1 / var(diff(DM.samp)))
# }

j_mod_rw <- jags.model(file = textConnection(dipo_mod_rw),
                       data = dipo_dat_jags,
                       n.chains = nchains)

#################################################-
## Sample from the posterior ----
#################################################-
j_mod_samp <- coda.samples(model = j_mod_rw,
                           variable.names = c("x", "tau_add"),
                           n.iter = 5000)

# check convergence diagnostics
effectiveSize(j_mod_samp)
# gelman.plot(j_mod_samp) 
gelman.diag(j_mod_samp)  # should be below 1.05

#################################################-
## Visualize the fit ----
#################################################-
time.rng <- c(1,length(period)) ## adjust to zoom in and out
out <- as.matrix(j_mod_samp)
x.cols <- grep("^x", colnames(out))
ci <- apply(out[, x.cols], 2, quantile,
            c(0.025, 0.5, 0.975)) ## model was fit on log scale

plot(period,
     ci[2, ], type = 'n',
     ylim = range(y, na.rm = TRUE),
     ylab = "DM Abundance",
     xlim = period[time.rng])
## adjust x-axis label to be monthly if zoomed
if(diff(time.rng) < 100){ 
  axis.Date(1, at=seq(period[time.rng[1]], period[time.rng[2]],by='month'), format = "%Y-%m")
}
ecoforecastR::ciEnvelope(period,ci[1,],ci[3,],col=ecoforecastR::col.alpha("lightBlue",0.75))
points(dipo_dat$period, dipo_dat$DM, pch = "+", cex = 0.5,
       col = "red")
points(period,y,pch="+",cex=0.5)

