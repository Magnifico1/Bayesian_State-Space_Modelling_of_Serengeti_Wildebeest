# 0. Data and Libraries and general wrangling ------------------------------------
## 0.1 Seed ----------------------------------------------------------------------
set.seed(141122) # today's date

## 0.2 Required libraries --------------------------------------------------------
library(tidyverse)
library(statsecol)
library(zoo)
library(jagsUI)
library(MCMCvis)
library(formattable)

# 0.3 loading and wrangling the data ---------------------------------------------
# load data - remember to set the working directory
data("wildebeest")
# Identify the observation data without NAs
non_na_indices <- which(!is.na(wildebeest$Nhat)) 
# Back filling the NAs
wildebeest.locf <- na.locf(wildebeest, fromLast = T)

# 1. Specifying the model --------------------------------------------------------
## 1. BUGS language --------------------------------------------------------------
sink("wildebeest.txt")
cat("
model{

  # Priors and constraints
  N1 ~ dunif(0, 2)
  N[1] <- N1

  b0 ~ dnorm(0, 0.01) # Prior for the Intercept
  b1 ~ dnorm(0, 0.01) # Prior for the slope

  sig.n ~ dunif(0, 1)     # Standard deviation for the population
  tau.n <- pow(sig.n, -2) # Precision for the state process

  # Likelihood - State process
  for(t in 1:(nyrs-1)) {
    log.r[t] <- b0 + b1*rain[t] # Linear model for the logarithm of the growth rate
    log(r[t]) <- log.r[t]       # Link function for the parameter
    N[t+1] ~ dnorm((r[t] * N[t]) - catch[t], tau.n)
  }

  # Likelihood - Observation process
  for (t in non_na_indices) {   # Excluded observation years with NAs
     y[t] ~ dnorm(N[t], pow(sig.y[t], -2))

  }
}
", fill = TRUE)
sink()

# 2. Model parameters -----------------------------------------------------------
# 2.1. Fitting the model --------------------------------------------------------
projection <- 5 # No. years to project

wildebeest.model <- list(
  rain = c(wildebeest$rain, rep(mean(wildebeest$rain), projection)),  # Vector of rain quantity observed each year + mean of rainfall in projected years
  catch = c(wildebeest$Catch, rep(0.04, projection)),                 # Catch - the amount harvested (or illegal poaching) + mean of catch in projected years
  y = c(wildebeest$Nhat, rep(NA, projection)),                        # Vector of population observations in millions
  nyrs = nrow(wildebeest) + projection,                               # Number of years to include in the simulation
  sig.y = c(wildebeest$sehat, sample(wildebeest$sehat, projection)),  # Standard deviation for the observed populations
  non_na_indices = non_na_indices                                     # Index of the actual observations (excluded NAs obervations)
)

# 2.2. Parameter initialisation -------------------------------------------------
initial_values <- function() {
  list(
    N = c(wildebeest.locf$Nhat, rep(NA, projection)),
    sig.n = runif(1, 0, 1),
    b0 = runif(1, -2, 2),
    b1 = runif(1, -2, 2)
  )
}

# 2.3. Parameters to log --------------------------------------------------------
monitored_parameters <- c("b0", "b1", "sig.n", "N")

# 2.4 MCMC Initialistion --------------------------------------------------------
nc <- 3            # No. of chains that we want to run
nb <- 10000        # No. of iterations for burn-in.
ni <- 100000 + nb  # No. of iterations
nt <- 1            # No. of thinning samples (try changing this)

# 3. The Model ------------------------------------------------------------------
# 3.1. Running the model --------------------------------------------------------
wildebeest.out <- jags(
  data = wildebeest.model,
  inits = initial_values,
  parameters.to.save = monitored_parameters,
  model.file = "wildebeest.txt",
  n.chains = nc,
  n.iter = ni,
  n.burnin = nb,
  n.thin = nt
)

# 3.2. Wrangling the output into a nice format for plotting --------------------
wildebeest.projection <- data.frame(
  Year = c(wildebeest$year, 1990:1994),
  Mean = wildebeest.out$mean$N,
  Lower = wildebeest.out$q2.5$N,
  Upper = wildebeest.out$q97.5$N,
  Obs = c(wildebeest$Nhat, rep(NA, projection)),
  proj = c(rep(NA, 30), wildebeest.out$mean$N[31:35])
)

# 4. MCMC trace and parameter summary ------------------------------------------
# 4.1. Trace plots and density plots for the monitored parameters --------------
MCMCtrace(
  wildebeest.out,                     # the fitted model
  params = monitored_parameters[1:3], # out parameters of interest
  type = "trace",
  iter = ni,                          # plot all iterations
  pdf = FALSE,                        # don't write to a PDF
  ind = FALSE                         # chain specific densities
)

MCMCtrace(
  wildebeest.out,                     # the fitted model
  params = monitored_parameters[1:3], # out parameters of interest
  type = "density",
  iter = ni,                          # plot all iterations
  pdf = FALSE,                        # don't write to a PDF
  ind = FALSE                         # chain specific densities
)

# 4.2. MCMC Parameter Summary --------------------------------------------------
MCMC_Param_Summary <- MCMCsummary(wildebeest.out,
                                  params = monitored_parameters[1:3]
)

# 5. Plots and tables ----------------------------------------------------------
## 5.1 Population with 95% bounds and projection--------------------------------
ggplot(data = wildebeest.projection) +
  geom_ribbon(aes(x = Year, y = Mean, ymin = Lower, ymax = Upper),
              fill = "grey", alpha = 0.25
  ) +
  geom_line(aes(x = Year, y = Mean), size = 1.2, color = "blue") +
  geom_line(aes(x = Year, y = proj), size = 1.2, color = "red") +
  
  geom_point(aes(x = Year, y = Obs), size = 1.2) +
  geom_vline(xintercept = 1990, col = "black", lty = 3) +  # line at 1990 to show where prediction begins
  geom_vline(xintercept = 1977, col = "black", lty = 3) +  # line at 1977 to show where catch begins
  theme_bw() +
  labs(
    title = "Population abundance (1969-89), and population projection (1990-94) with 95% CI bounds",
    y = "Mean Wildebeest Abundance (millions)",
    x = "Observation Year"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
# TO DO: change colour of the projected part, 1990:1994

## 5.2 Ensemble plot ------------------------------------------------------------
Years <- c(wildebeest$year, 1990:1994)
plot(Years, c(wildebeest.locf$Nhat, rep(NA, 5)),
     ylim = c(0, 2), xlim = c(1960, 1994), las = 1,
     ylab = "Abundance", xlab = "Year", pch = "", main = "Population simulated trajectories and projection"
)
matlines(
  x = Years, y = t(wildebeest.out$sims.list$N[sample(1:5000, 300, replace = TRUE), ]),
  lty = 1, col = adjustcolor("blue", 0.05)
)
lines(Years, apply(wildebeest.out$sims.list$N, 2, mean), col = "blue", lwd = 2)
lines(1990:1994, apply(wildebeest.out$sims.list$N, 2, mean)[31:35], col = "red", lwd = 2, lty = 2)
points(Years, c(wildebeest$Nhat, rep(NA, 5)), col = "black", pch = 16)
# This is based on Chris's plots - convert to GGPLOT???

## 5.3 Priors table ---------------------------------------------------------------
library(formattable)

Prior_df <- data.frame(
  Prior = c("N1 ~ uniform(0, 2)", "beta0 ~ normal(0, 100)", "beta1 ~ normal(0, 100)",
            "sig.n ~ uniform(0, 1)"))

formattable(Prior_df, 
            align =c("c")
)

## 5.4 Initial parameters ---------------------------------------------------------

## 5.5 Posterior summary ----------------------------------------------------------
colnames(MCMC_Param_Summary) <- c("Mean", "Sd", "Lower", "Medium", "Upper",
                                  "Rhat", "n.eff")
MCMC_Param_Summary <- rownames_to_column(MCMC_Param_Summary, var = "Parameter")
MCMC_Param_Summary[,2:6] <- MCMC_Param_Summary[,2:6] %>%
  mutate(across(where(is.numeric), round, digits = 3))

formattable(MCMC_Param_Summary, 
            align =c("l","c","c","c","c", "c", "c", "c", "r"), 
            list('Parameter' = formatter(
              "span", style = ~ style(color = "grey",font.weight = "bold")) 
            ))

## 5.6 Projected years summary ----------------------------------------------------
Proj_Summary <- data.frame(
  Mean = wildebeest.out$mean$N[31:35],
  Sd = wildebeest.out$sd$N[31:35],
  Lower = wildebeest.out$q2.5$N[31:35],
  Medium = wildebeest.out$q50$N[31:35],
  Upper= wildebeest.out$q97.5$N[31:35],
  row.names = 1990:1994)

Proj_Summary <- rownames_to_column(Proj_Summary, var = "Year")
Proj_Summary[,2:6] <- Proj_Summary[,2:6] %>%
  mutate(across(where(is.numeric), round, digits = 3))

formattable(Proj_Summary, 
            align =c("l","c","c","c","c", "c", "c", "c", "r"), 
            list('Year' = formatter(
              "span", style = ~ style(color = "grey",font.weight = "bold")) 
            ))

