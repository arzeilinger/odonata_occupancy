#### NIMBLE occupancy model for ZIB GLMM model of potato psyllid occupancy
#### Originally developed by Daniel Turek

# Clear workspace
rm(list = ls())
# Load packages
my.packages <- c("nimble", "coda", "lattice", "akima")
lapply(my.packages, require, character.only = TRUE)

source('R_functions/nimble_definitions.R')

#### Load data set for model fitting
inputData <- readRDS('data/odonata-data-for-nimble.rds')
str(inputData)
inputData$jdn2 <- inputData$jdn^2

# #### Increasing memory limit for R
# memory.limit()
# memory.limit(size = 7000) # Size in Mb
# memory.limit()


#####################################################
#### Define model in BUGS/NIMBLE language

code <- nimbleCode({
    mu_alpha ~ dnorm(0, 0.001)
    sigma_alpha ~ dunif(0, 1000)
    for(j in 1:nsite) { 
        alpha[j] ~ dnorm(mu_alpha, sd = sigma_alpha)  ## site random effect
    }
    for(i in 1:7) {
        beta[i] ~ dnorm(0, 0.001)
    }
    for(i in 1:N) {
        logit(p_occ[i]) <- alpha[siteID[i]] + beta[5]*year[i] + beta[6]*min_temp[i] + beta[7]*total_precip[i]
        logit(p_obs[i]) <- beta[1] + beta[2]*list_length[i] + beta[3]*jdn[i] + beta[4]*jdn2[i]
        y[i] ~ dOccupancy(p_occ[i], p_obs[i])
    }
})

constants <- with(inputData,
                  list(N=N, nsite=nsite, 
                       year=year, 
                       min_temp=min_temp,
                       total_precip=total_precip,
                       list_length=list_length, 
                       jdn=jdn,
                       jdn2=jdn2, 
                       siteID=siteID))

data <- with(inputData, list(y=Amphiagrion.abbreviatum))

inits <- list(mu_alpha=0, sigma_alpha=1, alpha=rep(0,inputData$nsite), beta=rep(0,7))

modelInfo_Amphriagrion <- list(code=code, constants=constants, data=data, inits=inits, name='Amphriagrion_model')


#### Set up model and samplers
Rmodel <- nimbleModel(modelInfo_Amphriagrion$code,
                      modelInfo_Amphriagrion$constants,
                      modelInfo_Amphriagrion$data,
                      modelInfo_Amphriagrion$inits)

Cmodel <- compileNimble(Rmodel)

spec <- configureMCMC(Rmodel)

#### Best configuration of samplers for random effect occupancy model
spec$removeSamplers('beta[1:7]')
spec$addSampler('beta[1:4]', 'RW_block') # detection sub-model sampler
spec$addSampler('beta[5:7]', 'RW_block') # occupancy sub-model sampler
spec$removeSamplers('sigma_alpha')
spec$addSampler('sigma_alpha', 'RW_log_shift', list(shiftNodes='alpha')) # random effect sampler
spec$getSamplers() # Check samplers
spec$addMonitors(c('p_occ')) # add a monitor to get p_occ in output
#spec$addMonitors(c('p_obs')) # add a monitor to get p_obs in output

#### Compile MCMC in R and C++
Rmcmc <- buildMCMC(spec)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)


#### Run MCMC with 150,000 iterations and 50,000 burn-in
niter <- 150000
burnin <- 50000

ti <- Sys.time()
samplesList <- lapply(1:3, mcmcClusterFunction)
tf <- Sys.time()

# The time it took to run MCMC
tf-ti

save(samplesList, file = 'output/MCMC_list_Amphriagrion_pocc.RData')

#### Loading saved MCMC run, sames as list, "samplesList"
#load(file = 'output/MCMC_month_list.RData')


#######################################################################
#### Assessing convergence and summarizing/plotting results
#### Assessing convergence of only covariates and mu.alpha. Memory requirements too great to assess for all p_occ[i]

# Make mcmc.list with only covariates and mu.alpha
mcmcs <- mcmc.list(lapply(1:length(samplesList), function(x) as.mcmc(samplesList[[x]][,1:12])))

## Rhat
coda::gelman.diag(mcmcs, autoburnin = FALSE)
## Effective sample size
effectiveSize(mcmcs)

## Posterior Density Plots
pdf("results/figures/trace_and_posterior_density_plots.pdf")
  plot(mcmcs[[1]], ask = FALSE)
dev.off()


#### Posterior Inferences
#### Mean and 95% Credible intervals
results <- as.data.frame(cbind(apply(samplesList[[1]], 2, mean),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.025)),
                               apply(samplesList[[1]], 2, function(x) quantile(x, 0.975))))
names(results) <- c("mean", "cil", "ciu")
results$params <- row.names(results)
print(results[1:15,]) # Coefficient results


##############################################################################################################
#### Plots

#### Plotting P(occupancy) against covariates
pocc <- results[grep("p_occ", results$params),]
# Load detection data set
detectData <- readRDS("output/potato_psyllid_detection_dataset.rds")
detectData$pocc <- pocc$mean

# Year vs P(occupancy)
tiff("results/figures/occupancy_year_vs_pocc.tif")
  plot(x = detectData$year, y = detectData$pocc)
  lines(smooth.spline(detectData$year, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()
  
# List length vs P(occupancy)
tiff("results/figures/occupancy_list_length_vs_pocc.tif")
  plot(x = detectData$list_length, y = detectData$pocc)
  lines(smooth.spline(detectData$list_length, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()

# List length vs month
tiff("results/figures/occupancy_month_vs_pocc.tif")
  plot(x = detectData$month, y = detectData$pocc)
  lines(smooth.spline(detectData$month, detectData$pocc, nknots = 4, tol = 1e-6), lwd = 2)
dev.off()


# trivariate plots with month and year
zz <- with(detectData, interp(x = year, y = month, z = pocc, duplicate = 'median'))
pdf("results/figures/year-month-occupancy_contourplot_nimble_occupancy.pdf")
  filled.contour(zz, col = topo.colors(32), xlab = "Year", ylab = "Month")
dev.off()



######################################################################################################
#### Extra code

# set.seed(1)
# Cmcmc$run(niter)
# samples1 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# set.seed(2)
# Cmcmc$run(niter)
# samples2 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# set.seed(3)
# Cmcmc$run(niter)
# samples3 <- as.matrix(Cmcmc$mvSamples)[(burnin+1):niter,]
# 
# samplesList <- list(samples1, samples2, samples3)
# 
# save(samples1, samples2, samples3, file = 'output/MCMC_season.RData')


# ##############################################################################
# #### Attempt at cluster
# sfInit(parallel = TRUE, cpus = 2)
# date()
# 
# sfLibrary(nimble)
# sfExport("Cmcmc", "Rmodel", "niter", "burnin", "mcmcClusterFunction")
# sfExport("spec")
# sfExport("Cmodel")
# sfExport("Rmcmc")
# sfLapply(1:2, mcmcClusterFunction)
# 
# date()
# sfStop() # Close the cluster
# #################################################################################
