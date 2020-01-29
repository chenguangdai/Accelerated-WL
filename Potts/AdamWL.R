##### Wang-Landau algorithm for 3D Potts model
rm(list = ls())
library(LaplacesDemon)
library(Rcpp)
### The C code contains two functions
### 1. Calculate the energy of a lattice
### 2. One-step of MCMC sweep, which cantains L^2 steps
### replication index
replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)
load("/n/home09/cdai/AccWL/revision/Potts_error/truth/L_80_truth.RData")
logdensity_truth <- data_save$logdensity
sourceCpp("/n/home09/cdai/AccWL/revision/code/Potts/AdamPotts.cpp")

### Algorithmic settings
num_mc <- 2*(10^6)
L <- 80
flatness_criterion <- 0.8
N <- L^2
logf <- 0.05
energy_range <- seq(from = 0, to = L^2*2, by = 1)
empty_index <- c(2, 3, 4, 6)
beta <- 0.9

### Initialization
### Initialize the logdensity. The second and the last to second are empty.
logdensity <- rep(0, L^2*2 + 1)  
### Initialize the histogram.
Hist <- rep(0, L^2*2 + 1)
Hist_index <- 0
### Keep track of error
error <- rep(0, num_mc)
### Initialize the momentum
momentum <- rep(0, L^2*2 + 1)
lattice <- matrix(1, nrow = L, ncol = L)
energy_level <- 0
### Last update
last_update <- rep(0, L^2*2 + 1)


### Wang-Landau algorithm
start_time <- Sys.time()
for(iter in 1:num_mc){
  ### Metropolis move
  MC_result <- MC_sweep(lattice, logdensity, Hist, energy_level, logf, beta, last_update, momentum, iter)
  lattice <- MC_result$lattice
  logdensity <- MC_result$logdensity
  energy_level <- MC_result$energy_level
  Hist <- MC_result$Hist
  momentum <- MC_result$momentum
  last_update <- MC_result$last_update
  
  ### Check the flatness criterion
  if(iter%%1000 == 0 && Hist_index == 0){
    Hist_min <- min(Hist[-empty_index])
    if(Hist_min > 0){
      ### Adjust the modification factor and reset the histogram
      logf <- logf / 2
      Hist <- rep(0, L^2*2 + 1)
    }
  }
  if(iter >= 1000 & logf < 1/iter){Hist_index <- 1}
  if(Hist_index == 1){logf <- 1/iter}
  
  ### Normalize the density
  maxlw <- max(logdensity[-empty_index])
  logconst <- log(sum(exp(logdensity[-empty_index] - maxlw))) + maxlw
  logdensity[-empty_index] <- logdensity[-empty_index] - logconst
  
  ### Calculate the error
  error[iter] <- sum(abs(1 - logdensity[-empty_index]/logdensity_truth[-empty_index]))/(N - 1)
  ### Calculate the free energy.
  if(iter%%1000 == 0) print(c(iter, logf, error[iter]))
}
end_time <- Sys.time()
running_time <- end_time - start_time

data_save <- list(error = error, running_time = running_time)
save(data_save, file = paste("/n/home09/cdai/AccWL/revision/Potts_error/result/AdamWL/replicate_", replicate, ".RData", sep = ""))


