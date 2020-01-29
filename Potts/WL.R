##### Wang-Landau algorithm for 2D Ising model
rm(list = ls())
library(LaplacesDemon)
library(Rcpp)
### The C code contains two functions
### 1. Calculate the energy of a lattice
### 2. One-step of MCMC sweep, which cantains L^2 steps
replicate <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(replicate)
sourceCpp("/n/home09/cdai/AccWL/revision/code/Potts/WL.cpp")
load("/n/home09/cdai/AccWL/revision/Potts_error/truth/L_80_truth.RData")
logdensity_truth <- data_save$logdensity

### Algorithmic settings
num_mc <- 2*(10^6)
L <- 80
flatness_criterion <- 0.8
N <- L^2
logf <- 0.05
min_logf <- 10^(-8)
energy_range <- seq(from = 0, to = L^2*2, by = 1)
empty_index <- c(2, 3, 4, 6)

### Initialization
### Initialize the logdensity. The second and the last to second are empty.
logdensity <- rep(1, L^2*2 + 1)  
logdensity[empty_index] <- 0
### Initialize the ground state. All spins are 1.
lattice <- matrix(1, nrow = L, ncol = L)
### Initialize the energy level, index by the C language.
energy_level <- 0
### Initialize the histogram.
Hist <- rep(0, L^2*2 + 1)
Hist_index <- 0
### Keep track of error
error <- rep(0, num_mc)
### Keep track of the free energy
free_energy <- rep(0, num_mc)

### Wang-Landau algorithm
start_time <- Sys.time()
for(iter in 1:num_mc){
  ### Metropolis move
  MC_result <- MC_sweep(lattice, logdensity, Hist, energy_level, logf)
  lattice <- MC_result$lattice
  logdensity <- MC_result$logdensity
  energy_level <- MC_result$energy_level
  Hist <- MC_result$Hist
  
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
  if(iter%%1000 == 0) print(c(iter, logf, error[iter]))
}
end_time <- Sys.time()
running_time <- end_time - start_time

data_save <- list(error = error, running_time = running_time)
save(data_save, file = paste("/n/home09/cdai/AccWL/revision/Potts_error/result/WL0.05/replicate_", replicate, ".RData", sep = ""))

