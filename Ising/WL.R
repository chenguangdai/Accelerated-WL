##### The Wang-Landau algorithm for the 2D Ising model
rm(list = ls())
library(LaplacesDemon)
library(Rcpp)
### The C code contains two functions
### 1. Calculate the energy of a lattice
### 2. One-step of MC sweep, which cantains LxL MC trial moves
sourceCpp("WL.cpp")

### load the true logdensity for L = 80
load("logdensityL80.RData")
logdensity_truth <- logdensity

### Algorithmic settings
num_mc <- 200000
L <- 80
flatness_criterion <- 0.8
N <- L^2
logf <- 0.05
min_logf <- 10^(-8)
energy_range <- seq(from = -2*L^2, to = 2*L^2, by = 4)
empty_index <- c(2, L^2)  ## non-existent energy level indexes

### Initialization
### Initialize the lattice
lattice <- matrix(2*rbern(N, prob = 0.5) - 1, nrow = L, ncol = L)
energy_level <- (get_energy(lattice) + 2*L^2)/4
### Initialize the logdensity. 
logdensity <- rep(0, L^2 + 1)  
### Initialize the histogram.
Hist <- rep(0, L^2 + 1)
Hist_index <- 0
### Keep track of the estimation error epsilon(t) defined as in Equation (24)
error <- rep(0, num_mc)

### The Wang-Landau algorithm
start_time <- Sys.time()
for(iter in 1:num_mc){
  ### Metropolis move
  MC_result <- MC_sweep(lattice, logdensity, Hist, energy_level, logf)
  lattice <- MC_result$lattice
  logdensity <- MC_result$logdensity
  energy_level <- MC_result$energy_level
  Hist <- MC_result$Hist
  
  ### Check the flatness criterion (1/t update)
  if(iter%%1000 == 0 && Hist_index == 0){
    Hist_min <- min(Hist[-empty_index])
    if(Hist_min > 0){
      ### Adjust the modification factor and reset the histogram
      logf <- logf / 2
      Hist <- rep(0, L^2 + 1)
    }
  }
  if(iter >= 1000 & logf < 1/iter){Hist_index <- 1}
  if(Hist_index == 1){logf <- 1/iter}
  
  ### Normalize the density
  maxlw <- max(logdensity[-empty_index])
  logconst <- log(sum(exp(logdensity[-empty_index] - maxlw))) + maxlw
  logdensity[-empty_index] <- logdensity[-empty_index] - logconst
  
  ### Calculate the estimation error epsilon(t)
  error[iter] <- sum(abs(1 - logdensity[-empty_index]/logdensity_truth[-empty_index]))/(N - 1)
  if(iter%%1000 == 0) print(c(iter, logf, error[iter]))
}
end_time <- Sys.time()

### save data
data_save  <- list(logdensity = logdensity, error = error, computation_time = end_time - start_time)

