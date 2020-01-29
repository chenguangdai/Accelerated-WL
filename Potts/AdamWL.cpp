#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
int get_energy(const IntegerMatrix & state){
  int size = state.rows();
  int s = 0;
  int i1;
  int j1;
  for (int i = 0; i < size; i++){
    if( i==(size-1) ) i1=0;
    else i1 = i+1;
    for (int j = 0; j < size; j++){
      if( j==(size-1) ) j1=0;
      else j1 = j+1;
      if(state(i, j) == state(i, j1)) s += 1;
      if(state(i, j) == state(i1, j)) s += 1;
    }
  }
  return s;
}

// [[Rcpp::export]]
int sample_int() {
    Rcpp::IntegerVector pool = Rcpp::seq(1, 10);
    std::random_shuffle(pool.begin(), pool.end());
    return pool[0];
} 

// [[Rcpp::export]]
List MC_sweep(IntegerMatrix & state, NumericVector & logdensity, NumericVector & Hist, int energy_level, double learning_rate, double beta, NumericVector & last_update_index, NumericVector & momentum, int sweep_index){
  // last_update_index: the last time step, in terms of MC sweeps, that each evergy level is updated.
  // sweep_index: the number of MC sweeps that the algorithm has gone through. Each MC sweep contains LxL MC trial moves.
  int size = state.rows();
  int current_energy_level = energy_level;
  int si, sj, site;
  int itop, ibottom, jright, jleft;
  int current_spin, propose_spin, energy_old, energy_new, proposed_energy_level;
  double log_acceptance_prob;
  double accum_beta = sqrt(beta) / (1.0 - sqrt(beta));
  double interval;

  for (double iter = 0; iter < size * size; iter++){
    // Choose a random site for flipping
    site = (int)(size * size * (runif(1))(0));
    si = site / size;
    sj = site % size;
    
    // Calculate the sum over neighboring spins 
    itop = (si + 1) % size;
    ibottom = (si + size - 1) % size;
    jright = (sj + 1) % size;
    jleft = (sj + size - 1) % size;
    propose_spin = sample_int();
    current_spin = state(si, sj);
    energy_old = (current_spin == state(itop, sj)) + (current_spin == state(ibottom, sj)) + (current_spin == state(si, jright)) + (current_spin == state(si, jleft));
    energy_new = (propose_spin == state(itop, sj)) + (propose_spin == state(ibottom, sj)) + (propose_spin == state(si, jright)) + (propose_spin == state(si, jleft));
    proposed_energy_level = current_energy_level + energy_old - energy_new;
    
    // Calculate the accumulated momentum
    interval = (sweep_index - 1 - last_update_index(proposed_energy_level)) * size * size + iter;
    logdensity(proposed_energy_level) += learning_rate * sqrt(-momentum(proposed_energy_level)) * (1 - pow(sqrt(beta), interval)) * accum_beta;
    momentum(proposed_energy_level) *= pow(beta, interval);
    last_update_index(proposed_energy_level) = sweep_index - 1 + iter / (size * size);

    // Calculate the acceptance probability
    log_acceptance_prob = logdensity(current_energy_level) - logdensity(proposed_energy_level);
    if (log(runif(1)(0)) < log_acceptance_prob){
    	// Accept the proposed state
    	state(si, sj) = propose_spin;
    	current_energy_level = proposed_energy_level;
    }
    
    // Update the momentum and the logdensity
    momentum(current_energy_level) = momentum(current_energy_level) * beta - (1.0 - beta);
    logdensity(current_energy_level) += sqrt(-momentum(current_energy_level)) * learning_rate;
    last_update_index(current_energy_level) = sweep_index - 1 + (iter + 1) / (size * size);
    
    // Update the histogram
    Hist(current_energy_level)++;
  }
  // return
  return List::create(Named("energy_level") = current_energy_level, Named("logdensity") = logdensity, 
                      Named("lattice") = state, Named("Hist") = Hist, Named("momentum") = momentum, 
                      Named("last_update_index") = last_update_index);
  
}








