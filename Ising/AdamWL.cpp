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
      s += state(i,j) * (state(i,j1) + state(i1,j));
    }
  }
  return -s;
}

// [[Rcpp::export]]
List MC_sweep(IntegerMatrix & state, NumericVector & logdensity, NumericVector & Hist, int energy_level, double learning_rate, double beta, NumericVector & last_update_index, NumericVector & momentum, int sweep_index){
  // last_update_index: the last time step, in terms of MC trial moves, that each evergy level is updated.
  // sweep_index: the number of MC sweeps that the algorithm has gone through. Each MC sweep contains LxL MC trial moves.
  // num_trial_moves: the number of MC trial moves that the algorithm has gone through.
  int size = state.rows();
  int current_energy_level = energy_level;
  int si, sj, site;
  int itop, ibottom, jright, jleft;
  int neighborsum, proposed_energy_level;
  double log_acceptance_prob;
  double accum_beta = sqrt(beta) / (1 - sqrt(beta));
  int interval;
  int num_trial_moves = (sweep_index - 1) * size * size;
	
  for (int iter = 0; iter < size * size; iter++){
    // Choose a random site for flipping
    site = (int)(size * size * (runif(1))(0));
    si = site / size;
    sj = site % size;
	  
    // Calculate the sum over neighboring spins 
    itop = (si + 1) % size;
    ibottom = (si + size - 1) % size;
    jright = (sj + 1) % size;
    jleft = (sj + size - 1) % size;
    neighborsum = state(itop, sj) + state(ibottom, sj) + state(si, jright) + state(si, jleft);
	  
    // Calculate the energy level of the proposed state
    proposed_energy_level = current_energy_level + state(si, sj) * neighborsum / 2;
	  
    // Calculate the accumulated momentum
    interval = num_trial_moves + iter - last_update_index(proposed_energy_level);
    logdensity(proposed_energy_level) += learning_rate * sqrt(-momentum(proposed_energy_level)) * (1 - pow(sqrt(beta), interval)) * accum_beta;
    momentum(proposed_energy_level) *= pow(beta, interval);
    last_update_index(proposed_energy_level) = num_trial_moves + iter;
    
    // Calculate the acceptance probability
    log_acceptance_prob = logdensity(current_energy_level) - logdensity(proposed_energy_level);
    if (log(runif(1)(0)) < log_acceptance_prob){
    	// Accept the probosed state
    	state(si, sj) *= -1;
    	current_energy_level = proposed_energy_level;
    }
    // Update momentum and logdensity
    momentum(current_energy_level) = momentum(current_energy_level) * beta - (1.0 - beta);
    logdensity(current_energy_level) += sqrt(-momentum(current_energy_level)) * learning_rate;
    last_update_index(current_energy_level) = num_trial_moves + iter + 1;
 
    // Update histogram
    Hist(current_energy_level)++;
  }
  // return
  return List::create(Named("energy_level") = current_energy_level, Named("logdensity") = logdensity, 
                      Named("lattice") = state, Named("Hist") = Hist, Named("last_update_index") = last_update_index,
                      Named("momentum") = momentum);
}








