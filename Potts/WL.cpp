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
int sample_int() {
    Rcpp::IntegerVector pool = Rcpp::seq(1, 10);
    std::random_shuffle(pool.begin(), pool.end());
    return pool[0];
} 

// [[Rcpp::export]]
List MC_sweep(IntegerMatrix & state, NumericVector & logdensity, NumericVector & Hist, int current_energy_level, double logf){
  int size = state.rows();
  int si, sj, site;
  int itop, ibottom, jright, jleft;
  int current_spin, propose_spin, energy_old, energy_new, proposed_energy_level;
  double log_acceptance_prob;
	
  for (int i = 0; i < size * size; i++){
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

    // Calculate the energy of the proposed state
    proposed_energy_level = current_energy_level + energy_old - energy_new;
    // Calculate the acceptance probability
    log_acceptance_prob = logdensity(current_energy_level) - logdensity(proposed_energy_level);
    if(log(runif(1)(0)) < log_acceptance_prob){
    	// Accept the probosed state
    	state(si, sj) = propose_spin;
    	current_energy_level = proposed_energy_level;
    }
	  
    // Update the logdensity and the Histogram
    logdensity(current_energy_level) += logf;
    Hist(current_energy_level)++;
  }
	
  // return
  return List::create(Named("lattice") = state, Named("energy_level") = current_energy_level, 
                      Named("logdensity") = logdensity, Named("Hist") = Hist);
}







