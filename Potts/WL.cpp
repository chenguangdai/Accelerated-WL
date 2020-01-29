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
      // int itop = (i+1) % size;
      // int ibottom = ((i + size - 1) % size);
      // int jright = (j+1) % size;
      // int jleft = (j + size - 1) % size;
      // s += state(i,j) * (state(itop, j) + state(ibottom, j) + state(i, jright) + state(i, jleft));
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
    // Calculate sum over neighbor spins 
    itop = (si + 1) % size;
    ibottom = (si + size - 1) % size;
    jright = (sj + 1) % size;
    jleft = (sj + size - 1) % size;
    propose_spin = sample_int();
    current_spin = state(si, sj);
    energy_old = (current_spin == state(itop, sj)) + (current_spin == state(ibottom, sj)) + (current_spin == state(si, jright)) + (current_spin == state(si, jleft));
    energy_new = (propose_spin == state(itop, sj)) + (propose_spin == state(ibottom, sj)) + (propose_spin == state(si, jright)) + (propose_spin == state(si, jleft));
    // Rcpp::Rcout << energy_old << std::endl;
    // Rcpp::Rcout << energy_new << std::endl;
    // Rcpp::Rcout << propose_spin << std::endl;
    // Calculate the energy of the proposed state
    proposed_energy_level = current_energy_level + energy_old - energy_new;
    // Calculate the acceptance probability
    log_acceptance_prob = logdensity(current_energy_level) - logdensity(proposed_energy_level);
    if(log(runif(1)(0)) < log_acceptance_prob){
    	// Accept the probosed state
    	state(si, sj) = propose_spin;
    	current_energy_level = proposed_energy_level;
    }
    // Update logdensity, Histogram
    // logdensity(current_energy_level) += (logf * 1.0 / pow(size, 2));
    logdensity(current_energy_level) += logf;
    // logdensity(current_energy_level) += -1.0 / (size * size - 1);
    Hist(current_energy_level)++;
    // Normalize logdensity
    // logdensity = logdensity - log(sum(exp(logdensity)) - 2);
    // logdensity(1) = 0;
    // logdensity(size * size - 1) = 0;
  }
  // return
  return List::create(Named("lattice") = state, Named("energy_level") = current_energy_level, 
                      Named("logdensity") = logdensity, Named("Hist") = Hist);
}







