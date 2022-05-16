// simulate positions with dependence 2013 09 22 Eli Levy Karin

/*
This code receives a tree file and simulates sequences accordingly using:
simulateTree st1(treeIn, *_sp, alph);
st1.generate_seq(num_pos_with_same_k);
which were written by another beloved group member.

Its feature is to simulate co-evolution between pairs of positions of binary data. Basic logic:
1. the basic concept is to use the regular independent model with 4 states to code a dependent model with 2 states.
   thus, all possible pairs of dada: 00, 01, 10, 11 are coded into A, C, G, T
2. dependency between possitions can be described as a tendency to have the same character (that is: 00 or 11). with this model we can accelerate
   the rate of evolution when an "unstable" state occures (rate increases when 01 or 10)

For more details, please see http://copap.tau.ac.il/benchmark.php and
Ofir Cohen, Haim Ashkenazy, Eli Levy Karin, David Burstein and Tal Pupko (2013)
CoPAP: Co-evolution of Presence-Absence Patterns.
Nucleic Acids Research 2013; doi: 10.1093/nar/gkt471

Eli Levy Karin, 2013

*/

#ifndef ___SIM_WITH_DEP
#define ___SIM_WITH_DEP

#include <string>
#include <iostream>
#include "tree.h"
#include "alphabet.h"
#include "nucleotide.h"
#include "simulateTree.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h" // distribution of rates across sites
#include "generalGammaDistributionPlusInvariant.h"
#include "generalGammaDistribution.h"
#include "fastaFormat.h"
#include <fstream>
#include "gtrModel.h"
#include <sstream>

namespace sim_with_dep {
	double simulate_with_dependence (string treeFile, double PI_1, double init_k, int total_positions, int num_pos_with_same_k, double k_increase, int is_gamma, double alpha, double beta, int num_cat);
	/*
	treeFile - newick format
	total_positions - number of positions to simulate (note that you'll get double this nuber in binary positions)
	num_pos_with_same_k - one can have a different k for parts of the pairs (a gradient, for example). 
	Make sure that: (num_pos_with_same_k <= total_positions) and (total_positions % num_pos_with_same_k = 0)
	k_increase - if you decide to simulate with different k's - set by how much k should increase
	is_gamma - if 0 uniDistribution, if 1 generalGammaDistribution
	alpha, beta, number of rate categories are only relevant if is_gamma=1 (otherwise, you can put whatever there)
	*/
	
};

#endif
