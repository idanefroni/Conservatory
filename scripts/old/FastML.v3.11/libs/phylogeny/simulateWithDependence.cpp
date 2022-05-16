#include "simulateWithDependence.h"

/*
This code receives a tree file and simulates sequences accordingly using:
simulateTree st1(treeIn, *_sp, alph);
st1.generate_seq(num_pos_with_same_k);
which were written by another beloved group member.

Its feature is to simulate co-evolution between pairs of positions of binary data. Basic logic:
1. the basic concept is to use the regular independent model with 4 states to code a dependent model with 2 states.
   thus, all possible pairs of dada: 00, 01, 10, 11 are coded into A, C, G, T
2. dependency between possitions can be described as a tendency to have the same character (that is: 00 or 11). with this model we can accelerate
   the rate of evolution when an "unstable" state occures (rate increases when 01 (C) or 10 (G))

For more details, please see http://copap.tau.ac.il/benchmark.php and
Ofir Cohen, Haim Ashkenazy, Eli Levy Karin, David Burstein and Tal Pupko (2013)
CoPAP: Co-evolution of Presence-Absence Patterns.
Nucleic Acids Research 2013; doi: 10.1093/nar/gkt471

Eli Levy Karin, 2013

----------------- usage example: ------------------
#include "simulateWithDependence.h"

using namespace sim_with_dep;

int main(int argc, char** argv) {
	
	string treeFile = argv[1];

	double exit_code;

	exit_code = simulate_with_dependence (treeFile, 0.5, 14, 500, 500, 0, 1, 0.893195, 1, 4);
	
	return 0;
}
-------------- end usage example: ---------------

*/
namespace sim_with_dep 
{

	double simulate_with_dependence (string treeFile, double PI_1, double init_k, int total_positions, int num_pos_with_same_k, double k_increase, int is_gamma, double alpha, double beta, int num_cat)
	{

		//read Newick format tree
		tree treeIn(treeFile);
	
		//four states alphabet A C G T (will later be rplaced to 00,01,10,11)
		alphabet* alph = new nucleotide;

		sequenceContainer SC_all; //this will contain all positions
	
		//parameters:
		double PI_0 = 1-PI_1;
		double k = init_k; //will be increased with each iteration

		//parameters:
		int jump_size = total_positions / num_pos_with_same_k;

		for(int i=0; i<jump_size; i++)
		{
			Vdouble freqs; //stationary probabilities PI_00, PI_01, PI_10, PI_11
			double TOTAL = k*PI_1*PI_1 + 2*PI_0*PI_1 + k*PI_0*PI_0;
			freqs.push_back(k*PI_0*PI_0 / TOTAL); //PI_00 = k*PI_0*PI_0 / TOTAL
			freqs.push_back(PI_0*PI_1 / TOTAL); //PI_01 = PI_0*PI_1 / TOTAL
			freqs.push_back(PI_0*PI_1 / TOTAL); //PI_10 = PI_0*PI_1 / TOTAL
			freqs.push_back(k*PI_1*PI_1 / TOTAL); //PI_11 = k*PI_1*PI_1 / TOTAL
	
			//Q matrix (partial values - the rest are calculated by gtrModel using freqs and these values)
			MDOUBLE a2c = PI_1; // --> c2a = freqs[a]*a2c/freqs[c] --> c2a = ((k*PI_0*PI_0 / TOTAL)*PI_1)/(PI_0*PI_1 / TOTAL) = k*PI_0
			MDOUBLE a2g = PI_1;
			MDOUBLE a2t = 0;
			MDOUBLE c2g = 0;
			MDOUBLE c2t = k*PI_1;
			MDOUBLE g2t = k*PI_1;

			//starting the evolutionary model
			distribution *currDist = NULL;
			if(is_gamma == 1)
			{
				currDist = new generalGammaDistribution(alpha,beta,num_cat); // ---> in the future we might want to turn these into param
			}
			else
			{
				currDist = new uniDistribution; // no among site rate variation
			}
		
			replacementModel *probMod = NULL;
			pijAccelerator *pijAcc = NULL;

			probMod = new gtrModel(freqs,a2c,a2g,a2t,c2g,c2t,g2t);
			pijAcc = new trivialAccelerator(probMod);
			stochasticProcess* _sp = new stochasticProcess(currDist, pijAcc);

			//simulate:
			simulateTree st1(treeIn, *_sp, alph);
			st1.generate_seq(num_pos_with_same_k); //simulate num_pos_with_same_k positions with the current k

			if(i == 0)
			{
				SC_all = st1.toSeqDataWithoutInternalNodes(); //first time
			}
			else
			{
				sequenceContainer SC = st1.toSeqDataWithoutInternalNodes(); //concatenate new positions to the ones you have
				SC_all.concatenate(SC);
			}
		
			delete currDist;
			delete probMod;
			delete pijAcc;
			delete _sp;

			k = k + k_increase; //k = 1 , 1.05 , 1.1 , ... , 5.5
		}

		//prepare out file name:
		std::stringstream sstm;
		if(is_gamma == 1)
		{
			sstm << treeFile << ".gammaRateNoInv.PI_1=" << PI_1 << ".init_k=" << init_k << ".k_group_size=" << num_pos_with_same_k << ".k_increase=" << k_increase << ".fas";
		}
		else
		{
			sstm << treeFile << ".NoRate.PI_1=" << PI_1 << ".init_k=" << init_k << ".k_group_size=" << num_pos_with_same_k << ".k_increase=" << k_increase << ".fas";
		}
		std::string seqOutputFile = sstm.str();

		//write out:
		ofstream seq_sim(seqOutputFile.c_str());
		fastaFormat::write(seq_sim,SC_all);
		seq_sim.close();
	
	
		delete alph;
		
		return 0;
	
	}

};