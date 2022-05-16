#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;

#include "nucJC.h"
#include "sequence.h"
#include "distribution.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "uniDistribution.h"
#include "trivialAccelerator.h"
#include "sequenceContainer.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "likelihoodComputation.h"
#include "bestHKYparam.h"
#include "evaluateCharacterFreq.h"
#include "bblEM.h"
#include "allTrees.h"


int main(int argc,char*argv []) {
	cout<<"exhaustive search"<<endl;
	
	// getting the data.
	string seqFile1 = "nuc7.phylip.txt";
	ifstream in1(seqFile1.c_str());
	if (!in1) {errorMsg::reportError("unable to open input sequence file");}
	nucleotide myAlph; 
	sequenceContainer original1 = phylipFormat::read(in1,&myAlph);
	in1.close();

	distribution *dist =  new uniDistribution; 
	replacementModel *probMod=new nucJC;
	pijAccelerator * pijAcc = new trivialAccelerator(probMod); 
	stochasticProcess sp(dist, pijAcc);

	allTrees allTrees1(false);
	allTrees1.recursiveFind(&original1,&sp);
	cout<<" Log likelihood for best tree = "<<allTrees1.getBestScore()<<endl;
	allTrees1.getBestTree().output(cout);

	delete dist;
	delete probMod;
	return 0;
}
