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
// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {

	
	cout<<"This program computes for the HKY model, the ML estimate of a given tree (when the branch lengths are given)."<<endl;
	string seqFile = "nuc7.phylip.txt";
//	distribution *dist =  new gammaDistribution(1,8); 
	ifstream in(seqFile.c_str());
	if (!in) {errorMsg::reportError("unable to open input sequence file");}
	nucleotide myAlph; 
	sequenceContainer original = phylipFormat::read(in,&myAlph);
	
	// check

	vector<MDOUBLE> myFreq = evaluateCharacterFreq(original);
	for (int j=0; j < myFreq.size(); ++j) {
		cout<<" the freq of nuc "<<j<<" is: "<<myFreq[j]<<endl;
	}
	
	distribution *dist =  new uniDistribution; 
	replacementModel *probMod=new hky(myFreq[0],myFreq[1],myFreq[2],myFreq[3],2);
	pijAccelerator * pijAcc = new trivialAccelerator(probMod); 
	stochasticProcess sp(dist, pijAcc);

	

	// end check
	const string treeFileName = "startTree.txt";
	tree myT(treeFileName);
	cout<<"computing the log likelihood of the tree..."<<endl;

	MDOUBLE  resL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(myT,original,sp);
	cout<<" starting L = "<<resL<<endl;
	bestHkyParamAndBBL bestHkyParamAndBBL1(myT,original,sp);
	
	resL = bestHkyParamAndBBL1.getBestL();
	MDOUBLE resAlpha = bestHkyParamAndBBL1.getBestHkyParam();
	cout<<"final likelihood: "<<resL<<endl;
	cout<<"best HKY parameter: "<<resAlpha<<endl;
	//	ofstream out("tAfterBBL.txt");
	myT.output(cout);

	delete dist;
	delete probMod;
	return 0;
}

//C:\tal\semphyCheck\nuc7.phylip.txt
//C:\tal\semphyCheck\startTree.txt
