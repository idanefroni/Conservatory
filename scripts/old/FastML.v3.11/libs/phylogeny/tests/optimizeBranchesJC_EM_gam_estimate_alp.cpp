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
#include "trivialAccelerator.h"
#include "sequenceContainer.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "likelihoodComputation.h"
#include "bestAlpha.h"

// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {
	cout<<"This program computes for the JC model, the likelihood of a given tree (when the branch lengths are given)."<<endl;
	string seqFile = "nuc7.phylip.txt";
	distribution *dist =  new gammaDistribution(1,8); 
	replacementModel *probMod=new nucJC;
	pijAccelerator * pijAcc = new trivialAccelerator(probMod); 
	stochasticProcess sp(dist, pijAcc);
	ifstream in(seqFile.c_str());
	if (!in) {errorMsg::reportError("unable to open input sequence file");}
	nucleotide myAlph; 
	sequenceContainer original = phylipFormat::read(in,&myAlph);
	
	const string treeFileName = "startTree.txt";
	tree myT(treeFileName);
	cout<<"computing the log likelihood of the tree..."<<endl;

	MDOUBLE  resL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(myT,original,sp);
	cout<<" starting L = "<<resL<<endl;
	bestAlphaAndBBL bestAlphaAndBBL1(myT,original,sp);
	
	resL = bestAlphaAndBBL1.getBestL();
	MDOUBLE resAlpha = bestAlphaAndBBL1.getBestAlpha();
	cout<<"final likelihood: "<<resL<<endl;
	cout<<"best Alpha: "<<resAlpha<<endl;
	//	ofstream out("tAfterBBL.txt");
	myT.output(cout);

	delete dist;
	delete probMod;
	return 0;
}
