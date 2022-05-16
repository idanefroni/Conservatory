// $Id: optimizeBranchesJC_EM.cpp 370 2005-05-25 20:46:34Z ninio $ 

#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;

#include "nucJC.h"
#include "sequence.h"
#include "distribution.h"
#include "stochasticProcess.h"
#include "uniDistribution.h"
#include "trivialAccelerator.h"
#include "sequenceContainer.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "likelihoodComputation.h"
#include "bblEM.h"

// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {
	cout<<"This program computes for the JC model, the likelihood of a given tree (when the branch lengths are given)."<<endl;
	string seqFile = "nuc7.phylip.txt";
	distribution *dist =  new uniDistribution; 
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

	bblEM bblEM1(myT,original,sp,NULL,1000,0.05);
	resL = bblEM1.getTreeLikelihood();

	cout<<" end L, after bbl = "<<setprecision(12)<<resL<<endl;
									
	//	ofstream out("tAfterBBL.txt");
	myT.output(cout);

	delete dist;
	delete probMod;
	return 0;
}
