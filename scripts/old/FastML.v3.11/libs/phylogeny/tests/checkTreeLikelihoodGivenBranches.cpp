#include <fstream>
#include <iostream>
#include <string>
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
	//	const MDOUBLE myToll = 0.0001;
	
	const string treeFileName = "sevenTaxaTree.txt";
	tree myT(treeFileName);
	cout<<"computing the log likelihood of the tree..."<<endl;

	MDOUBLE resL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(myT,original,sp);
	//tree njTree = nj1.computeNJtree(disTab,vNames);
	//ofstream out("njTreeRes.txt");
	//njTree.output(out);



	//MDOUBLE resL = 0;
	//MDOUBLE resD = likeDist1.giveDistance(s1,s2,NULL,&resL);
	cout<<" the likelihood of the tree is:"<<resL<<endl;
	//cout<<" the ML distance between these 2 sequences is:"<<resD<<endl;

	delete dist;
	delete probMod;
	return 0;
}
