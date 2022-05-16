#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include "sequence.h"
#include "distribution.h"
#include "stochasticProcess.h"
#include "uniDistribution.h"
#include "nucJC.h"
#include "trivialAccelerator.h"
#include "sequenceContainer.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "likeDist.h"
// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {
	cout<<"This program computes for the JC model, when two sequences are given, and the distance between these two sequences is known, the likelihood."<<endl;
	string seqFile = "s2l4_DNA.txt";
	distribution *dist =  new uniDistribution; 
	replacementModel *probMod=new nucJC;
	pijAccelerator * pijAcc = new trivialAccelerator(probMod); 
	stochasticProcess sp(dist, pijAcc);
	ifstream in(seqFile.c_str());
	if (!in) {errorMsg::reportError("unable to open input sequence file");}
	nucleotide myAlph; 
	sequenceContainer original = phylipFormat::read(in,&myAlph);
	const MDOUBLE myToll = 0.0001;
	if (original.numberOfSeqs() != 2) {
		errorMsg::reportError("for this check, there suppose to be only 2 sequences",1);
	}
	sequence s1 = original[0];
	sequence s2 = original[1];
	likeDist likeDist1(sp,myToll);
	MDOUBLE resL = 0;
	MDOUBLE resD = likeDist1.giveDistance(s1,s2,NULL,&resL);
	cout<<" the likelihood of these 2 sequences is:"<<resL<<endl;
	cout<<" the ML distance between these 2 sequences is:"<<resD<<endl;

	delete dist;
	delete probMod;
	return 0;
}
