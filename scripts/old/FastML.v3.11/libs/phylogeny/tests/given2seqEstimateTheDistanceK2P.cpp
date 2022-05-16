#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include "hky.h"
#include "sequence.h"
#include "distribution.h"
#include "stochasticProcess.h"
#include "uniDistribution.h"
#include "trivialAccelerator.h"
#include "sequenceContainer.h"
#include "nucleotide.h"
#include "phylipFormat.h"
#include "likeDist.h"
// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {
	cout<<"This program computes for the K2P model, when two sequences are given, the ML distance and its likelihood."<<endl;
	string seqFile = "s2l4_DNA.txt";
	distribution *dist =  new uniDistribution; 
	replacementModel *probMod1=new hky(0.25,0.25,0.25,0.25,0.5);
	pijAccelerator * pijAcc1 = new trivialAccelerator(probMod1); 
	replacementModel *probMod2=new hky(0.25,0.25,0.25,0.25,10);
	pijAccelerator * pijAcc2 = new trivialAccelerator(probMod2); 
	stochasticProcess sp1(dist, pijAcc1);
	stochasticProcess sp2(dist, pijAcc2);
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

	MDOUBLE resL1 = 0;
	likeDist likeDist1(sp1,myToll);
	MDOUBLE resD1 = likeDist1.giveDistance(s1,s2,NULL,&resL1);
	cout<<endl<<"For Tr/Tv = 0.5" <<endl;
	cout<<" the likelihood of these 2 sequences is:"<<resL1<<endl;
	cout<<" the ML distance between these 2 sequences is:"<<resD1<<endl;

	MDOUBLE resL2 = 0;
	likeDist likeDist2(sp2,myToll);
	MDOUBLE resD2 = likeDist2.giveDistance(s1,s2,NULL,&resL2);
	cout<<endl<<"For Tr/Tv = 10" <<endl;
	cout<<" the likelihood of these 2 sequences is:"<<resL2<<endl;
	cout<<" the ML distance between these 2 sequences is:"<<resD2<<endl;

	delete dist;
	delete probMod1;
	delete probMod2;
	return 0;
}
