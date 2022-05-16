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
#include "jcDistance.h"
#include "distanceTable.h"
#include "nj.h"
// NOTE: YOU MUST CHANGE THE NAME OF THE string seqFile TO MATCH YOUR OWN LOCATION OF THE SEQUENCE FILE NAME!

int main(int argc,char*argv[]) {
	cout<<"This program computes for the JC model, the NJ tree."<<endl;
	string seqFile = "nuc7.phylip.txt";
	if (argc>1)
	  seqFile=argv[1];
	distribution *dist =  new uniDistribution; 
	replacementModel *probMod=new nucJC;
	pijAccelerator * pijAcc = new trivialAccelerator(probMod); 
	stochasticProcess sp(dist, pijAcc);
	ifstream in(seqFile.c_str());
	if (!in) {errorMsg::reportError("unable to open input sequence file");}
	nucleotide myAlph; 
	sequenceContainer original = phylipFormat::read(in,&myAlph);
	
	
	//const MDOUBLE myToll = 0.0001;
	
	cout<<"computing the NJ tree..."<<endl;
	jcDistance likeDist1(myAlph.size());
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&likeDist1,
						original,
						disTab,
						vNames);
	NJalg nj1;
	tree njTree = nj1.computeTree(disTab,vNames);
	//	ofstream out("njTreeRes.txt");
	njTree.output(cout);



	//MDOUBLE resL = 0;
	//MDOUBLE resD = likeDist1.giveDistance(s1,s2,NULL,&resL);
	//cout<<" the likelihood of these 2 sequences is:"<<resL<<endl;
	//cout<<" the ML distance between these 2 sequences is:"<<resD<<endl;

	delete dist;
	delete probMod;
	return 0;
}
