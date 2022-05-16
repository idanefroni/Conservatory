/*
Copyright (C) 2011 Tal Pupko  TalP@tauex.tau.ac.il.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "definitions.h"
#include "simulateOnePos.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"
#include "simulateTree.h"
#include "threeStateAlphabet.h"
#include "recognizeFormat.h"
#include "evaluateCharacterFreq.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h"
#include "sequence.h"
#include "simulateChangesAlongTree.h"
#include "treeIt.h"
#include "fastaFormat.h"
#include "gainLoss.h"

#include <fstream>
#include <string>
using namespace std;

/********************************************************************************************
*********************************************************************************************/
simulateOnePos::simulateOnePos(string simSeqFile, ostream* resFile, ostream* simulatedEvents, int simNum, string treeFile
							   , MDOUBLE sumGainLoss, MDOUBLE theta							   
							   , bool is3states, stochasticProcess* sp, tree* pTree
							   , Vdouble* init_cpN_vals, Vdouble* freq_cpN)
: _pAlph(NULL),_simulateNullModel(false),_simNum(simNum),_theta(theta),_sumGainLoss(sumGainLoss)
	,_is3states(is3states), _init_cpN_vals(init_cpN_vals),_freq_cpN(freq_cpN),_simulatedEvents(simulatedEvents),_resFile(resFile)  {
	if(pTree!=NULL)
		init(pTree);		
	else
		init(treeFile);
		

	if (_simulateNullModel)
		simulateOnePos_cpN_Model(simSeqFile);
	else
		simulateOnePosLGT(sp,simSeqFile);
}
/********************************************************************************************
*********************************************************************************************/
simulateOnePos::~simulateOnePos()
{
	if (_sp)
		delete _sp;
	if (_pAlph)
		delete _pAlph;
	//if (_out)
	//	delete _out;
	//if (_res)		
	//	delete _res;
	//if (_outTree)		
	//	delete _outTree;
}

/********************************************************************************************
*********************************************************************************************/
void simulateOnePos::init(string strTree)
{
	_tree = tree(strTree);
	if (!(_rootAt =="")){
		tree::nodeP myroot = _tree.findNodeByName(_rootAt); //returns NULL if not found
		if (myroot){
			_tree.rootAt(myroot);
			//*_res<<"# tree rooted at "<<myroot->name()<<" id, "<<myroot->id()<<endl;
			*_simulatedEvents<<"# tree rooted at "<<myroot->name()<<" id, "<<myroot->id()<<endl;
		}
		else {
			//*_res<<"# default rooting used "<<endl;
			*_simulatedEvents<<"# default rooting used "<<endl;
		}
	}
	if(_is3states)
		_pAlph = new threeStateAlphabet();
	else
		_pAlph = new gainLossAlphabet();
	_alphVecDist.resize(_pAlph->size());
	_changesOccurred.resize(_tree.getNodesNum());
	for (int i=0; i<_tree.getNodesNum(); ++i)
		resizeMatrix(_changesOccurred[i], _pAlph->size(), _pAlph->size());

}
/********************************************************************************************
*********************************************************************************************/
void simulateOnePos::init(tree* pTree)
{
	_tree = *pTree;
	if(_is3states)
		_pAlph = new threeStateAlphabet();
	else
		_pAlph = new gainLossAlphabet();
	_alphVecDist.resize(_pAlph->size());
	_changesOccurred.resize(_tree.getNodesNum());
	for (int i=0; i<_tree.getNodesNum(); ++i)
		resizeMatrix(_changesOccurred[i], _pAlph->size(), _pAlph->size());

}

/********************************************************************************************
*********************************************************************************************/
void simulateOnePos::simulateOnePos_cpN_Model(string strOutFile) {
	Vdouble freq(2,0.0);/// FILL IN!!!
	freq[0]= 0.6;
	freq[1]= 0.4;
	
	MDOUBLE init_gain = 0.0; // No HGT
	MDOUBLE init_loss = 3.23;
	bool _isHGT_normal_Pij = true;
	bool _isHGT_with_Q = true;
	//gainLossModel  glm(init_gain,freq,_isHGT_normal_Pij,_isHGT_with_Q);
	gainLossModelNonReversible  glm(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,_isHGT_normal_Pij,_isHGT_with_Q);
	trivialAccelerator pijAcc(&glm);
	uniDistribution uniDistr;
	_sp = new stochasticProcess(&uniDistr,&pijAcc,false);

	// simulate:
	simulateTree st1(_tree, *_sp, _pAlph);
	Vdouble rates(1,1.0);
    st1.generate_seqWithRateVector(rates,1);

	_sc = st1.toSeqDataWithoutInternalNodes();	
	ofstream seq_sim(strOutFile.c_str());
	seq_sim.precision(PRECISION);
	fastaFormat::write(seq_sim,_sc);
	seq_sim.close();
}

/********************************************************************************************
*********************************************************************************************/
void simulateOnePos::simulateOnePosLGT(stochasticProcess* sp, string strOutFile) 
{
	if(!sp){
		if(_is3states){
			Vdouble init_cpN_vals(4);
			if(_init_cpN_vals){
				init_cpN_vals = *_init_cpN_vals;
			}
			else{
				init_cpN_vals[0]=0.25; //gain (0->1)
				init_cpN_vals[1]=1; //more (1->more)
				init_cpN_vals[2]=1; // less (more->1) 
				init_cpN_vals[3]=0.5; // loss (1->0)
			}
			if(_simNum==0)// printed once only
				LOGnOUT(3,<<"Rate values: gain (0->1)="<<init_cpN_vals[0]<<" more (1->more)="<<init_cpN_vals[1]
									<<" less (more->1)="<<init_cpN_vals[2]<<" loss (1->0)="<<init_cpN_vals[3]<<"\n");
			Vdouble freq_cpN(3);
			if(_freq_cpN)
				freq_cpN = *_freq_cpN;
			else{
				freq_cpN[0]=0.5;
				freq_cpN[1]=0.2;
				freq_cpN[2]=1 - (freq_cpN[0] + freq_cpN[1]);
			}
			bool useMarkovLimiting = false;
			if(_simNum==0){
				LOGnOUT(3,<<"Freq values: 0="<<freq_cpN[0]<<" 1="<<freq_cpN[1]<<" more="<<freq_cpN[2]<<" loss (1->0)="<<init_cpN_vals[3]<<"\n");
				LOGnOUT(3,<<"Freq useMarkovLimiting="<<useMarkovLimiting<<"\n");}

			oneTwoMoreModel glm_cpN(init_cpN_vals[0],init_cpN_vals[1],
				init_cpN_vals[2],init_cpN_vals[3],freq_cpN,useMarkovLimiting);
			trivialAccelerator pijAcc_cpN(&glm_cpN);
			uniDistribution uniDistr_cpN;
			bool isRevers = false;
			_sp = new stochasticProcess(&uniDistr_cpN,&pijAcc_cpN,false);
			MDOUBLE sumQii=(static_cast<oneTwoMoreModel*>(_sp->getPijAccelerator()->getReplacementModel()))->sumPijQij();
			(static_cast<oneTwoMoreModel*>(_sp->getPijAccelerator()->getReplacementModel()))->norm(1/sumQii);
			//cout<<" sumQii before norm="<<sumQii<<"\n";
		}
		else{
			// frequencies taken as estimated from the stationary distribution of the stochastic process
			LOGnOUT(3,<<"ERROR: simulateOnePosLGT with no stochastic process. Use constant default parameters\n");
			Vdouble freq(2,0.0);
			freq[0]= 0.6;
			freq[1]= 0.4;
			MDOUBLE init_gain = 0.942; // taken from original runs of COG data under the lgt model
			MDOUBLE init_loss = 5.23;
			bool _isHGT_normal_Pij = true;
			bool _isHGT_with_Q = true;
			//gainLossModel  glm(init_gain,freq,_isHGT_normal_Pij,_isHGT_with_Q);
			gainLossModelNonReversible  glm(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,_isHGT_normal_Pij,_isHGT_with_Q);
			trivialAccelerator pijAcc(&glm);
			uniDistribution uniDistr;
			_sp = new stochasticProcess(&uniDistr,&pijAcc,false);
			MDOUBLE sumQii = 1.0;
			sumQii = normalizeQ(_sp);
		}
	}
	else{
		_sp = sp->clone();
	}

	simulateChangesAlongTree sim(_tree,*_sp,_pAlph);
	_sc = sim.simulatePosition();
	_alphVecDist = _sc.getAlphabetDistribution();
	bool isFinishOneRun = false;
	do{		
		if(isFinishOneRun)
			LOGnOUT(6,<<"The number of 1s simulated "<< _alphVecDist[1]<<" was less than "<<gainLossOptions::_minNumOfOnes<<"\n");
		for(int alph = 0 ; alph<_pAlph->size() ;++alph){
			LOGnOUT(6,<<_alphVecDist[alph]<<" ");
		}
		LOGnOUT(6,<<"\n");
		sim.removeAllSequnces();
		_sc = sim.simulatePosition();
		_alphVecDist = _sc.getAlphabetDistribution();
		isFinishOneRun = true;
	}
	while(_alphVecDist[1]< gainLossOptions::_minNumOfOnes);
	_occurFraction = (float)_alphVecDist[1]/(float)_sc.numberOfSeqs();

	ofstream seq_sim(strOutFile.c_str());
	seq_sim.precision(PRECISION);
	fastaFormat::write(seq_sim,_sc);
	seq_sim.close();



	treeIterTopDownConst tit(_tree);
	int totalNumChangesInTree = 0;
	//*_res<<"# print values by simulations "<<endl;
	//*_res<<"G/L"<<"\t"<<"SIM"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"sumGainLoss"<<"\t"<<"rootFq"<<"\t"<<"occur"<<"\t"<<"events"<<endl;
	if(_simNum+1==1){ // print only at first position
		*_simulatedEvents<<"# print values by simulations "<<endl;
		*_simulatedEvents<<"G/L"<<"\t"<<"SIM"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"sumGainLoss"<<"\t"<<"rootFq"<<"\t"<<"occur"<<"\t"<<"events"<<endl;
		*_resFile<<"branch"<<"\t"<<"positions"<<"\t"<<"state"<<endl;
	}	
	for (tree::nodeP myN = tit.first();myN!=tit.end(); myN = tit.next()) {
		*_resFile<<myN->name()<<"\t"<<_simNum+1<<"\t"<< sim.getNodeContent(myN->id())<<endl;
		if(myN->isRoot())
			continue;
		VVint changesInNode = sim.getChangesForBranch(myN->id());
			_changesOccurred[myN->id()] = changesInNode;
			//*_res<<"Node id="<<myN->id()<<" name="<<myN->name()<< " content=" << sim.getNodeContent(myN->id()) << endl;
			for (int i=0; i<changesInNode.size(); ++i) {
				for (int j=0; j<changesInNode.size(); ++j) {
					totalNumChangesInTree+=changesInNode[i][j];
					//if(changesInNode[i][j]>0) // DEBUG
					//	cout<<"total number of changes: "<<totalNumChangesInTree<<" "<<myN->name()<<" "<<i<<" "<<j<<"\n";

					//*_res<<changesInNode[i][j]<<" ";
					if((i==0)&&(j==1)&&(changesInNode[i][j]>0)){
						//*_res<<"gain"<<"\t"<<_simNum+1<<"\t"<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()
						//	<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()
						//	<<"\t"<<_sumGainLoss<<"\t"<<_theta<<"\t"<<occur<<"\t"<<changesInNode[i][j]<<endl;
						*_simulatedEvents<<"gain"<<"\t"<<_simNum+1<<"\t"<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()
							<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()
							<<"\t"<<_sumGainLoss<<"\t"<<_theta<<"\t"<<_occurFraction<<"\t"<<changesInNode[i][j]<<endl;
					}
					if((i==1)&&(j==0)&&(changesInNode[i][j]>0)){ //NOTE: in both gain and loss use changesInNode[i][j] for event indication
						//*_res<<"loss"<<"\t"<<_simNum+1<<"\t"<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()
						//	<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()
						//	<<"\t"<<_sumGainLoss<<"\t"<<_theta<<"\t"<<_occurFraction<<"\t"<<changesInNode[i][j]<<endl;
						*_simulatedEvents<<"loss"<<"\t"<<_simNum+1<<"\t"<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()
							<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()
							<<"\t"<<_sumGainLoss<<"\t"<<_theta<<"\t"<<_occurFraction<<"\t"<<changesInNode[i][j]<<endl;
					}
				}
				//*_res<<endl;
			}

		//*_res<<"TOTAL Number of changes along the tree were "<<totalNumChangesInTree<<endl;
	}
	//printTreeWithNodeIdBPStyle(*_outTree);  // WARN - the removal of this print was not tested
}

// copied from the original covarion code , from the file checkov.cpp
void simulateOnePos::printTreeWithNodeIdBPStyle(ostream &out) const{
	recursivePrintTree(out,_tree.getRoot());
	out<<";";
}

// similar to the file checkov.cpp from the original covarion code
// The bootstrap values is the nodes id.
void simulateOnePos::recursivePrintTree(ostream &out,const tree::nodeP &myNode) const {
	if (myNode->isLeaf()) {
		out << myNode->name() << "_" << myNode->id();
		out << ":"<< myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			recursivePrintTree(out, myNode->getSon(i));
		}
		out <<")";
		if (myNode->isRoot()==false) {
			out<<":"<< myNode->dis2father();
			out << "["<<myNode->id()<<"]";
		}
	}
}


VVint simulateOnePos::getChangesForBranch(int nodeID){
	if (nodeID>_changesOccurred.size())
		errorMsg::reportError("error in simulateChangesAlongTree::getChangesForBranch, nodeID doesn't exist");
	return _changesOccurred[nodeID];
}
