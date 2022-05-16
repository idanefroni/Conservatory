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
#include "ancestralReconstructStates.h"
#include <cmath>

using namespace std;

/********************************************************************************************
ancestralReconstructStates
*********************************************************************************************/
ancestralReconstructStates::ancestralReconstructStates(const tree &tr, const sequenceContainer &sc, stochasticProcess *sp):
_tr(tr), _sc(sc){
	if(!sp){
		errorMsg::reportError("error in the constructor ancestralReconstructStates sp argument is NULL");
	}
	else{
		_sp = sp;
	}
	_statesV.resize(_sc.seqLen());
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		initializeStatesVector(pos);
	}
}
void ancestralReconstructStates::initializeStatesVector(int pos){
	_statesV[pos].resize(_tr.getNodesNum(),-1000);
	checkThatNamesInTreeAreSameAsNamesInSequenceContainer(_tr,_sc);
	seqContainerTreeMap scTreeMap(_sc,_tr);
	vector <tree::nodeP> leaves;
	_tr.getAllLeaves(leaves,_tr.getRoot());
	for (unsigned int i=0; i< leaves.size();i++){
		int myleafId = (leaves[i])->id();
		int mySeqId = scTreeMap.seqIdOfNodeI(myleafId);
		_statesV[pos][myleafId] = _sc[mySeqId][pos];
	}
}
/********************************************************************************************
upL[node][letter] = max(letter_here){P(letter->letter_here)*upL[son1][letter_here]*upL[son2][letter_here]} for letter at father node.
backtrack[node][letter] = argmax of above
*********************************************************************************************/
void ancestralReconstructStates::traverseUpML(VVVdouble &upL, VVVint &backtrack){ // input as empty vector to be filled
	LOGnOUT(4,<<"traverseUpML..."<<endl);
	upL.resize(_sc.seqLen());
	backtrack.resize(_sc.seqLen());
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		traverseUpML(upL[pos], backtrack[pos], pos);
	}
}
/********************************************************************************************
*********************************************************************************************/
void ancestralReconstructStates::traverseUpML(VVdouble &upL, VVint &backtrack, int pos){ // input as empty vector to be filled
	computePijGam pi;
	pi.fillPij(_tr,*_sp);
	upL.resize(_tr.getNodesNum());
	for (unsigned int i = 0; i < upL.size(); i++)
		upL[i].resize(_sp->alphabetSize());
	backtrack.resize(_tr.getNodesNum());
	for (unsigned int i = 0; i < backtrack.size(); i++)
		backtrack[i].resize(_sp->alphabetSize());
	treeIterDownTopConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int father_state = 0;
		if (mynode->isLeaf()) {
			for (father_state=0; father_state<_sp->alphabetSize();father_state++){ // looping over states at father
				int myState = _statesV[pos][mynode->id()];
				if(myState == _sc.getAlphabet()->unknown()){
					myState = father_state; // same as relations=1, for missing data
				}
				for (int i=0; i < _sp->categories();++i) {
					upL[mynode->id()][father_state] += pi.getPij(i,mynode->id(),father_state,myState)*_sp->ratesProb(i);
				}
				backtrack[mynode->id()][father_state]=myState;
			}
		}
		else if (!(mynode->isRoot())) {
			for (father_state=0; father_state<_sp->alphabetSize();father_state++){ // looping over states at father
				MDOUBLE myMax = -1;
				int myArgMax=-1;
				for (int my_state=0;my_state<_sp->alphabetSize();my_state++){ // loop to find max over current node
					//MDOUBLE val=_sp->Pij_t(father_state,my_state,mynode->dis2father());
					MDOUBLE val=0;
					for (int i=0; i < _sp->categories();++i) {
						val += pi.getPij(i,mynode->id(),father_state,my_state)*_sp->ratesProb(i);
					}
					for (int son=0;son<mynode->getNumberOfSons();son++)
						val*=upL[mynode->getSon(son)->id()][my_state];
					if (val>myMax){
						myMax=val;
						myArgMax=my_state;
					}
				}
				if ((myMax<0) || (myArgMax<0))
					errorMsg::reportError("Error in traverseUpML: cannot find maximum");
				upL[mynode->id()][father_state]=myMax;
				backtrack[mynode->id()][father_state]=myArgMax;
			}
		}
		else {// root
			for (int root_state=0; root_state<_sp->alphabetSize();root_state++){
				MDOUBLE val=_sp->freq(root_state);
				for (int son=0;son<mynode->getNumberOfSons();son++)
					val*=upL[mynode->getSon(son)->id()][root_state];
				upL[mynode->id()][root_state]=val;
			}
		}
	}
}

/********************************************************************************************
return likelihood of max joint reconstruction
*********************************************************************************************/
Vdouble ancestralReconstructStates::traverseDownML(VVVdouble &upL, VVVint &backtrack,VVVint &transitionTypeCount) { // input as already filled vector
	LOGnOUT(4,<<"traverseDownML..."<<endl);
	Vdouble LofJointV;
	LofJointV.resize(_sc.seqLen());
	transitionTypeCount.resize(_sc.seqLen());

	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		LofJointV[pos] = traverseDownML(upL[pos], backtrack[pos],transitionTypeCount[pos], pos);
	}
	return LofJointV;
}

/********************************************************************************************
fill _statesV, transitionTypeCount
*********************************************************************************************/
MDOUBLE ancestralReconstructStates::traverseDownML(VVdouble &upL, VVint &backtrack,VVint &transitionTypeCount, int pos) { // input as already filled vector
	if (backtrack.size() == 0)
		errorMsg::reportError("error in ancestralReconstruct::traverseDownML, input vector backtrack must be filled (call traverseUpML() first)");
	MDOUBLE LofJoint;
	int stateOfRoot;
	findMaxInVector(upL[(_tr.getRoot())->id()], LofJoint, stateOfRoot);
	_statesV[pos][(_tr.getRoot())->id()] = stateOfRoot;
	transitionTypeCount.resize(_sp->alphabetSize());
	for (unsigned int i = 0; i < transitionTypeCount.size(); i++)
		transitionTypeCount[i].resize(_sp->alphabetSize(),0);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot()) continue;
		int myId = mynode->id();
		int stateAtFather = _statesV[pos][mynode->father()->id()];
		int myState = _statesV[pos][mynode->id()];
		if(myState == _sc.getAlphabet()->unknown()){
			myState = stateAtFather; // same as relations=1, for missing data
		}
		if (mynode->isLeaf()) {
			transitionTypeCount[stateAtFather][myState]++;
			if ((_statesV[pos][mynode->id()]!=stateAtFather))
				LOG(7,<<"switch from "<<mynode->father()->name()<<"("<<stateAtFather<<") to "<<mynode->name()<<"("<<_statesV[pos][mynode->id()]<<")"<<endl);
			continue;
		}
		if(_statesV[pos][mynode->id()] == -2)
			cout<<_statesV[pos][mynode->id()]<<" unKnown at pos="<<pos<<" node="<<mynode->id()<<endl;
		_statesV[pos][mynode->id()]=backtrack[myId][stateAtFather];
		transitionTypeCount[stateAtFather][_statesV[pos][mynode->id()]]++;
	}
	return log(LofJoint);
}



/********************************************************************************************
compute Prob(letter at Node N is x|Data): the posterior probabilities at ancestral states 
Use the pre-calculated joint posterior probability P(N=x, father(N)=y|D) and just sum over these probs:
Prob(N=x|Data) = sum{fatherState}[P(N=x, father(N)=y|D)]}
stores results in member VVVdouble[pos][node][state] _ancestralProbs
use VVVVdouble _probChanges_PosNodeXY == jointPost[pos][nodeID][fatherLetter][letter]- after computePosteriorOfChangeGivenTerminals
*********************************************************************************************/
void ancestralReconstructStates::computeAncestralPosterior(const VVVVdouble& jointPost)
{
	LOGnOUT(4,<<"computeAncestralPosterior (take into acount joint probabilty)..."<<endl);
	int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	//int alphabetSizeForProbsSize = alphabetSize;
	//bool isThereMissingData = _sc.getAlphabetDistribution(true)[2]>0;
	//if(isThereMissingData)
	//	alphabetSizeForProbsSize++; // resize for one more

	_ancestralProbs.resize(_sc.seqLen());
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		resizeMatrix(_ancestralProbs[pos], numNodes, alphabetSize);
		treeIterTopDownConst tIt(_tr);
		int letter;
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
			if (mynode->isRoot()) {
				//for(letter = 0; letter<alphabetSize; ++letter) // itay's version - the computed vals are all 0
				//	_ancestralProbs[pos][mynode->id()][letter] = jointPost[pos][mynode->id()][0][letter];
				for(letter = 0; letter < alphabetSize; ++letter) {
					MDOUBLE sum = 0.0;
					for(int sonLetter = 0; sonLetter < alphabetSize; ++sonLetter) {
						sum += jointPost[pos][mynode->getSon(0)->id()][letter][sonLetter];	// sum over the son joint prob (instead of father)
					}
					_ancestralProbs[pos][mynode->id()][letter] = sum;
				}
				continue;
			}
			for(letter = 0; letter < alphabetSize; ++letter) {
				MDOUBLE sum = 0.0;
				for(int fatherLetter = 0; fatherLetter < alphabetSize; ++fatherLetter) {
					sum += jointPost[pos][mynode->id()][fatherLetter][letter];
				}
				_ancestralProbs[pos][mynode->id()][letter] = sum;
			}
		}
	}	
}

