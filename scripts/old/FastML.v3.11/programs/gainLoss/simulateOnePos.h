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


#ifndef ___SIMULATE_1POS__
#define ___SIMULATE_1POS__

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "alphabet.h"
#include "threeStateModel.h"
#include "oneTwoMoreModel.h"


using namespace std;

/******************************************************************
Simulate one position using the 3stateLGT stochastic process
*******************************************************************/
class simulateOnePos{
public:
	//simulateOnePos();
	simulateOnePos(string simSeqFile, ostream* resFile, ostream* simulatedEvents, int simNum, string treeFile
		, MDOUBLE sumGainLoss, MDOUBLE theta
		, bool is3states=false, stochasticProcess* sp=NULL, tree* pTree=NULL
		, Vdouble* init_cpN_vals=NULL, Vdouble* freq_cpN=NULL);
	virtual ~simulateOnePos();
	
	VVint getChangesForBranch(int nodeID);
	sequenceContainer getSequenceContainer(){return _sc;};
	MDOUBLE getOccurFraction(){return _occurFraction;};


private:
	void init(string strTree);
	void init(tree* pTree);
	void simulateOnePosLGT(stochasticProcess* sp, string strOutFile);

	void simulateOnePos_cpN_Model(string strOutFile);
	void printTreeWithNodeIdBPStyle(ostream &out) const;
	void recursivePrintTree(ostream &out,const tree::nodeP &myNode) const;

private:
	tree _tree;
	stochasticProcess *_sp;
	sequenceContainer _sc; // as simulated
	int _simNum;
	MDOUBLE _sumGainLoss;
	MDOUBLE _theta;
	MDOUBLE _occurFraction;
	alphabet* _pAlph;
	vector<int> _alphVecDist;
	ostream *_simulatedEvents;
	ostream *_resFile;

	bool _simulateNullModel;
	bool _is3states;
	Vdouble* _init_cpN_vals;
	Vdouble* _freq_cpN;

	string _rootAt;
	VVVint _changesOccurred; // number of times changes from i to j occurred , for each branch
};


#endif

