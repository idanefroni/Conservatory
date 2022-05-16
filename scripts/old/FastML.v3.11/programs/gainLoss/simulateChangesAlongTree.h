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


#ifndef ___SIMULATE_CHANGES__
#define ___SIMULATE_CHANGES__

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"
#include "sequenceContainer.h"

#include <map>
#include <vector>
using namespace std;

/******************************************************************
This class simulates jumps (events) along  a 
given tree, with the aim of creating a dataset (seqeunceContainer) 
in which we know the exact number of transitions along the tree
*******************************************************************/

class simulateChangesAlongTree  {
public:
	simulateChangesAlongTree(const tree& inTree, const stochasticProcess& sp, alphabet* pAlph);
	virtual ~simulateChangesAlongTree();
	sequenceContainer simulatePosition(); 
	VVint getChangesForBranch(int nodeID);
	int getNodeContent(int nodeId) {return _nodesContent[nodeId];}
	void removeAllSequnces(){
		_sc.removeAll();
	};

private:
	void init();
	void simulateOnce(tree::nodeP curNode, MDOUBLE disFromNode, int previousContent, int whichSon = 0);
	

private:
	tree _tree;
	stochasticProcess _sp;
	alphabet* _pAlph;

	Vdouble _waitingTimeParams;//each entry is the lambda parameter of the exponential distribution modeling the waiting time for "getting out" of state i
	//_jumpProbs[i][j] is the probability of jumping from state i to state j (given that a change has ocured).
	VVdouble _jumpProbs; 

	VVVint _changesOccurred; // number of times changes from i to j occurred , for each branch
	Vint _nodesContent; // the actual state at each node, retrieval according to node id
	sequenceContainer _sc;

};

#endif
