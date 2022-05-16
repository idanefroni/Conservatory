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
#include "simulateChangesAlongTree.h"
#include "talRandom.h"
#include "matrixUtils.h"
#include "gainLoss.h"

#include <algorithm>


simulateChangesAlongTree::simulateChangesAlongTree(const tree& inTree, const stochasticProcess& sp, alphabet* pAlph)
: _tree(inTree), _sp(sp), _pAlph(pAlph)	
{
}

simulateChangesAlongTree::~simulateChangesAlongTree()
{
}




void simulateChangesAlongTree::init()
{
	//init the vector of waiting times. 
	_waitingTimeParams.clear();
	_waitingTimeParams.resize(_pAlph->size());
	int i, j;
	for (i = 0; i < _pAlph->size(); ++i)
	{
		_waitingTimeParams[i] = -_sp.dPij_dt(i, i, 0.0);
		
	}

	//init _jumpProbs.
	//_jumpProbs[i][j] = Q[i][j] / -Q[i][i]
	_jumpProbs.clear();
	_jumpProbs.resize(_pAlph->size());
	for (i = 0; i < _pAlph->size(); ++i)
	{
		MDOUBLE sum = 0.0;
		_jumpProbs[i].resize(_pAlph->size());
		for (j = 0; j < _pAlph->size(); ++j)
		{
			if (i == j)
				_jumpProbs[i][j] = 0.0;
			else
			{
				_jumpProbs[i][j] = _sp.dPij_dt(i, j, 0.0) / _waitingTimeParams[i];
			}
			sum += _jumpProbs[i][j];
		}
		if (! DEQUAL(sum, 1.0)){
			string err = "error in simulateJumps::init(): sum probabilities is not 1 and equal to ";
			err+=double2string(sum);
			errorMsg::reportError(err);
		}
	}
	int nodesNum = _tree.getNodesNum();
	_changesOccurred.clear();
	_changesOccurred.resize(nodesNum);
	for (int i=0; i<nodesNum; ++i)
		resizeMatrix(_changesOccurred[i], _pAlph->size(), _pAlph->size());
	_nodesContent.clear();
	_nodesContent.resize(nodesNum, 0);
}

sequenceContainer simulateChangesAlongTree::simulatePosition(){
	init();
	Vdouble freqs(_pAlph->size(),0.0);
	for (int i = 0; i< freqs.size(); ++i)
		freqs[i]=_sp.freq(i);
	int rootState = giveRandomState(_pAlph->size(), freqs);
	//int rootState = giveRandomState(_pAlph, freqs);

	_nodesContent[_tree.getRoot()->id()] = rootState;
	simulateOnce(_tree.getRoot(),0,rootState,0);
	simulateOnce(_tree.getRoot(),0,rootState,1);
	if (_tree.getRoot()->getNumberOfSons() > 2)
		simulateOnce(_tree.getRoot(),0,rootState,2);
	return _sc;
}

void simulateChangesAlongTree::simulateOnce(tree::nodeP curNode, 
											MDOUBLE disFromNode, 
											int previousContent, int whichSon){
		tree::nodeP sonNode = curNode->getSon(whichSon);
		MDOUBLE avgWaitingTime = 1.0 / _waitingTimeParams[previousContent];
		MDOUBLE timeTillChange = talRandom::rand_exp(avgWaitingTime);
		disFromNode += timeTillChange;
		//int nextContent = giveRandomState(_pAlph, previousContent, _jumpProbs);
		int nextContent = giveRandomState(_pAlph->size(), previousContent, _jumpProbs);

		while (disFromNode < sonNode->dis2father()) {
			_changesOccurred[sonNode->id()][previousContent][nextContent]++;
			previousContent=nextContent;
			MDOUBLE avgWaitingTime = 1.0 / _waitingTimeParams[previousContent];
			MDOUBLE timeTillChange = talRandom::rand_exp(avgWaitingTime);
			disFromNode += timeTillChange;
			//nextContent = giveRandomState(_pAlph, nextContent, _jumpProbs);
			nextContent = giveRandomState(_pAlph->size(), nextContent, _jumpProbs);

		}
		while (disFromNode >= sonNode->dis2father()) {
			_nodesContent[sonNode->id()] = previousContent;
			if (sonNode->isLeaf()) {
				//string name = "leaf_" + int2string(sonNode->id()) + "_" + sonNode->name();
				string name = sonNode->name();
				sequence seq(int2string(previousContent),name, "", sonNode->id(), _pAlph);
				_sc.add(seq);
				return;
			}
			simulateOnce(sonNode, 0, previousContent, 1);
			disFromNode-=sonNode->dis2father();
			curNode = sonNode;
			sonNode = curNode->getSon(0);
		}
		_changesOccurred[sonNode->id()][previousContent][nextContent]++;
		simulateOnce(curNode, disFromNode, nextContent, 0);
}

VVint simulateChangesAlongTree::getChangesForBranch(int nodeID){
	if (nodeID>_changesOccurred.size())
		errorMsg::reportError("error in simulateChangesAlongTree::getChangesForBranch, nodeID doesn't exist");
	return _changesOccurred[nodeID];
}