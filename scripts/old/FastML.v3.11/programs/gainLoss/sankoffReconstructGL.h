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


#ifndef ___SANKOFF__GL__H
#define ___SANKOFF__GL__H


#include "tree.h"
#include "logFile.h"
#include "someUtil.h"
#include "definitions.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "gainLossUtils.h"
#include <map>


class sankoffReconstructGL {

public:
	explicit sankoffReconstructGL(sequenceContainer& sc, tree& tr, string& outDir, MDOUBLE costMatrixGainLossRatio,  MDOUBLE distanceFromRootForRecent);
	virtual ~sankoffReconstructGL() ;
	void traverseUpMP(VVdouble &upCosts, vector <VVint> &backtrack); // input as empty vector to be filled
	MDOUBLE traverseDownMP(VVdouble &upCosts, vector <VVint> &backtrack, VVint &transitionTypeCount, VVdouble &totalCosts); // input as already filled vector
	Vdouble getGainMPPerPos(){return _gainMPPerPos;}
	Vdouble getLossMPPerPos(){return _lossMPPerPos;}
	VVVdouble getMPPerPos(){return _MPPerPos;}	
	VVVdouble getMPPerBranch(){return _MPPerBranch;}
	VVVVdouble getMPPerPosPerNode(){return _MPPerPosPerNode;}
	int getNumOfGainEvnetsMP(){return _numOfGains;}
	int getNumOfLossEvnetsMP(){return _numOfLosses;}



private:
	void initialize();
	void run();
	void startTree();
	void startSequenceContainer();
	void startCostMatrix();
	MDOUBLE runPosition(int pos, ofstream& gainLossMPPerPosPerBranchStream, ofstream& MPprints, ofstream& gainLossMPAncestralReconstructStream);
	void preparePrintData(Vstring &data);//prepares the data to be printed as BP data on the tree

	void printMPPerBranch(ostream& out);
	void printMPPerPos(ostream& out);


public:


private:
	VVdouble _costMatrix;
	Vint _states; // the vector with the states of the leaves, to be filled with reconstructed states
	alphabet * _alph; 
	tree _tr;
	sequenceContainer _sc;
	MDOUBLE _costOfTree;
	int _numOfGains;
	int _numOfLosses;

	Vdouble _lossMPPerPos;
	Vdouble _gainMPPerPos;
	VVVdouble _MPPerPos;
	VVVdouble _MPPerBranch;
	VVVVdouble _MPPerPosPerNode;


	MDOUBLE _distanceFromRootForRecent;
	MDOUBLE _costMatrixGainLossRatio;
	string _outDir;
};


#endif
