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



#ifndef ___COMPUTE_POSTERIOR_EXPECTATION_OF_CHANGE
#define ___COMPUTE_POSTERIOR_EXPECTATION_OF_CHANGE

#include "definitions.h"
#include "simulateJumps.h"
#include "computeJumps.h"

#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "suffStatComponent.h"
#include "computePijComponent.h"

class computePosteriorExpectationOfChange {

public:
	explicit computePosteriorExpectationOfChange(const tree &tr, const sequenceContainer &sc, stochasticProcess *sp);
	virtual ~computePosteriorExpectationOfChange(){};


	VVdouble computeExpectationAcrossTree(simulateJumps &sim,  //input given from simulation studies
		const VVVdouble &posteriorProbs, VVVdouble &expForBranch);
	
	VVdouble computeExpectationAcrossTree(computeJumps &computeJumpsObj,  //Suchard
		const VVVdouble &posteriorProbs,VVVdouble &expForBranch);

	VVdouble computePosteriorAcrossTree(simulateJumps &sim, //input given from simulation studies
		const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch);

	VVdouble computePosteriorAcrossTree(computeJumps &computeJumpsObj, //Suchard
		const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch);

	
	void computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos);

private:
	MDOUBLE computePosteriorOfChangePerBranch(
		simulateJumps &sim, //input given from simulation studies
		const  VVVdouble &posteriorProbs,
		tree::nodeP node,
		int fromState, int toState);

	MDOUBLE computePosteriorOfChangePerBranch(
		computeJumps &computeJumpsObj, //Suchard
		const  VVVdouble &posteriorProbs,
		tree::nodeP node,
		int fromState, int toState);

	
	MDOUBLE computeExpectationOfChangePerBranch(
		simulateJumps &sim, //input given from simulation studies
		const VVVdouble &posteriorProbsGivenTerminals,
		tree::nodeP node,
		int fromState, int toState);
	
	MDOUBLE computeExpectationOfChangePerBranch(	//Suchard
		computeJumps &computeJumpsObj, //Suchard
		const VVVdouble &posteriorProbsGivenTerminals,
		tree::nodeP node,int fromState, int toState);

	
	MDOUBLE computePosterioGivenTerminalsPerBranch	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
		suffStatGlobalGamPos &sscDown,computePijHom &pi, doubleRep &LData, const string nodeName);


private:
	const tree &_tr;
	const sequenceContainer &_sc;
	stochasticProcess *_sp;
};


#endif
