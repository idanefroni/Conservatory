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
#include "computePosteriorExpectationOfChange.h"
#include "definitions.h"
#include "computeDownAlg.h"
#include "computeUpAlg.h"
#include "matrixUtils.h"
#include "treeIt.h"
#include "likelihoodComputation.h"
#include "gainLossOptions.h"
#include "gainLossModel.h"
#include "definitions.h"

using namespace std;

/********************************************************************************************
computePosteriorExpectationOfChange
*********************************************************************************************/
computePosteriorExpectationOfChange::computePosteriorExpectationOfChange(const tree &tr, const sequenceContainer &sc, stochasticProcess *sp):
_tr(tr), _sc(sc){
	if(!sp){
		errorMsg::reportError("error in the constructor computePosteriorExpectationOfChange sp argument is NULL");
	}
	else{
		_sp = sp;
	}
}
/********************************************************************************************
Expectation of number of changes from character u to v --- =
sum over all changes x,y:
Posterior(Node=x,Father=y|D)*Exp(changes u to v|Node=x,Father=y)
The second term is given to the function as input (can be obtained via simulations)
*********************************************************************************************/
VVdouble computePosteriorExpectationOfChange::computeExpectationAcrossTree(
	simulateJumps &sim,  //input given from simulation studies
	const VVVdouble &posteriorProbs,
	VVVdouble &expForBranch)
{
	int alphabetSize = _sp->alphabetSize();
	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				expForBranch[mynode->id()][fromState][toState] = computeExpectationOfChangePerBranch(sim,posteriorProbs,mynode,fromState,toState);
				res[fromState][toState] +=expForBranch[mynode->id()][fromState][toState];

			}
		}
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfChange::computeExpectationOfChangePerBranch(
	simulateJumps &sim, //input given from simulation studies
	const VVVdouble &posteriorProbsGivenTerminals,
	tree::nodeP node,int fromState, int toState)
{
	int alphabetSize = _sp->alphabetSize();

	MDOUBLE nodeExpectation = 0;
	for (int x = 0; x<alphabetSize; ++x){
		for (int y = 0; y<alphabetSize; ++y){
			nodeExpectation+=(posteriorProbsGivenTerminals[node->id()][x][y]*
				sim.getExpectation(node->name(),x,y,fromState,toState));
			
			if(node->name()=="A" && x==0){	//DEBUG
				LOG(9,<<"node "<<node->name()<<" from "<<fromState<<" to "<<toState<<" given "<<x<<" and "<<y
					<<" sim= "<< sim.getExpectation(node->name(),x,y,fromState,toState)
					<<" sim*'Pij'= "<< sim.getExpectation(node->name(),x,y,fromState,toState)
					*sim.getTotalTerminal(node->name(),x,y)/gainLossOptions::_numOfSimulationsForPotExp
					//*sim.getTotalTerminal(node->name(),x,y)/(sim.getTotalTerminal(node->name(),x,1)+sim.getTotalTerminal(node->name(),x,0))
					//<<" terminals x and y= "<< sim.getTotalTerminal(node->name(),x,y)<<" terminal start x= "<< (sim.getTotalTerminal(node->name(),x,1)+sim.getTotalTerminal(node->name(),x,0))
					<<endl);
			}
		}
	}
	if(node->name()=="A" ){ // DEBUG 
		//LOGnOUT(9,<<"Sim node A "<<node->dis2father()<<" from "<<fromState<<" to "<<toState<<" exp="<<expGivenStart0nodeA<<endl);	// DEBUG
		LOGnOUT(9,<<"nodeExpectation fromState "<<fromState<<" toState "<<toState<<" = "<<nodeExpectation<<endl);
	}
	return nodeExpectation;
}

/********************************************************************************************
Posterior probabilities computed across entire tree, for all changes from character u to v 
*********************************************************************************************/
VVdouble computePosteriorExpectationOfChange::computePosteriorAcrossTree(
	simulateJumps &sim, //input given from simulation studies
	const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch)
{
	int alphabetSize = _sp->alphabetSize();
	// N: resized before 
	//probsForBranch.resize(numNodes);
	//for (int n=0;n<numNodes;++n)
	//	resizeMatrix(probsForBranch[n],alphabetSize,alphabetSize);

	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				probsForBranch[mynode->id()][fromState][toState]= computePosteriorOfChangePerBranch(sim,posteriorProbsGivenTerminals,mynode,fromState,toState);
				res[fromState][toState] +=probsForBranch[mynode->id()][fromState][toState];
			}
		}
	}
	return res;
}


/********************************************************************************************
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfChange::computePosteriorOfChangePerBranch(simulateJumps &sim, //input given from simulation studies
	const VVVdouble &posteriorProbs,
	tree::nodeP node,
	int fromState, int toState)
{
	int alphabetSize = _sp->alphabetSize();
	MDOUBLE nodeProbability = 0;

	for (int x=0;x<alphabetSize;++x)
	{
		for (int y=0;y<alphabetSize;++y)
		{
			nodeProbability+=sim.getProb(node->name(),x,y,fromState,toState)*posteriorProbs[node->id()][x][y];
			if(node->name()=="A" ){ //// DEBUG  && x==0 && y==0
				LOGnOUT(9,<<"Sim nodeProbability, given start "<<x<<" end "<<y<<" fromState "<<fromState<<" toState "<<toState<<" = "
					<<sim.getProb(node->name(),x,y,fromState,toState)
					<<endl);
			}
		}
	}
	if(node->name()=="A"){ //// DEBUG 
		LOGnOUT(8,<<"Sim nodeProbability fromState "<<fromState<<" toState "<<toState<<" = "<<nodeProbability<<endl);
	}
	return nodeProbability;
}


/********************************************************************************************
Posterior of observing a certain state change along a branch (joint):
P(Node=x,Father=y|D) = P(D,Node=x,Father=y)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
void computePosteriorExpectationOfChange::computePosteriorOfChangeGivenTerminals(VVVdouble &posteriorPerNodePer2States, int pos){
	int numNodes = _tr.getNodesNum();
	int alphabetSize = _sp->alphabetSize();
	posteriorPerNodePer2States.resize(numNodes);
	for (int n=0;n<posteriorPerNodePer2States.size();++n)
		resizeMatrix(posteriorPerNodePer2States[n],alphabetSize,alphabetSize);
	suffStatGlobalHomPos sscUp;
	suffStatGlobalGamPos sscDownNonRev;	// The "Gam" is used for the letter at father - sscGivenRoot
	//suffStatGlobalHomPos sscDown;
	sscUp.allocatePlace(numNodes,alphabetSize);
	computePijHom pi;
	pi.fillPij(_tr,*_sp); 

	computeUpAlg comp_Up;
	computeDownAlg comp_Down;
	comp_Up.fillComputeUp(_tr,_sc,pos,pi,sscUp);
	//if(!_sp->isReversible())
		comp_Down.fillComputeDownNonReversible(_tr,_sc,pos,pi,sscDownNonRev,sscUp);
	//else
		//comp_Down.fillComputeDown(_tr,_sc,pos,pi,sscDown,sscUp);
		//errorMsg::reportError("error @computePosteriorExpectationOfChange::computePosteriorOfChangeGivenTerminals - Reversible not implemented\n");

	treeIterTopDownConst tIt(_tr);
	doubleRep Ldata = likelihoodComputation::getLofPos(pos,_tr,_sc,pi,*_sp);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int sonState = 0; sonState<alphabetSize; ++sonState){
			for (int fatherState = 0; fatherState<alphabetSize; ++fatherState){
				posteriorPerNodePer2States[mynode->id()][fatherState][sonState]= computePosterioGivenTerminalsPerBranch(mynode->id(),sonState,fatherState,sscUp,sscDownNonRev, pi,Ldata,mynode->name());
				LOGnOUT(7,<<"mynode"<<"\t"<<"fatherState"<<"\t"<<"sonState"<<"\t"<<"posterior(joint)"<<endl);
				LOGnOUT(7,<<mynode->name()<<"\t"<<fatherState<<"\t"<<sonState<<"\t"<<posteriorPerNodePer2States[mynode->id()][fatherState][sonState]<<endl);
			}
		}
	}
}
/********************************************************************************************
Posterior of observing a certain state change along a branch:
P(Node=sonState,Father=fatherState|D) = P(D,Node=sonState,Father=fatherState)/P(D)
usage: posteriorPerNodePer2States[mynode->id()][fatherState][sonState]
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfChange::computePosterioGivenTerminalsPerBranch
	(int nodeId,int sonState, int fatherState,suffStatGlobalHomPos &sscUp,
	suffStatGlobalGamPos &sscDown,computePijHom &pi, doubleRep &Ldata, const string nodeName)
{
	doubleRep res=0.0;
	doubleRep resDXY, Down, Up;
	MDOUBLE pij;
	for (int stateAtRoot = 0; stateAtRoot<_sp->alphabetSize(); ++stateAtRoot){
		Down = sscDown.get(stateAtRoot,nodeId,fatherState);
		Up = sscUp.get(nodeId,sonState);
		pij = pi.getPij(nodeId,fatherState,sonState);

		res+=(_sp->freq(stateAtRoot)*
			Down*
			Up*
			pij);
	}
	resDXY = res;
	res/=Ldata;
	if(gainLossOptions::_printDEBUGinfo)
		LOG(3,<<nodeName<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<Ldata<<" prob "<<res<<endl);

	if (res >1+1e-4){
		LOGnOUT(2,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<Ldata<<" prob "<<res<<endl);
		res = 1;
	}
	if (res<-1e-4){
		LOGnOUT(2,<<nodeId<<" son "<<sonState<<" Down "<<Down<<" father "<<fatherState<<" Up "<<Up<<" pij "<<pij<<" resDXY "<<resDXY<<" LLData "<<Ldata<<" prob "<<res<<endl);
		res = 0;
	}
	if ((res > 1 +0.01) || (res< -0.01)){
		string err = "Error in computePosteriorExpectationOfChange::computePosterioGivenTerminalsPerBranch, non probability value ";
		err+=double2string(convert(res));
		err+=" at node ";
		err+=int2string(nodeId);
		err+=  " sonState ";
		err+= int2string(sonState);
		err+= " fatherState ";
		err+= int2string(fatherState);
		errorMsg::reportError(err);
	}
	return convert(res);
}






/********************************************************************************************
Suchard - Analytic solution - Expectation
*********************************************************************************************/
/********************************************************************************************
Expectation of number of changes from character u to v --- =
Suchard...
*********************************************************************************************/
VVdouble computePosteriorExpectationOfChange::computeExpectationAcrossTree(
	computeJumps &computeJumpsObj,  // object for Analytical computation
	const VVVdouble &posteriorProbs,
	VVVdouble &expForBranch)	// 2 be filled
{
	int alphabetSize = _sp->alphabetSize();
	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				expForBranch[mynode->id()][fromState][toState] = computeExpectationOfChangePerBranch(computeJumpsObj,posteriorProbs,mynode,fromState,toState);
				res[fromState][toState] +=expForBranch[mynode->id()][fromState][toState];
			}
		}		
	}
	return res;
}
/********************************************************************************************
computeExpectationOfChangePerBranch - Analytic...
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfChange::computeExpectationOfChangePerBranch(
	computeJumps &computeJumpsObj, // object for analytical computation
	const VVVdouble &posteriorProbsGivenTerminals,
	tree::nodeP node,int fromState, int toState)
{
	MDOUBLE nodeExpectation = 0;
	//MDOUBLE expGivenStart0nodeA = 0;	// DEBUG
	//LOG(6,<<"\n analytic "<<endl);


	if(node->dis2father()<0) // ROOT
		return nodeExpectation;
	int alphabetSize = _sp->alphabetSize();
	for (int x = 0; x<alphabetSize; ++x){
		for (int y = 0; y<alphabetSize; ++y){
			nodeExpectation+=(posteriorProbsGivenTerminals[node->id()][x][y]*
				computeJumpsObj.getExpectation(node->dis2father(),x,y,fromState,toState));
			if(node->name()=="A" && x==0){	//// DEBUG
				LOG(9,<<"node "<<node->name()<<" All transitions "<<" given "<<x<<" and "<<y
					<<" m= "<< computeJumpsObj.getTotalExpectation(node->dis2father(),x,y)
					<<" m/pij= "<< computeJumpsObj.getTotalExpectation(node->dis2father(),x,y)/static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->Pij_t(x,y,node->dis2father())
					<<" Pij="  << static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->Pij_t(x,y,node->dis2father())<<endl);
			}
			LOG(7,<<"node "<<node->name()<<" given St="<<x<<" and Ed="<<y<<" from "<<fromState<<" to "<<toState
				<<" JointPost= "<< posteriorProbsGivenTerminals[node->id()][x][y]
				<<" Exp= "<< computeJumpsObj.getExpectation(node->dis2father(),x,y,fromState,toState)
				<<endl);
		}
	}
	if(node->name()=="A"){ //// DEBUG 
		LOGnOUT(9,<<"ComG node A All "<<node->dis2father()<<" exp="<<computeJumpsObj.gFunc_dr(node->dis2father(),0)<<endl);
		LOGnOUT(9,<<"nodeExpectation fromState "<<fromState<<" toState "<<toState<<" = "<<nodeExpectation<<endl);
	}
	return nodeExpectation;
}
/********************************************************************************************
Suchard - Analytic solution - Probability
*********************************************************************************************/
/********************************************************************************************
computePosteriorAcrossTree...
*********************************************************************************************/
VVdouble computePosteriorExpectationOfChange::computePosteriorAcrossTree(
	computeJumps &computeJumpsObj, //input given from simulation studies
	const VVVdouble &posteriorProbsGivenTerminals,VVVdouble &probsForBranch)
{
	int alphabetSize = _sp->alphabetSize();
	// N: resized before 
	//probsForBranch.resize(numNodes);
	//for (int n=0;n<numNodes;++n)
	//	resizeMatrix(probsForBranch[n],alphabetSize,alphabetSize);

	VVdouble res;
	resizeMatrix(res,alphabetSize,alphabetSize);
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int fromState=0;fromState<alphabetSize;++fromState)
		{
			for (int toState=0;toState<alphabetSize;++toState)
			{
				if (fromState==toState) 
					continue;
				probsForBranch[mynode->id()][fromState][toState]= computePosteriorOfChangePerBranch(computeJumpsObj,posteriorProbsGivenTerminals,mynode,fromState,toState);
				res[fromState][toState] +=probsForBranch[mynode->id()][fromState][toState];
			}
		}
	}
	return res;
}
/********************************************************************************************
Suchard
*********************************************************************************************/
MDOUBLE computePosteriorExpectationOfChange::computePosteriorOfChangePerBranch(
	computeJumps &computeJumpsObj, //input given from simulation studies
	const VVVdouble &posteriorProbs,
	tree::nodeP node,int fromState, int toState)
{
	int alphabetSize = _sp->alphabetSize();
	MDOUBLE nodeProbability = 0;

	for (int x=0;x<alphabetSize;++x)
	{
		for (int y=0;y<alphabetSize;++y)
		{
			nodeProbability+=computeJumpsObj.getProb(node->dis2father(),x,y,fromState,toState)*posteriorProbs[node->id()][x][y];
			if(node->name()=="A" ){ //// DEBUG && x==0 && y==0
				LOGnOUT(9,<<"Anal nodeProbability, given start "<<x<<" end "<<y<<" fromState "<<fromState<<" toState "<<toState<<" = "
					<<computeJumpsObj.getProb(node->dis2father(),x,y,fromState,toState)
					<<endl);
			}
		}
	}
	if(node->name()=="A"){ //// DEBUG 
		LOGnOUT(8,<<"Anal nodeProbability fromState "<<fromState<<" toState "<<toState<<" = "<<nodeProbability<<endl);
	}

	return nodeProbability;
}




