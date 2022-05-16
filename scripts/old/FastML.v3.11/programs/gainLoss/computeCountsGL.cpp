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
#include "computeCountsGL.h"
#include "gainLossUtils.h"
#include "gainLossAlphabet.h"
#include "computePosteriorExpectationOfChange.h"
#include "computeJumps.h"



/********************************************************************************************
computeCountsGL
*********************************************************************************************/
computeCountsGL::computeCountsGL(sequenceContainer& sc, tree& tr, stochasticProcess* sp,  string& outDir, VVdouble& logLpostPerCatPerPos, MDOUBLE distanceFromNearestOTUForRecent, bool isSilent):
_tr(tr),_sp(sp),_sc(sc),_outDir(outDir),_postProbPerCatPerPos(logLpostPerCatPerPos),_distanceFromNearestOTUForRecent(distanceFromNearestOTUForRecent), _isSilent(isSilent)
{
	_alphabetSize = _sp->alphabetSize();
}
computeCountsGL::computeCountsGL(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution* gainDist, distribution* lossDist, string& outDir, VVVdouble& logLpostPerSpPerCatPerPos, MDOUBLE distanceFromNearestOTUForRecent, bool isSilent):
_tr(tr),_spVVec(spVVec), _gainDist(gainDist), _lossDist(lossDist),_sc(sc),_outDir(outDir),_postProbPerSpPerCatPerPos(logLpostPerSpPerCatPerPos),_distanceFromNearestOTUForRecent(distanceFromNearestOTUForRecent), _isSilent(isSilent)
{
	_alphabetSize = _spVVec[0][0]->alphabetSize();
}

computeCountsGL::~computeCountsGL(){
	//clearVVVV(_jointProb_PosNodeXY);
}

computeCountsGL& computeCountsGL::operator=(const computeCountsGL &other){
	if (this != &other) {              // Check for self-assignment
	}
	return *this;
}


/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::run()
{
	LOGnOUT(4, <<endl<<"Computation stochastic mapping"<<endl);
	time_t t1,t2;
	time(&t1);

	_expV01.resize(_sc.seqLen());
	_expV10.resize(_sc.seqLen());
	_probV01.resize(_sc.seqLen());
	_probV10.resize(_sc.seqLen());
	resizeVVV(_sc.seqLen(),_alphabetSize,_alphabetSize,_expV);
	resizeVVV(_sc.seqLen(),_alphabetSize,_alphabetSize,_probV);
	resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_jointProb_PosNodeXY);
	resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_probChanges_PosNodeXY);
	resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_alphabetSize,_alphabetSize,_expChanges_PosNodeXY);

	if(!gainLossOptions::_gainLossDist){
		computePosteriorOfChangeGivenTerminalsPerCat();
	}
	else
		computePosteriorOfChangeGivenTerminalsPerSpPerCat();	// GLM - multiple SPs

	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}
/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::computePosteriorOfChangeGivenTerminalsPerCat()
{
	// Per RateCategory -- All the computations is done while looping over rate categories
	for (int rateIndex=0 ; rateIndex< _sp->categories(); ++rateIndex)
	{
		tree copy_et = _tr;
		MDOUBLE rateVal = _sp->rates(rateIndex);
		MDOUBLE  minimumRate = 0.000000001;	//0.0000001
		MDOUBLE rate2multiply = max(rateVal,minimumRate);
		if(rateVal<minimumRate){
			LOGnOUT(4, <<" >>> NOTE: the rate category "<<rateVal<<" is too low for computePosteriorExpectationOfChangePerSite"<<endl);	}
		copy_et.multipleAllBranchesByFactor(rate2multiply);
		if(!_isSilent) 
			LOGnOUT(4, <<"Computation performed analytically for rate "<<rate2multiply<<endl);		
		//gainLossAlphabet alph;	// needed for Alphabet size
		//int alphSize = ;
		simulateJumps simPerRateCategory(copy_et,*_sp,_alphabetSize);
		// Per POS		
		for (int pos = 0; pos <_sc.seqLen(); ++pos)
		{
			LOG(9,<<"pos "<<pos+1<<endl);
			// I) computeJoint "computePosteriorOfChangeGivenTerminals" (posteriorPerNodePer2States[mynode->id()][fatherState][sonState])
			VVVdouble posteriorsGivenTerminalsPerRateCategoryPerPos;			
			computePosteriorExpectationOfChange cpecPerRateCategoryPerPos(copy_et,_sc,_sp);	// Per POS,CAT
			cpecPerRateCategoryPerPos.computePosteriorOfChangeGivenTerminals(posteriorsGivenTerminalsPerRateCategoryPerPos,pos);

			// Exp vars - allocate
			VVVdouble expChangesForBranchPerRateCategoryPerPos;	// Sim+Exp
			resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),expChangesForBranchPerRateCategoryPerPos);			
			VVdouble expVV;	// Per POS

			// Prob vars - allocate
			VVVdouble probChangesForBranchPerRateCategoryPerPos;	// Sim+Prob
			resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),probChangesForBranchPerRateCategoryPerPos);
			VVdouble probVV;

			////////////////////////////////////////////////////////////////////////// Analytical
			if(gainLossOptions::_isAnaliticComputeJumps){
				MDOUBLE Lambda1 = static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getMu1();
				MDOUBLE Lambda2 = static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getMu2();
				if(Lambda1 == Lambda2)
					Lambda2 += 0.000000000000001; //NOTE: this is required for analyticComputeSimulateion, to avoid Lambda1=Lambda2
				computeJumps computeJumpsObj(Lambda1,Lambda2);
				
				// II) PostExp:  take in account both: 1) Analytical equations 2) posteriorsGivenTerminal 
				VVVdouble expChangesForBranchPerRateCategoryPerPosAnal;	
				resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),expChangesForBranchPerRateCategoryPerPosAnal);
				VVdouble expVVAnal = cpecPerRateCategoryPerPos.computeExpectationAcrossTree(computeJumpsObj,posteriorsGivenTerminalsPerRateCategoryPerPos,expChangesForBranchPerRateCategoryPerPosAnal);
				expVV = expVVAnal;
				expChangesForBranchPerRateCategoryPerPos = expChangesForBranchPerRateCategoryPerPosAnal;

				// III) PostProbChange: take in account both: 1) Analytical equations 2) posteriorsGivenTerminal 
				VVVdouble probChangesForBranchPerRateCategoryPerPosAnal;	
				resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),probChangesForBranchPerRateCategoryPerPosAnal);
				VVdouble probVVAnal = cpecPerRateCategoryPerPos.computePosteriorAcrossTree(computeJumpsObj,posteriorsGivenTerminalsPerRateCategoryPerPos,probChangesForBranchPerRateCategoryPerPosAnal);
				probVV = probVVAnal;
				probChangesForBranchPerRateCategoryPerPos = probChangesForBranchPerRateCategoryPerPosAnal;
			}
			else{
				if(!_isSilent) 
					LOGnOUT(4, <<"running "<<gainLossOptions::_numOfSimulationsForPotExp<<" simulations for rate "<<rate2multiply<<endl);
				simPerRateCategory.runSimulation(gainLossOptions::_numOfSimulationsForPotExp);	
				if(!_isSilent )
					LOGnOUT(4,<<"finished simulations"<<endl);

				// II) PostExp:  take in account both: 1) simulations 2) posteriorsGivenTerminal 
				expVV = cpecPerRateCategoryPerPos.computeExpectationAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
					expChangesForBranchPerRateCategoryPerPos);	
				// III) PostProbChange: take in account both: 1) simulations 2) posteriorsGivenTerminal 
				probVV = cpecPerRateCategoryPerPos.computePosteriorAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
					probChangesForBranchPerRateCategoryPerPos);

			}
			//////////////////////////////////////////////////////////////////////////

			MDOUBLE exp01 = expVV[0][1];
			MDOUBLE exp10 = expVV[1][0];
			_expV01[pos]+=exp01*_postProbPerCatPerPos[rateIndex][pos];
			_expV10[pos]+=exp10*_postProbPerCatPerPos[rateIndex][pos];
			_expV[pos][0][1]+=exp01*_postProbPerCatPerPos[rateIndex][pos];
			_expV[pos][1][0]+=exp10*_postProbPerCatPerPos[rateIndex][pos];
			
			MDOUBLE prob01 = probVV[0][1];
			MDOUBLE prob10 = probVV[1][0];
			_probV01[pos]+=prob01*_postProbPerCatPerPos[rateIndex][pos];
			_probV10[pos]+=prob10*_postProbPerCatPerPos[rateIndex][pos];
			_probV[pos][0][1]+=prob01*_postProbPerCatPerPos[rateIndex][pos];
			_probV[pos][1][0]+=prob10*_postProbPerCatPerPos[rateIndex][pos];


//	Store all information PerCat,PerPOS
			for(int i=0;i<_probChanges_PosNodeXY[pos].size();++i){	// nodeId
				for(int j=0;j<_probChanges_PosNodeXY[pos][i].size();++j){	// fatherState
					for(int k=0;k<_probChanges_PosNodeXY[pos][i][j].size();++k){	// sonState
						_probChanges_PosNodeXY[pos][i][j][k] += probChangesForBranchPerRateCategoryPerPos[i][j][k]*_postProbPerCatPerPos[rateIndex][pos];
						_expChanges_PosNodeXY[pos][i][j][k] += expChangesForBranchPerRateCategoryPerPos[i][j][k]*_postProbPerCatPerPos[rateIndex][pos];
						_jointProb_PosNodeXY[pos][i][j][k] += posteriorsGivenTerminalsPerRateCategoryPerPos[i][j][k]*_postProbPerCatPerPos[rateIndex][pos];
					}
				}
			}
		}
	}
}

/********************************************************************************************
spVV
*********************************************************************************************/
void computeCountsGL::computePosteriorOfChangeGivenTerminalsPerSpPerCat()
{	
	int numOfSPs = _gainDist->categories()*_lossDist->categories();

	// per Sp
	for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
		int gainIndex =fromIndex2gainIndex(spIndex,_gainDist->categories(),_lossDist->categories());
		int lossIndex =fromIndex2lossIndex(spIndex,_gainDist->categories(),_lossDist->categories());
		_sp = _spVVec[gainIndex][lossIndex];
		if(!_isSilent){
			LOGnOUT(4,<<"computePosteriorOfChangeGivenTerminalsPerSpPerCat with sp:\n Gain= "<<static_cast<gainLossModel*>((*_sp).getPijAccelerator()->getReplacementModel())->getMu1() <<endl);
			if(!gainLossOptions::_isReversible)LOGnOUT(4,<<" Loss= "<<static_cast<gainLossModelNonReversible*>((*_sp).getPijAccelerator()->getReplacementModel())->getMu2() <<endl);
		}
		// Per RateCategory -- All the computations is done while looping over rate categories
		int numOfRateCategories = _spVVec[gainIndex][lossIndex]->categories();	// same for all SPs
		for (int rateIndex=0 ; rateIndex< numOfRateCategories; ++rateIndex)
		{
			tree copy_et = _tr;
			MDOUBLE rateVal = _sp->rates(rateIndex);
			MDOUBLE  minimumRate = 0.000000001;	//0.0000001
			MDOUBLE rate2multiply = max(rateVal,minimumRate);
			if(rateVal<minimumRate){
				LOGnOUT(4, <<" >>> NOTE: the rate category "<<rateVal<<" is too low for computePosteriorExpectationOfChangePerSite"<<endl);	}
			copy_et.multipleAllBranchesByFactor(rate2multiply);
			
			
			//if(!_isSilent) 
			//	LOGnOUT(4, <<"running "<<gainLossOptions::_numOfSimulationsForPotExp<<" simulations for rate "<<rate2multiply<<endl);
			////gainLossAlphabet alph;	// needed for Alphabet size
			//simulateJumps simPerRateCategory(copy_et,*_sp,_alphabetSize);	
			//simPerRateCategory.runSimulation(gainLossOptions::_numOfSimulationsForPotExp);
			//if(!_isSilent) 
			//	LOGnOUT(4,<<"finished simulations"<<endl);
			
			
			simulateJumps simPerRateCategory(copy_et,*_sp,_alphabetSize);	
			// Per POS		
			for (int pos = 0; pos <_sc.seqLen(); ++pos)
			{
				LOG(7,<<"pos "<<pos+1<<endl);
				// I) computeJoint "computePosteriorOfChangeGivenTerminals" (posteriorPerNodePer2States[mynode->id()][fatherState][sonState])
				VVVdouble posteriorsGivenTerminalsPerRateCategoryPerPos;			
				computePosteriorExpectationOfChange cpecPerRateCategoryPerPos(copy_et,_sc,_sp);	// Per POS,CAT
				cpecPerRateCategoryPerPos.computePosteriorOfChangeGivenTerminals(posteriorsGivenTerminalsPerRateCategoryPerPos,pos);

				// Exp vars - allocate
				VVVdouble expChangesForBranchPerRateCategoryPerPos;	// Sim+Exp
				resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),expChangesForBranchPerRateCategoryPerPos);			
				VVdouble expVV;	// Per POS

				// Prob vars - allocate
				VVVdouble probChangesForBranchPerRateCategoryPerPos;	// Sim+Prob
				resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),probChangesForBranchPerRateCategoryPerPos);
				VVdouble probVV;

				////////////////////////////////////////////////////////////////////////// Analytical
				if(gainLossOptions::_isAnaliticComputeJumps){
					MDOUBLE Lambda1 = static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getMu1();
					MDOUBLE Lambda2 = static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getMu2();
					computeJumps computeJumpsObj(Lambda1,Lambda2);

					// II) PostExp:  take in account both: 1) Analytical equations 2) posteriorsGivenTerminal 
					VVVdouble expChangesForBranchPerRateCategoryPerPosAnal;	
					resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),expChangesForBranchPerRateCategoryPerPosAnal);
					VVdouble expVVAnal = cpecPerRateCategoryPerPos.computeExpectationAcrossTree(computeJumpsObj,posteriorsGivenTerminalsPerRateCategoryPerPos,expChangesForBranchPerRateCategoryPerPosAnal);
					expVV = expVVAnal;
					expChangesForBranchPerRateCategoryPerPos = expChangesForBranchPerRateCategoryPerPosAnal;

					// III) PostProbChange: take in account both: 1) Analytical equations 2) posteriorsGivenTerminal 
					VVVdouble probChangesForBranchPerRateCategoryPerPosAnal;	
					resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),probChangesForBranchPerRateCategoryPerPosAnal);
					VVdouble probVVAnal = cpecPerRateCategoryPerPos.computePosteriorAcrossTree(computeJumpsObj,posteriorsGivenTerminalsPerRateCategoryPerPos,probChangesForBranchPerRateCategoryPerPosAnal);
					probVV = probVVAnal;
					probChangesForBranchPerRateCategoryPerPos = probChangesForBranchPerRateCategoryPerPosAnal;
				}
				else{
					if(!_isSilent) 
						LOGnOUT(4, <<"running "<<gainLossOptions::_numOfSimulationsForPotExp<<" simulations for rate "<<rate2multiply<<endl);
					simPerRateCategory.runSimulation(gainLossOptions::_numOfSimulationsForPotExp);	
					if(!_isSilent )
						LOGnOUT(4,<<"finished simulations"<<endl);

					// II) PostExp:  take in account both: 1) simulations 2) posteriorsGivenTerminal 
					expVV = cpecPerRateCategoryPerPos.computeExpectationAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
						expChangesForBranchPerRateCategoryPerPos);	
					// III) PostProbChange: take in account both: 1) simulations 2) posteriorsGivenTerminal 
					probVV = cpecPerRateCategoryPerPos.computePosteriorAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
						probChangesForBranchPerRateCategoryPerPos);

				}
				//////////////////////////////////////////////////////////////////////////

				

				//// II) Exp - take in account both: 1) simulations 2) posteriorsGivenTerminal
				//VVVdouble expChangesForBranchPerRateCategoryPerPos;	// Sim+Exp
				//resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),expChangesForBranchPerRateCategoryPerPos);

				//VVdouble expVV = cpecPerRateCategoryPerPos.computeExpectationAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,
				//	expChangesForBranchPerRateCategoryPerPos);	// Per POS
				MDOUBLE exp01 = expVV[0][1];
				MDOUBLE exp10 = expVV[1][0];
				_expV01[pos]+=exp01*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_expV10[pos]+=exp10*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_expV[pos][0][1]+=exp01*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_expV[pos][1][0]+=exp10*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];

				//// III) Sim - take in account both: 1) simulations 2) posteriorsGivenTerminal
				//VVVdouble probChangesForBranchPerRateCategoryPerPos;	// Sim+Prob
				//resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),probChangesForBranchPerRateCategoryPerPos);
				//VVdouble probVV = cpecPerRateCategoryPerPos.computePosteriorAcrossTree(simPerRateCategory,posteriorsGivenTerminalsPerRateCategoryPerPos,probChangesForBranchPerRateCategoryPerPos);
				MDOUBLE prob01 = probVV[0][1];
				MDOUBLE prob10 = probVV[1][0];
				_probV01[pos]+=prob01*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_probV10[pos]+=prob10*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_probV[pos][0][1]+=prob01*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
				_probV[pos][1][0]+=prob10*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];

				//	Store all information PerCat,PerPOS
				for(int i=0;i<_probChanges_PosNodeXY[pos].size();++i){	// nodeId
					for(int j=0;j<_probChanges_PosNodeXY[pos][i].size();++j){	// fatherState
						for(int k=0;k<_probChanges_PosNodeXY[pos][i][j].size();++k){	// sonState
							_jointProb_PosNodeXY[pos][i][j][k] += posteriorsGivenTerminalsPerRateCategoryPerPos[i][j][k]*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
							_probChanges_PosNodeXY[pos][i][j][k] += probChangesForBranchPerRateCategoryPerPos[i][j][k]*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
							_expChanges_PosNodeXY[pos][i][j][k] += expChangesForBranchPerRateCategoryPerPos[i][j][k]*_postProbPerSpPerCatPerPos[spIndex][rateIndex][pos];
						}
					}
				}
			}
			// Per POS
		}
		// per rateCat
	}
	// Per Sp
}



/********************************************************************************************
printProbExp()
print perPos (over all branches)
use the members _expV01, _expV10 for basic 
*********************************************************************************************/
void computeCountsGL::printProbExp()
{

	string posteriorExpectationOfChangeString = _outDir + "//" + "PosteriorExpectationOfChange.txt";
	ofstream posteriorExpectationStream(posteriorExpectationOfChangeString.c_str());
	posteriorExpectationStream.precision(PRECISION);
	string posteriorProbabilityOfChangeString = _outDir + "//" + "PosteriorProbabilityOfChange.txt";
	ofstream posteriorProbabilityStream(posteriorProbabilityOfChangeString.c_str());
	posteriorProbabilityStream.precision(PRECISION);

	posteriorExpectationStream<<"POS"<<"\t"<<"exp01"<<"\t"<<"exp10"<<endl;
	posteriorProbabilityStream<<"POS"<<"\t"<<"prob01"<<"\t"<<"prob10"<<endl;
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		posteriorExpectationStream<<pos+1<<"\t"<<_expV01[pos]<<"\t"<<_expV10[pos]<<endl;
		posteriorProbabilityStream<<pos+1<<"\t"<<_probV01[pos]<<"\t"<<_probV10[pos]<<endl;

	}
}


/********************************************************************************************
printProbabilityPerPosPerBranch 1
produce 2 print files:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeCountsGL::printProbabilityPerPosPerBranch()
{
	string gainLossProbabilityPerPosPerBranch = gainLossOptions::_outDir + "//" + "ProbabilityPerPosPerBranch.txt"; 
	ofstream gainLossProbabilityPerPosPerBranchStream(gainLossProbabilityPerPosPerBranch.c_str());
	gainLossProbabilityPerPosPerBranchStream.precision(PRECISION);

	gainLossProbabilityPerPosPerBranchStream<<"# print values over probCutOff "<<gainLossOptions::_probCutOffPrintEvent<<endl;
	gainLossProbabilityPerPosPerBranchStream<<"G/L"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"probability"<<endl;
	string gainLossCountProbPerPos = _outDir + "//" + "ProbabilityPerPos.txt"; 
	ofstream gainLossCountProbPerPosStream(gainLossCountProbPerPos.c_str());
	gainLossCountProbPerPosStream.precision(PRECISION);

	//gainLossCountProbPerPosStream<<"# print values over probCutOff "<<gainLossOptions::_probCutOffSum<<endl;
	gainLossCountProbPerPosStream<<"POS"<<"\t"<<"prob01"<<"\t"<<"prob10"<<endl;	
	
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		printGainLossProbabilityPerPosPerBranch(pos, gainLossOptions::_probCutOffPrintEvent, _probChanges_PosNodeXY[pos],gainLossProbabilityPerPosPerBranchStream,gainLossCountProbPerPosStream);
	}
}
/********************************************************************************************
printGainLossProbabilityPerPosPerBranch 1.1
*********************************************************************************************/
void computeCountsGL::printGainLossProbabilityPerPosPerBranch(int pos, MDOUBLE probCutOff, VVVdouble& probChanges, ostream& out, ostream& outCount)
{	
	MDOUBLE count01 =0;
	MDOUBLE count10 =0;
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (probChanges[mynode->id()][0][1] >= probCutOff){
			out<<"gain"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][0][1]<<endl;
		}
		count01+= probChanges[mynode->id()][0][1];
		if (probChanges[mynode->id()][1][0] >= probCutOff){
			out<<"loss"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][1][0]<<endl;
		}
		count10+= probChanges[mynode->id()][1][0];
	}
	outCount<<pos+1<<"\t"<<count01<<"\t"<<count10<<endl;
}


/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::produceExpectationPerBranch(){
	resizeVVV(_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),_expChanges_NodeXY);
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		for(int i=0;i<_expChanges_PosNodeXY[pos].size();++i){
			for(int j=0;j<_expChanges_PosNodeXY[pos][i].size();++j){
				for(int k=0;k<_expChanges_PosNodeXY[pos][i][j].size();++k){
					_expChanges_NodeXY[i][j][k] += _expChanges_PosNodeXY[pos][i][j][k];
				}
			}
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::printExpectationPerBranch()
{
	string gainLossExpectationPerBranch = _outDir + "//" + "ExpectationPerBranch.txt"; 
	ofstream gainLossExpectationPerBranchStream(gainLossExpectationPerBranch.c_str());
	gainLossExpectationPerBranchStream.precision(PRECISION);	
	printGainLossExpectationPerBranch(_expChanges_NodeXY,gainLossExpectationPerBranchStream);
}
/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::printGainLossExpectationPerBranch(VVVdouble& expectChanges, ostream& out)
{
	treeIterTopDownConst tIt(_tr);
	out<<"# Gain and Loss"<<"\n";
	out<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"exp01"<<"\t"<<"exp10"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())
			continue;
		out<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<expectChanges[mynode->id()][0][1]<<"\t"<<expectChanges[mynode->id()][1][0]<<endl;
	}
}

/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::updateTreeByGainLossExpectationPerBranch(tree& tr, int from, int to)
{
	tr = _tr;
	treeIterTopDownConst tIt(tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())
			continue;
		mynode->setDisToFather(_expChanges_NodeXY[mynode->id()][from][to]);
	}
}


/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::printTreesWithExpectationValuesAsBP()
{
	// ExpectationPerPosPerBranch - Print Trees
	Vstring Vnames;
	fillVnames(Vnames,_tr);
	createDir(gainLossOptions::_outDir, "TreesWithExpectationValuesAsBP");
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		string strTreeNum = _outDir + "//" + "TreesWithExpectationValuesAsBP" + "//" + "expTree" + int2string(pos+1) + ".ph";
		ofstream tree_out(strTreeNum.c_str());
		tree_out.precision(PRECISION);
		printTreeWithValuesAsBP(tree_out,_tr,Vnames,&_expChanges_PosNodeXY[pos]);
	}
}

/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::printTreesWithProbabilityValuesAsBP()
{
	// ProbabilityPerPosPerBranch - Print Trees
	Vstring Vnames;
	fillVnames(Vnames,_tr);
	createDir(_outDir, "TreesWithProbabilityValuesAsBP");
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		string strTreeNum = _outDir + "//" + "TreesWithProbabilityValuesAsBP"+ "//" + "probTree" + int2string(pos+1) + ".ph";
		ofstream tree_out(strTreeNum.c_str());
		printTreeWithValuesAsBP(tree_out,_tr,Vnames,&_probChanges_PosNodeXY[pos]);
	}
}

/********************************************************************************************
printProbExpPerPosPerBranch 1
produce 2 print files:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeCountsGL::printProbExpPerPosPerBranch(MDOUBLE probCutOff, MDOUBLE countsCutOff)
{
	string gainLossProbExpPerPosPerBranch = _outDir + "//" + "gainLossProbExpPerPosPerBranch.txt"; 
	ofstream gainLossProbExpPerPosPerBranchStream(gainLossProbExpPerPosPerBranch.c_str());
	gainLossProbExpPerPosPerBranchStream.precision(PRECISION);
	gainLossProbExpPerPosPerBranchStream<<"# print values over probCutOff "<<probCutOff<<endl;
	gainLossProbExpPerPosPerBranchStream<<"G/L"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"probability"<<"\t"<<"expectation"<<endl;
	string gainLossProbExpPerPos = _outDir + "//" + "gainLossProbExpCountPerPos.txt"; 
	ofstream gainLossCountProbPerPosStream(gainLossProbExpPerPos.c_str());
	gainLossCountProbPerPosStream.precision(PRECISION);
	gainLossCountProbPerPosStream<<"# print count over countsCutOff "<<countsCutOff<<endl;
	gainLossCountProbPerPosStream<<"POS"<<"\t"<<"prob01"<<"\t"<<"prob10"<<"\t"<<"exp01"<<"\t"<<"exp10"<<"\t"<<"count01"<<"\t"<<"count10"<<endl;
	
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		printGainLossProbExpPerPosPerBranch(pos, probCutOff,countsCutOff, _probChanges_PosNodeXY[pos],_expChanges_PosNodeXY[pos],gainLossProbExpPerPosPerBranchStream,gainLossCountProbPerPosStream);
	}
}
/********************************************************************************************
 PrintExpPerPosPerBranchMatrix (CoMap input)
 NOTE!!! this version only consist of gain or loss values
 Alternatively, (1) abs(gain+loss) (2) gain-loss (3) separate gain and loss matrices
*********************************************************************************************/
void computeCountsGL::printExpPerPosPerBranchMatrix(const int from, const int to)
{
	int numOfpositions = _sc.seqLen();
	int numOfbranches = _tr.getNodesNum()-1; // minus the root node
	
	string expPerPosPerBranchMatrix = _outDir + "//" + "expPerPosPerBranchMatrix."+ int2string(from)+int2string(to)+".txt"; 
	ofstream expPerPosPerBranchMatrixStream(expPerPosPerBranchMatrix.c_str());
	expPerPosPerBranchMatrixStream.precision(6);
	expPerPosPerBranchMatrixStream<<"Name\tLength\tBranches\tMean";
	for (int pos = 0; pos <numOfpositions; ++pos){
		expPerPosPerBranchMatrixStream<<"\tSite"<<pos+1;
	}
	expPerPosPerBranchMatrixStream<<"\n";
	treeIterTopDownConst tIt(_tr);
	int branchNum = 0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())
			continue;
		expPerPosPerBranchMatrixStream<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<branchNum<<"\t"<<_expChanges_NodeXY[mynode->id()][from][to]/numOfbranches;
		for (int pos = 0; pos <numOfpositions; ++pos){
			expPerPosPerBranchMatrixStream<<"\t"<<_expChanges_PosNodeXY[pos][mynode->id()][from][to];
		}		
		expPerPosPerBranchMatrixStream<<"\n";
		++branchNum;
	}
	expPerPosPerBranchMatrixStream.close();
}


///********************************************************************************************
//*********************************************************************************************/
//void computeCountsGL::fillCorrPerSelectedSites(Vdouble& correlationPerPos,VVdouble& expEventsPerPosPerBranch,VVdouble& expEventsPerPosPerBranch_B,const int selectedSite, const bool isPearson){
//	int numOfpositions = expEventsPerPosPerBranch_B.size();
//	//correlationPerPos.resize(numOfpositions);
//	for (int pos = 0; pos <numOfpositions; ++pos){
//		MDOUBLE correlation = 0;
//		if(isMinEQMaxInVector(expEventsPerPosPerBranch[selectedSite]) || isMinEQMaxInVector(expEventsPerPosPerBranch_B[pos]))
//			correlationPerPos[pos]=-99; // can't compute correlation
//		else{
//			if(isPearson)
//				correlation = calcPearsonCorrelation(expEventsPerPosPerBranch[selectedSite], expEventsPerPosPerBranch_B[pos]);
//			else
//				correlation = calcRankCorrelation(expEventsPerPosPerBranch[selectedSite], expEventsPerPosPerBranch_B[pos]);
//				correlationPerPos[pos]=correlation;
//		}
//	}		
//}
/********************************************************************************************
Compute the Pearson / Spearman correlation among sites.
*********************************************************************************************/
//void computeCountsGL::computedCorrelations(const Vint& selectedPositions, const bool isNormalizeForBranch)
//{	
//	int numOfpositions = _sc.seqLen();
//	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node
//
//	//// Mapping vectors
//	LOGnOUT(6, <<"Copy events vectors"<<endl);
//	// Expectation
//	fillMapValPerPosPerBranch(_expPerPosPerBranch01,0,1,_expChanges_PosNodeXY,isNormalizeForBranch);
//	fillMapValPerPosPerBranch(_expPerPosPerBranch10,1,0,_expChanges_PosNodeXY,isNormalizeForBranch);
//	_expPerPosPerBranch = _expPerPosPerBranch01; // gain and loss appended (double size vector)
//	appendVectors(_expPerPosPerBranch, _expPerPosPerBranch10);
//
//	//// correlation vectors, filled below
//	LOGnOUT(6, <<"Resize correlation vectors vectors"<<endl);
//	resizeMatrix(_correlationPerSitePerPosGainGainSpearman, selectedPositions.size(), numOfpositions);
//	resizeMatrix(_correlationPerSitePerPosLossLossSpearman, selectedPositions.size(), numOfpositions);
//	resizeMatrix(_correlationPerSitePerPosBothSpearman, selectedPositions.size(), numOfpositions);
//
//	resizeMatrix(_correlationPerSitePerPosGainGainPearson, selectedPositions.size(), numOfpositions);
//	resizeMatrix(_correlationPerSitePerPosLossLossPearson, selectedPositions.size(), numOfpositions);
//	resizeMatrix(_correlationPerSitePerPosBothPearson, selectedPositions.size(), numOfpositions);
//
//	for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
//		int selectedSite = selectedPositions[selectedSiteIndex];
//		LOGnOUT(6, <<"Compute pearson for G-G, L-L, both site"<<selectedSiteIndex<<endl);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosGainGainPearson[selectedSiteIndex],_expPerPosPerBranch01,selectedSite,true);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosLossLossPearson[selectedSiteIndex],_expPerPosPerBranch10,selectedSite,true);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosBothPearson[selectedSiteIndex],_expPerPosPerBranch,selectedSite,true);
//
//		LOGnOUT(6, <<"Compute spearman for G-G, L-L site"<<selectedSiteIndex<<endl);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosGainGainSpearman[selectedSiteIndex],_expPerPosPerBranch01,selectedSite,false);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosLossLossSpearman[selectedSiteIndex],_expPerPosPerBranch10,selectedSite,false);
//		fillCorrPerSelectedSites(_correlationPerSitePerPosBothSpearman[selectedSiteIndex],_expPerPosPerBranch,selectedSite,false);
//	}	
//}


/********************************************************************************************
PrintExpPerPosPerBranchMatrix (CoMap input)
NOTE!!! this version only consist of gain or loss values
Alternatively, (1) abs(gain+loss) (2) gain-loss (3) separate gain and loss matrices
*********************************************************************************************/
//void computeCountsGL::printComputedCorrelations(const Vint& selectedPositions, const bool isNormalizeForBranch, const bool correlationForZscore)
//{	
//	bool isTransform = false;
//	bool isMinForPrint = true;
//	bool isPearson = false;
//	int precisionCorr = 8;
//	MDOUBLE minForPrint = 0.1; // max =1
//
//	int numOfpositions = _sc.seqLen();
//	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node
//	
//	//// Mapping vectors
//	LOGnOUT(6, <<"Copy events vectors"<<endl);
//
//	//////////////////////////////////////////////////////////////////////////	 
//	if(!gainLossOptions::_printComputedCorrelationsAllSites){
//		for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
//			int selectedSite = selectedPositions[selectedSiteIndex];
//
//			MDOUBLE meanCorrBoth = computeAverage(_correlationPerSitePerPosBothPearson[selectedSiteIndex]);
//			MDOUBLE stdCorrBoth = computeStd(_correlationPerSitePerPosBothPearson[selectedSiteIndex]);
//			MDOUBLE meanCorrGainGain = computeAverage(_correlationPerSitePerPosGainGainPearson[selectedSiteIndex]);
//			MDOUBLE stdCorrGainGain = computeStd(_correlationPerSitePerPosGainGainPearson[selectedSiteIndex]);
//			MDOUBLE meanCorrLossLoss = computeAverage(_correlationPerSitePerPosLossLossPearson[selectedSiteIndex]);
//			MDOUBLE stdCorrLossLoss = computeStd(_correlationPerSitePerPosLossLossPearson[selectedSiteIndex]);
//
//
//			// for each selectedSite a new file is created
//			LOGnOUT(4, <<"Correlations with site="<<selectedSite<<" With NormalizeForBranch "<<isNormalizeForBranch<<" With correlationForZscore "<<correlationForZscore<<endl);
//			string corrPerSite = _outDir + "//" + "selectedCorr.Site"+ int2string(selectedSite+1)+ ".isNormForBr."+int2string(isNormalizeForBranch)/*+ ".isCorrForZ."+int2string(correlationForZscore)*/+ ".txt";
//			//string corrPerSite = _outDir + "//" + "selectedCorr.Site"+ int2string(selectedSite+1)+".txt";
//
//			ofstream corrPerSiteStream(corrPerSite.c_str());
//			corrPerSiteStream.precision(precisionCorr);
//			corrPerSiteStream<<"# "<<selectedSite+1<<"\n";
//			corrPerSiteStream<<"# Both(gain N loss concat) correlation(Pearson): Ave= "<<meanCorrBoth<<" Std= "<<stdCorrBoth<<"\n";
//			corrPerSiteStream<<"# Gain correlation(Pearson): Ave= "<<meanCorrGainGain<<" Std= "<<stdCorrGainGain<<"\n";
//			corrPerSiteStream<<"# Loss correlation: Ave= "<<meanCorrLossLoss<<" Std= "<<stdCorrLossLoss<<"\n";
//			corrPerSiteStream<<"pos"<<"\t"<<"bothPearson"<<"\t"<<"bothSpearman"<<"\t"<<"ExpGainGainPearson"<<"\t"<<"ExpLossLossPearson"<<"\t"<<"ExpGainGainSpearman"<<"\t"<<"ExpLossLossSpearman"<<"\n";
//
//			for (int pos = 0; pos<numOfpositions; ++pos){
//				if(selectedSite == pos)	// since selectedSite starts from 1
//					continue;
//				bool isPosOneOfSelectedSites = false;
//				if(gainLossOptions::_isIgnoreCorrelationAmongSelectedSites){
//					for (int selectedSiteI = 0; selectedSiteI <selectedPositions.size(); ++selectedSiteI){
//						int selectedS = selectedPositions[selectedSiteI];
//						if(selectedS == pos){
//							isPosOneOfSelectedSites = true;
//							continue;
//						}				
//					}
//					if(isPosOneOfSelectedSites)
//						continue;
//				}
//				corrPerSiteStream<<pos+1
//					<<"\t"<<_correlationPerSitePerPosBothPearson[selectedSiteIndex][pos]<<"\t"<<_correlationPerSitePerPosBothSpearman[selectedSiteIndex][pos]
//					<<"\t"<<_correlationPerSitePerPosGainGainPearson[selectedSiteIndex][pos]<<"\t"<<_correlationPerSitePerPosLossLossPearson[selectedSiteIndex][pos]
//					<<"\t"<<_correlationPerSitePerPosGainGainSpearman[selectedSiteIndex][pos]<<"\t"<<_correlationPerSitePerPosLossLossSpearman[selectedSiteIndex][pos]<<"\n";
//			}
//		}
//	}
//	//////////////////////////////////////////////////////////////////////////	All-against-all different format
//	else{
//		string corrAllSites = _outDir + "//" + "allCorrelations.isNormForBr."+int2string(isNormalizeForBranch)+ ".isCorrForZ."+int2string(correlationForZscore)+ ".txt";
//		ofstream* corrAllStream_p;
//		corrAllStream_p = new ofstream(corrAllSites.c_str());
//		corrAllStream_p->precision(precisionCorr);
//		*corrAllStream_p<<"#COGA"<<"\t"<<"COGB"<<"\t"<<"posGainGain"<<"\t"<<"posLossLoss"<<"\t"<<"negGainGain"<<"\t"<<"negLossLoss"<<"\n";
//		for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
//			int selectedSite = selectedPositions[selectedSiteIndex];
//
//			MDOUBLE meanCorrGainGain = computeAverage(_correlationPerSitePerPosGainGainPearson[selectedSiteIndex]);
//			MDOUBLE stdCorrGainGain = computeStd(_correlationPerSitePerPosGainGainPearson[selectedSiteIndex]);
//			MDOUBLE meanCorrLossLoss = computeAverage(_correlationPerSitePerPosLossLossPearson[selectedSiteIndex]);
//			MDOUBLE stdCorrLossLoss = computeStd(_correlationPerSitePerPosLossLossPearson[selectedSiteIndex]);
//
//			for (int pos = 0; pos<numOfpositions; ++pos){
//				if(selectedSite == pos)
//					continue;
//				MDOUBLE correlationGainGain = _correlationPerSitePerPosGainGainPearson[selectedSiteIndex][pos];
//				MDOUBLE correlationLossLoss = _correlationPerSitePerPosLossLossPearson[selectedSiteIndex][pos];
//
//				if(correlationForZscore){
//					correlationGainGain = (correlationGainGain - meanCorrGainGain)/stdCorrGainGain;
//					correlationLossLoss = (correlationLossLoss - meanCorrLossLoss)/stdCorrLossLoss;
//				}
//				if(isMinForPrint && max(abs(correlationGainGain),abs(correlationLossLoss))<minForPrint)
//					continue;
//				MDOUBLE posCorrelationGainGain = (correlationGainGain >=0) ? correlationGainGain*1000-1 : 0;
//				MDOUBLE negCorrelationGainGain = (correlationGainGain < 0) ? correlationGainGain*1000-1 : 0;
//				MDOUBLE posCorrelationLossLoss = (correlationLossLoss >=0) ? correlationLossLoss*1000-1 : 0;
//				MDOUBLE negCorrelationLossLoss = (correlationLossLoss < 0) ? correlationLossLoss*1000-1 : 0;
//				if(isTransform){
//					posCorrelationGainGain = pow(posCorrelationGainGain/10,2)/10;
//					negCorrelationGainGain = pow(negCorrelationGainGain/10,2)/10;
//					posCorrelationLossLoss = pow(posCorrelationLossLoss/10,2)/10;
//					negCorrelationLossLoss = pow(negCorrelationLossLoss/10,2)/10;
//				}
//				*corrAllStream_p<<selectedSiteIndex+1<<"\t"<<pos+1<<"\t"<<(int)posCorrelationGainGain<<"\t"<<(int)posCorrelationLossLoss<<"\t"<<(int)negCorrelationGainGain<<"\t"<<(int)negCorrelationLossLoss<<"\n";
//			}
//		}
//	}
//}
//
///********************************************************************************************
//*********************************************************************************************/
//void computeCountsGL::fillMapValPerPosPerBranch(VVdouble& expEventsPerPosPerBranch,const int from, const int to, VVVVdouble& map_PosNodeXY
//												   ,const bool isNormalizeForBranch, MDOUBLE* cutOff_p){
//
//	int numOfpositions = _sc.seqLen();
//	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node
//
//	expEventsPerPosPerBranch.resize(numOfpositions);
//	treeIterTopDownConst tIt(_tr);
//	for (int pos = 0; pos <numOfpositions; ++pos){
//		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
//		{
//			if(mynode->isRoot())
//				continue;
//			MDOUBLE val = 0;
//			if(isNormalizeForBranch){
//				MDOUBLE normalizationFactor =  _expChanges_NodeXY[mynode->id()][from][to]/numOfbranches; // _expChanges_NodeXY[mynode->id()][from][to]/numOfbranches
//				val = (map_PosNodeXY[pos][mynode->id()][from][to] ) / normalizationFactor;
//			}else{
//				val = map_PosNodeXY[pos][mynode->id()][from][to]; 
//			}
//
//			if(cutOff_p){
//				if(val>= *cutOff_p)
//					expEventsPerPosPerBranch[pos].push_back(1);
//				else
//					expEventsPerPosPerBranch[pos].push_back(0);
//			}
//			else
//				expEventsPerPosPerBranch[pos].push_back(val);			
//		}
//	}
//}




/********************************************************************************************
printGainLossProbExpPerPosPerBranch 1.1
Get pos, and iterate over all branches:
1. print detailed file (out)
2. print summary over all branches (outSum)
*********************************************************************************************/
void computeCountsGL::printGainLossProbExpPerPosPerBranch(int pos, MDOUBLE probCutOff, MDOUBLE countCutOff, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outSum)
{
	MDOUBLE prob01 =0;
	MDOUBLE prob10 =0;
	MDOUBLE exp01 =0;
	MDOUBLE exp10 =0;
	MDOUBLE count01 =0;
	MDOUBLE count10 =0;

	countCutOff = floorf(countCutOff * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example

	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())  continue;
		if (probChanges[mynode->id()][0][1] >= probCutOff || probCutOff == 0) // only per branch print must exceed cutoff
			out<<"gain"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][0][1]<<"\t"<<expChanges[mynode->id()][0][1]<<endl;		
		if (probChanges[mynode->id()][0][1] > countCutOff)
			count01+= 1;
		prob01+= probChanges[mynode->id()][0][1];
		exp01+= expChanges[mynode->id()][0][1];
		if (probChanges[mynode->id()][1][0] >= probCutOff || probCutOff == 0) // only per branch print must exceed cutoff
			out<<"loss"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][1][0]<<"\t"<<expChanges[mynode->id()][1][0]<<endl;		
		if (probChanges[mynode->id()][1][0] > countCutOff)
			count10+= 1;            
		prob10+= probChanges[mynode->id()][1][0];
		exp10+= expChanges[mynode->id()][1][0];
	}
	outSum<<pos+1<<"\t"<<prob01<<"\t"<<prob10<<"\t"<<exp01<<"\t"<<exp10<<"\t"<<count01<<"\t"<<count10<<endl;
}



/********************************************************************************************
FewCutOffs
*********************************************************************************************/
void computeCountsGL::printProbExpPerPosPerBranchFewCutOffs(MDOUBLE probCutOff)
{
	MDOUBLE countCutOff;
	MDOUBLE countCutOffLow = 0.1;
	MDOUBLE countCutOffIncrem = 0.05;
	MDOUBLE countCutOffHigh = 0.9;
	string count01 = "count01_";
	string count10 = "count10_";

    //Math::Round(3.44, 1);

	
	string gainLossProbExpPerPosPerBranch = _outDir + "//" + "gainLossProbExpPerPosPerBranch.txt"; 
	ofstream gainLossProbExpPerPosPerBranchStream(gainLossProbExpPerPosPerBranch.c_str());
	gainLossProbExpPerPosPerBranchStream.precision(PRECISION);
	
	gainLossProbExpPerPosPerBranchStream<<"# print values over probCutOff "<<probCutOff<<endl;
	gainLossProbExpPerPosPerBranchStream<<"G/L"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"probability"<<"\t"<<"expectation"<<endl;
	string gainLossProbExpPerPos = _outDir + "//" + "gainLossProbExpCountPerPos.txt"; 
	ofstream gainLossCountProbPerPosStream(gainLossProbExpPerPos.c_str());
	gainLossCountProbPerPosStream.precision(PRECISION);


	gainLossCountProbPerPosStream<<"# print count over countCutOffLow="<<countCutOffLow<<" to countCutOffHigh="<<countCutOffHigh<<" with increment="<<countCutOffIncrem<<endl;
	gainLossCountProbPerPosStream<<"POS"<<"\t"<<"prob01"<<"\t"<<"prob10"<<"\t"<<"exp01"<<"\t"<<"exp10"<<"\t"<<"prob01_Rec"<<"\t"<<"prob10_Rec"<<"\t"<<"exp01_Rec"<<"\t"<<"exp10_Rec"<<"\t"<<"prob01_Anc"<<"\t"<<"prob10_Anc"<<"\t"<<"exp01_Anc"<<"\t"<<"exp10_Anc"<<"\t";
	// print all cut-offs
	for(countCutOff=countCutOffLow; countCutOff<=countCutOffHigh ;countCutOff+=countCutOffIncrem){
		countCutOff = floorf(countCutOff * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example
        gainLossCountProbPerPosStream<<count01+double2string(countCutOff)<<"\t"<<count10+double2string(countCutOff)<<"\t";
	}
	gainLossCountProbPerPosStream<<endl;

	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		printGainLossProbExpPerPosPerBranchFewCutOffs(pos, probCutOff,countCutOffLow,countCutOffIncrem,countCutOffHigh, _probChanges_PosNodeXY[pos],_expChanges_PosNodeXY[pos],gainLossProbExpPerPosPerBranchStream,gainLossCountProbPerPosStream);
	}
}
/********************************************************************************************
*********************************************************************************************/
void computeCountsGL::printGainLossProbExpPerPosPerBranchFewCutOffs(int pos, MDOUBLE probCutOff, 
				MDOUBLE countCutOffLow,MDOUBLE countCutOffIncrem, MDOUBLE countCutOffHigh, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outSum)
{
	MDOUBLE prob01 =0;
	MDOUBLE prob10 =0;
	MDOUBLE exp01 =0;
	MDOUBLE exp10 =0;

	MDOUBLE prob01_R =0;
	MDOUBLE prob10_R =0;
	MDOUBLE exp01_R =0;
	MDOUBLE exp10_R =0;

	MDOUBLE prob01_Anc =0;
	MDOUBLE prob10_Anc =0;
	MDOUBLE exp01_Anc =0;
	MDOUBLE exp10_Anc =0;

	int FewCutOffsSize = (int)ceil((countCutOffHigh-countCutOffLow)/countCutOffIncrem)+1;	
	Vdouble count01(FewCutOffsSize);
	Vdouble count10(FewCutOffsSize);

	MDOUBLE countCutOff;
	int i;

	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (probChanges[mynode->id()][0][1] >= probCutOff)	// only per branch print must exceed cutoff
			out<<"gain"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][0][1]<<"\t"<<expChanges[mynode->id()][0][1]<<endl;		
		prob01+= probChanges[mynode->id()][0][1];
		exp01+= expChanges[mynode->id()][0][1];
//		if(mynode->isLeaf() || (mynode->getDistance2ROOT()<_distanceFromRootForRecent) ){
		if(mynode->isLeaf() || (mynode->getMinimalDistance2OTU()<_distanceFromNearestOTUForRecent) ){
			prob01_R+= probChanges[mynode->id()][0][1];
			exp01_R+= expChanges[mynode->id()][0][1];
		}
		else{
			prob01_Anc+= probChanges[mynode->id()][0][1];
			exp01_Anc+= expChanges[mynode->id()][0][1];
		}
		i = 0;
		for( countCutOff=countCutOffLow; countCutOff<=countCutOffHigh ; countCutOff+=countCutOffIncrem){
			countCutOff = floorf(countCutOff * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example

			if (probChanges[mynode->id()][0][1] > countCutOff)
				count01[i]+= 1;
			++i;
		}
		if (probChanges[mynode->id()][1][0] >= probCutOff)	// only per branch print must exceed cutoff
			out<<"loss"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<probChanges[mynode->id()][1][0]<<"\t"<<expChanges[mynode->id()][1][0]<<endl;		
		prob10+= probChanges[mynode->id()][1][0];
		exp10+= expChanges[mynode->id()][1][0];
//		if(mynode->isLeaf() || mynode->getDistance2ROOT() < _distanceFromRootForRecent){
		if(mynode->isLeaf() || mynode->getMinimalDistance2OTU() < _distanceFromNearestOTUForRecent){
			prob10_R+= probChanges[mynode->id()][1][0];
			exp10_R+= expChanges[mynode->id()][1][0];
		}
		else{
			prob10_Anc+= probChanges[mynode->id()][1][0];
			exp10_Anc+= expChanges[mynode->id()][1][0];
		}
		i = 0;
		for(countCutOff=countCutOffLow; countCutOff<=countCutOffHigh ; countCutOff+=countCutOffIncrem){
			countCutOff = floorf(countCutOff * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example
			if (probChanges[mynode->id()][1][0] > countCutOff)
				count10[i]+= 1;
			++i;				
		}
	}
	outSum<<pos+1<<"\t"<<prob01<<"\t"<<prob10<<"\t"<<exp01<<"\t"<<exp10
		<<"\t"<<prob01_R<<"\t"<<prob10_R<<"\t"<<exp01_R<<"\t"<<exp10_R
		<<"\t"<<prob01_Anc<<"\t"<<prob10_Anc<<"\t"<<exp01_Anc<<"\t"<<exp10_Anc<<"\t";
	
	// print all cut-offs
	i = 0;
	for(countCutOff=countCutOffLow; countCutOff<=countCutOffHigh ; countCutOff+=countCutOffIncrem){
		countCutOff = floorf(countCutOff * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example
		outSum<<count01[i]<<"\t"<<count10[i]<<"\t";
		++i;
	}
	outSum<<endl;
}

/********************************************************************************************
*********************************************************************************************/
//void computeCountsGL::computeMeanAndSdPerBranch(Vdouble& meanEventsPerBranch01,	Vdouble& meanEventsPerBranch10,	Vdouble& sdEventsPerBranch01,Vdouble& sdEventsPerBranch10){
//	int numOfpositions = _sc.seqLen();
//	Vdouble eventsAllPos01(numOfpositions);
//	Vdouble eventsAllPos10(numOfpositions);
//
//	treeIterTopDownConst tIt(_tr);
//	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
//	{
//		if(mynode->isRoot())
//			continue;
//		for (int pos = 0; pos <numOfpositions; ++pos){
//			eventsAllPos01[pos] = _expChanges_PosNodeXY[pos][mynode->id()][0][1];
//			eventsAllPos10[pos] = _expChanges_PosNodeXY[pos][mynode->id()][1][0];
//		}
//		meanEventsPerBranch01[mynode->id()]= computeAverage(eventsAllPos01);
//		meanEventsPerBranch10[mynode->id()]= computeAverage(eventsAllPos10);
//		sdEventsPerBranch01[mynode->id()] = computeStd(eventsAllPos01);
//		sdEventsPerBranch10[mynode->id()] = computeStd(eventsAllPos10);	
//	}
//}


