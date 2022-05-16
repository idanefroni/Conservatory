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
#include "gainLossOptimizer.h"
#include "bblEMfixRoot.h"
#include "bblEM.h"
#include "bblLS.h"

/********************************************************************************************
gainLossOptimizer
*********************************************************************************************/
gainLossOptimizer::gainLossOptimizer(tree& tr, stochasticProcess* sp, const sequenceContainer &sc,
									 const MDOUBLE epsilonOptimization, const int numIterations,
									 const MDOUBLE epsilonOptimizationModel, const int numIterationsModel,
									 const MDOUBLE epsilonOptimizationBBL, const int numIterationsBBL,
									 Vdouble*  weights,
									 unObservableData* unObservableData_p, bool performOptimizationsBBL, bool isbblLSWhenbblEMdontImprove):
_tr(tr),_sp(sp),_sc(sc),
_epsilonOptimization(epsilonOptimization),_maxNumOfIterations(numIterations),
_epsilonOptimizationModel(epsilonOptimizationModel),_maxNumOfIterationsModel(numIterationsModel),
_epsilonOptimizationBBL(epsilonOptimizationBBL),_maxNumOfIterationsBBL(numIterationsBBL),
_weightsUniqPatterns(weights),
_unObservableData_p(unObservableData_p),_performOptimizationsBBL(performOptimizationsBBL),
_isbblLSWhenbblEMdontImprove(isbblLSWhenbblEMdontImprove)
{
	//gainLossOptions::distributionType rateDistributionType = getRateDistributionType(sp->distr());
	//_weights = gainLossOptions::_weights;			// since - no weights are used over positions
	_isReversible = !dynamic_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel());
	_isSkipBblEM = false; // will change to T if like is not improved by BBL-EM
	_freq.resize(_sc.alphabetSize());
	optimizations();
}


/********************************************************************************************
*********************************************************************************************/
gainLossOptimizer::gainLossOptimizer(tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist,
									 const sequenceContainer &sc,
									 const MDOUBLE epsilonOptimization, const int numIterations,
									 const MDOUBLE epsilonOptimizationModel, const int numIterationsModel,
									 const MDOUBLE epsilonOptimizationBBL, const int numIterationsBBL,
									 Vdouble*  weights,
									 unObservableData* unObservableData_p, bool performOptimizationsBBL, bool isbblLSWhenbblEMdontImprove):
_tr(tr),_spVVec(spVVec),_gainDist(gainDist),_lossDist(lossDist),
_sc(sc),//_spSimple(spSimple), // ignore sent model, make new one
_epsilonOptimization(epsilonOptimization),_maxNumOfIterations(numIterations),
_epsilonOptimizationModel(epsilonOptimizationModel),_maxNumOfIterationsModel(numIterationsModel),
_epsilonOptimizationBBL(epsilonOptimizationBBL),_maxNumOfIterationsBBL(numIterationsBBL),
_weightsUniqPatterns(weights),
_unObservableData_p(unObservableData_p),_performOptimizationsBBL(performOptimizationsBBL),
_isbblLSWhenbblEMdontImprove(isbblLSWhenbblEMdontImprove)
{
	//_sp = _spVVec[0][0];	//used for reference (Alpha and such)
	//_weights = gainLossOptions::_weights;			// since - no weights are used over positions
	_spSimple = NULL;
	_isSkipBblEM = false; // will change to T if like is not improved by BBL-EM
	_freq.resize(_sc.alphabetSize());
	_bestGainBeta = 1;
	_bestLossBeta = 1;
	optimizationsSPvv();
}

/********************************************************************************************
optimizations
*********************************************************************************************/
void gainLossOptimizer::optimizations(){
	time_t t1;
	time(&t1);
	time_t t2;

	LOGnOUT(4,<<"-------------------------------"<<endl
		<<"Starting optimizations: maxNumIterations="<<_maxNumOfIterations<<endl);	

	_bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
	MDOUBLE currBestL=VERYSMALL;
	MDOUBLE previousL;
	bool noLikeImprovmentAtBBL = false; // if BBL did not produce a new tree, end. (no point in another iteration model+BBL)

	bool isSkipParamsOptimization = gainLossOptions::_isSkipFirstParamsOptimization;
	LOGnOUT(3,<<endl<<"#########"<<" optimization starting epsilonCycle="<<_epsilonOptimization<<" maxNumIterations="<<_maxNumOfIterations<<endl);
	LOGnOUT(3,<<"start optimization with L= "<<_bestL<<endl);
	int iter;
	for (iter=1;iter<=_maxNumOfIterations;iter++){
		LOGnOUT(4,<<endl<<"------"<<" Model+BBL iter="<<iter<<endl);
		previousL = _bestL;		// breaking out of loop when no (>epsilon) improvement is made by comparing to previousL
// model optimization
		if(!isSkipParamsOptimization){
			currBestL = optimizeParameters();}
		else{
			LOGnOUT(4,<<"Optimize Params - Skipped"<<endl);}
		
		if (currBestL>_bestL) {
			_bestL = currBestL;
		}
		else if(!isSkipParamsOptimization && currBestL<_bestL){
			LOGnOUT(4,<<" !!! Warning !!!: after model optimization likelihood went down"<< currBestL<<" "<<_bestL<<endl);
		}		
		//_bestL = max(currBestL,_bestL);
		
		isSkipParamsOptimization = false;	// only first iteration skipped
// BBL optimization
		if (gainLossOptions::_performOptimizationsBBL && _performOptimizationsBBL) // we use the && 2 enable optimizationsManyStarts not to perform BBL
		{
			//LOGnOUT(4,<<"Start BBL... with epsilonOptimizationBBL= "<<_epsilonOptimizationBBL<<endl);
			if(gainLossOptions::_performOptimizationsBBLOnlyOnce)	// Next iteration - no BBL
				_performOptimizationsBBL = false;
			currBestL = optimizeBranchLengths(iter);
			if (currBestL>_bestL) {
				_bestL = currBestL;
			}
			else{
				noLikeImprovmentAtBBL = true;
				LOGnOUT(4,<<" !!! Warning !!!: after BBL likelihood did not improve"<< currBestL<<" "<<_bestL<<endl);
			}			
			//_bestL = max(currBestL,_bestL);
			string treeINodes = gainLossOptions::_outDir + "//" + "TheTree.INodes.iter"+ int2string(iter) + ".ph"; 
			printTree (_tr, treeINodes);
		}

// ROOT optimization
		if (gainLossOptions::_performOptimizationsROOT) 
		{
			currBestL = optimizeRoot();
			if (currBestL>_bestL) {
				_bestL = currBestL;
			}
			else{
				LOGnOUT(4,<<" !!! Warning !!!: after Root likelihood did not improve"<< currBestL<<" "<<_bestL<<endl);
			}			
			//_bestL = max(currBestL,_bestL);
		}		
		if ( (_bestL-previousL) < max(_epsilonOptimization, abs(_bestL/10000))  || noLikeImprovmentAtBBL) // stop Opt for less than epsilon likelihood point
		{
			LOGnOUT(3,<<" OverAll optimization converged. Iter= "<<iter<<" Likelihood="<<_bestL<<endl);
			break;
		}
		if(gainLossOptions::_simulatedAnnealing){
			_epsilonOptimization = max(_epsilonOptimization*gainLossOptions::_simulatedAnnealingCoolingFactor,0.3*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
			_epsilonOptimizationModel = max(_epsilonOptimizationModel*gainLossOptions::_simulatedAnnealingCoolingFactor,0.1*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
			_epsilonOptimizationBBL = max(_epsilonOptimizationBBL*gainLossOptions::_simulatedAnnealingCoolingFactor,0.2*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
		}
	}
	if (iter>_maxNumOfIterations)
		LOGnOUT(4,<<" Too many="<<iter-1<<" iterations in Model+BBL. Last optimized parameters are used."<<endl);
	time(&t2);
	LOGnOUT(4,<<"Optimization RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}
/********************************************************************************************
optimizationsSPvv
*********************************************************************************************/
void gainLossOptimizer::optimizationsSPvv(){
	time_t t1;
	time(&t1);
	time_t t2;
	
	LOGnOUT(4,<<"-------------------------------"<<endl
		<<"Starting optimizations: maxNumIterations="<<_maxNumOfIterations<<endl);	
	_bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	MDOUBLE currBestL=VERYSMALL;
	MDOUBLE previousL;
	bool noLikeImprovmentAtBBL = false; // if BBL did not produce a new tree, end. (no point in another iteration model+BBL)

	bool isSkipParamsOptimization = gainLossOptions::_isSkipFirstParamsOptimization;
	int iter;
	LOGnOUT(3,<<endl<<"#########"<<" optimization starting epsilonCycle="<<_epsilonOptimization<<" maxNumIterations="<<_maxNumOfIterations<<endl);
	LOGnOUT(3,<<"start optimization with L= "<<_bestL<<endl);
	for (iter=1;iter<=_maxNumOfIterations;iter++){
		LOGnOUT(4,<<endl<<"------"<<" Model+BBL iter="<<iter<<endl);
		previousL = _bestL;		// breaking out of loop when no (>epsilon) improvement is made by comparing to previousL
// model optimization
		if(!isSkipParamsOptimization){
			currBestL = optimizeParametersSPvv();}
		else{
			LOGnOUT(4,<<"Optimize Params - Skipped"<<endl);
		}		
		if (currBestL>_bestL) {
			_bestL = currBestL;
		}
		else if(!isSkipParamsOptimization && currBestL<_bestL){
			LOGnOUT(4,<<" !!! Warning !!!: after model optimization likelihood went down"<< currBestL<<" "<<_bestL<<endl);
		}		
		//_bestL = max(currBestL,_bestL);
		isSkipParamsOptimization = false;	// only first iteration skipped
// ROOT optimization
		if (gainLossOptions::_performOptimizationsROOT) 
		{
			currBestL = optimizeRootSPvv();
			if (currBestL>_bestL) {
				_bestL = currBestL;
			}
			else{
				LOGnOUT(4,<<" !!! Warning !!!: after Root likelihood did not improve"<< currBestL<<" "<<_bestL<<endl);
			}		
			//_bestL = max(currBestL,_bestL);
		}
// BBL optimization
		if (gainLossOptions::_performOptimizationsBBL && _performOptimizationsBBL){
			//LOGnOUT(4,<<"Start BBL... with epsilonOptimizationBBL= "<<_epsilonOptimizationBBL<<endl);
			if(gainLossOptions::_performOptimizationsBBLOnlyOnce)	// Next iteration - no BBL
				_performOptimizationsBBL = false;
			currBestL = optimizeBranchLengthsvv(iter);
			if (currBestL>_bestL) {
				_bestL = currBestL;
			}
			else{
				noLikeImprovmentAtBBL = true;
				LOGnOUT(4,<<" !!! Warning !!!: after BBL likelihood did not improve"<< currBestL<<" "<<_bestL<<endl);
			}		
			//_bestL = max(currBestL,_bestL);
			string treeINodes = gainLossOptions::_outDir + "//" + "TheTree.INodes.iter"+ int2string(iter) + ".ph"; 
			printTree (_tr, treeINodes);
		}
		if ((_bestL-previousL) < max(_epsilonOptimization, abs(_bestL/10000)) || noLikeImprovmentAtBBL ) // stop Opt for less than 2 likelihood point
		{
			LOGnOUT(3,<<" OverAll optimization converged. Iter= "<<iter<<" Likelihood="<<_bestL<<endl);
			break;
		}
		if(gainLossOptions::_simulatedAnnealing){
			_epsilonOptimization = max(_epsilonOptimization*gainLossOptions::_simulatedAnnealingCoolingFactor,0.3*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
			_epsilonOptimizationModel = max(_epsilonOptimizationModel*gainLossOptions::_simulatedAnnealingCoolingFactor,0.1*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
			_epsilonOptimizationBBL = max(_epsilonOptimizationBBL*gainLossOptions::_simulatedAnnealingCoolingFactor,0.2*gainLossOptions::_simulatedAnnealingMinEpsilonFactor); // simulated annealing
		}
	}
	if (iter>_maxNumOfIterations)
		LOGnOUT(4,<<" Too many="<<iter-1<<" iterations in Model+BBL. Last optimized parameters are used."<<endl);
	time(&t2);
	LOGnOUT(4,<<"Optimization RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}


/********************************************************************************************
optimizeParameters
*********************************************************************************************/
MDOUBLE  gainLossOptimizer::optimizeParameters(){
	//LOGnOUT(4,<<"Starting optimizeParameters with: numOfIterations="<<_maxNumOfIterations<<" and epsilonOptimization="<<_epsilonOptimizationModel<<endl);
	time_t t1;
	time(&t1);
	time_t t2;

	optimizeGainLossModel* opt = new optimizeGainLossModel(_tr,*_sp,_sc,
		_isReversible,_epsilonOptimizationModel,_maxNumOfIterationsModel,_weightsUniqPatterns,_unObservableData_p);
	//optimizeGainLossModel* opt = new optimizeGainLossModel(_tr,*_sp,_sc,
	//	_isReversible,_epsilonOptimizationModel,_maxNumOfIterationsModel,_plogLforMissingData);


	LOGnOUT(4,<<"-------------------------------"<<endl<<"Model optimization over with: "<<endl);	
	_bestGain=opt->getBestMu1();
	LOGnOUT(4,<<"Gain "<<_bestGain<<endl);

	//if (!gainLossOptions::_isReversible) {
		//MDOUBLE bestM2=opt->getBestMu2();
		_bestLoss = static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->getMu2();
		LOGnOUT(4,<<"Loss "<<_bestLoss<<endl);
		LOGnOUT(4,<<"	Gain/Loss  ratio= "<< _bestGain/_bestLoss<<endl);

	//}	
	if(isAlphaOptimization((*_sp).distr())){
		_bestAlphaRate=opt->getBestAlpha();
		LOGnOUT(4,<<"AlphaRateRate "<<_bestAlphaRate<<endl);
	}
	if(isBetaOptimization((*_sp).distr())){
		_bestBetaRate=opt->getBestBeta();
		LOGnOUT(4,<<"BetaRate "<<_bestBetaRate<<endl);
		LOGnOUT(4,<<" Rate Expectancy = "<< _bestAlphaRate/_bestBetaRate<<endl);
		LOGnOUT(4,<<" Rate Standatd Deviation = "<< sqrt(_bestAlphaRate/(_bestAlphaRate*_bestBetaRate))<<endl);
	}
	if(isInvariantOptimization((*_sp).distr())){
		_bestRateProbInvariant =opt->getBestRateProbInvariant();
		LOGnOUT(4,<<"ProbInvariantRate "<<_bestRateProbInvariant<<endl);
	}
	if(dynamic_cast<mixtureDistribution*>((*_sp).distr())){
		printMixtureParams();	
	}
	if (isThetaOptimization() && !gainLossOptions::_isRootFreqEQstationary) {
		_bestTheta=opt->getBestTheta();
		LOGnOUT(4,<<"Theta "<<_bestTheta<<endl);
		_freq[1] = _bestTheta;
		_freq[0] = 1-_freq[1];
	}
	else{
		_freq[1] = _bestGain/(_bestGain+_bestLoss);
		_freq[0] = 1-_freq[1];
	}
	MDOUBLE bestL = opt->getBestL();
	
	MDOUBLE currentlogL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
	if(!DEQUAL(currentlogL,bestL)){ //
		LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-bestL <<"\n");
	}
	LOGnOUT(4,<<"updated likelihood (after optimizeParameters)= "<<bestL<<endl);
	time(&t2);
	LOGnOUT(4,<<"Model optimization RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	if(opt) delete opt;
	return bestL;	
}
/********************************************************************************************
optimizeParametersSPvv
*********************************************************************************************/
MDOUBLE  gainLossOptimizer::optimizeParametersSPvv(){
	//LOGnOUT(4,<<"Starting optimizeParametersSPvv with: numOfIterations="<<_maxNumOfIterationsModel<<" and epsilonOptimization="<<_epsilonOptimizationModel<<endl);
	time_t t1;
	time(&t1);
	time_t t2;	

	optimizeGainLossModelVV* opt = new optimizeGainLossModelVV(_tr,_spVVec,_sc,_gainDist,_lossDist,
		gainLossOptions::_isReversible,_epsilonOptimizationModel,_maxNumOfIterationsModel,_weightsUniqPatterns,_unObservableData_p);
	LOGnOUT(4,<<"-------------------------------"<<endl<<"Model optimization over with: "<<endl);

	_bestGainAlpha=opt->getBestGainAlpha();
	LOGnOUT(4,<<"AlphaGain "<<_bestGainAlpha<<endl);
	if(isBetaOptimization(_gainDist)){
		_bestGainBeta=opt->getBestGainBeta();
		LOGnOUT(4,<<"BetaGain "<<_bestGainBeta<<endl);
		_gainExp = _bestGainAlpha/_bestGainBeta;
		_gainSTD = sqrt(_bestGainAlpha/(_bestGainBeta*_bestGainBeta));
		LOGnOUT(4,<<" Gain Expectation = "<< rateExpectation(_gainDist)<<endl);
		LOGnOUT(4,<<" Gain Expectancy = "<< _gainExp<<endl);
		LOGnOUT(4,<<" Gain Standard Deviation = "<<_gainSTD<<endl);
	}
	if(isInvariantOptimization(_gainDist)){
		_bestGainProbInvariant=opt->getBestGainProbInvariant();
		LOGnOUT(4,<<"ProbInvariantGain "<<_bestGainProbInvariant<<endl);
	}
	if (!gainLossOptions::_isReversible) {
		_bestLossAlpha=opt->getBestLossAlpha();
		LOGnOUT(4,<<"AlphaLoss "<<_bestLossAlpha<<endl);
		if(isBetaOptimization(_lossDist)){
			_bestLossBeta=opt->getBestLossBeta();
			_lossExp = _bestLossAlpha/_bestLossBeta;
			_lossSTD = sqrt(_bestLossAlpha/(_bestLossBeta*_bestLossBeta));
			LOGnOUT(4,<<"BetaLoss "<<_bestLossBeta<<endl);
			LOGnOUT(4,<<" Loss Expectation = "<< rateExpectation(_lossDist)<<endl);
			LOGnOUT(4,<<" Loss Expectancy = "<< _lossExp<<endl);
			LOGnOUT(4,<<" Loss Standard Deviation = "<<_lossSTD<<endl);
			LOGnOUT(4,<<"	Gain/Loss Expectancy ratio= "<< (_bestGainAlpha/_bestGainBeta)/(_bestLossAlpha/_bestLossBeta)<<endl);
			LOGnOUT(4,<<"	Expectancy(Gain)/Expectancy(Loss)  by computation = "<< computeExpOfGainByExpOfLossRatio(_gainDist, _lossDist)<<endl);
		}
		if(isInvariantOptimization(_lossDist)){
			_bestLossProbInvariant=opt->getBestLossProbInvariant();
			LOGnOUT(4,<<"ProbInvariantLoss "<<_bestLossProbInvariant<<endl);
		}
	}
	if(isAlphaOptimization((*_spVVec[0][0]).distr())){
		MDOUBLE bestAlpha=opt->getBestRateAlpha();
		LOGnOUT(4,<<"AlphaRate "<<bestAlpha<<endl);
	}
	if(isInvariantOptimization((*_spVVec[0][0]).distr())){
		MDOUBLE bestRateProbInvariant =opt->getBestRateProbInvariant();
		LOGnOUT(4,<<"ProbInvariantRate "<<bestRateProbInvariant<<endl);
	}
	if (isThetaOptimization() && !gainLossOptions::_isRootFreqEQstationary) {
		_bestTheta=opt->getBestTheta();
		LOGnOUT(4,<<"Theta "<<_bestTheta<<endl);
		_freq[1] = _bestTheta;
		_freq[0] = 1-_freq[1];
	}
	else{
		_freq[1] = _gainExp/(_gainExp+_lossExp);
		_freq[0] = 1-_freq[1];
	}
	MDOUBLE bestL = opt->getBestL();

	//if(_unObservableData_p) _unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
	MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	if(!DEQUAL(currentlogL,bestL)){ //DEQUAL(currentlogL,bestL)
		LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-bestL <<"\n");
	}
	LOGnOUT(4,<<"updated likelihood (after optimizeParameters)= "<<bestL<<endl);
	time(&t2);
	LOGnOUT(4,<<"Model optimization RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	if(opt) delete opt;
	// update spSimple based on latest est parameters - deleted at the end
	if(_spSimple) 	delete _spSimple;
	_spSimple = startStochasticProcessSimpleGamma(_gainExp, _lossExp, _freq);
	return bestL;	
}

/********************************************************************************************
optimizeRoot
*********************************************************************************************/
MDOUBLE gainLossOptimizer::optimizeRootSPvv()
{
	time_t t1;
	time(&t1);
	time_t t2;	

	MDOUBLE oldL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	MDOUBLE newL = VERYSMALL;
	tree tempTree = _tr;
	LOGnOUT(4,<<"*** Starting optimizeRoot="<<oldL<<endl);

	treeIterDownTopConst tIt(tempTree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isLeaf())
			continue;

		tempTree.rootAt(mynode);
		if(_unObservableData_p)	_unObservableData_p->setLforMissingData(tempTree,_spVVec,_gainDist,_lossDist);
		MDOUBLE newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tempTree,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
		if((newL>oldL+_epsilonOptimizationBBL*gainLossOptions::_epsilonForReRootFactor)){ // only for substantial improvement the tree will be re-rooted
			_tr = tempTree;
			LOGnOUT(4,<<"tree rooted at "<<_tr.getRoot()->name()<<" with lL="<<newL<<endl);
			LOGnOUT(4,<<"sons of root are "<<_tr.getRoot()->getSon(0)->name()<<" , "<<_tr.getRoot()->getSon(1)->name()<<" , "<<_tr.getRoot()->getSon(2)->name()<<endl);
		}
		else{
			if(_unObservableData_p)	_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist); 
		}

	}
	LOGnOUT(4,<<"*** After optimizeRoot="<<max(newL,oldL)<<endl);
	time(&t2);
	LOGnOUT(4,<<"optimizeRoot RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return newL;
}
/********************************************************************************************
optimizeRoot
*********************************************************************************************/
MDOUBLE gainLossOptimizer::optimizeRoot()
{
	time_t t1;
	time(&t1);
	time_t t2;

	MDOUBLE oldL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
	MDOUBLE newL = VERYSMALL;
	tree tempTree = _tr;
	LOGnOUT(4,<<"*** Starting optimizeRoot= "<<oldL<<endl);

	treeIterDownTopConst tIt(tempTree);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isLeaf())
			continue;
		tempTree.rootAt(mynode);

		if(_unObservableData_p){
			_unObservableData_p->setLforMissingData(tempTree,_sp);
		}
		MDOUBLE newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tempTree,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
		if((newL>oldL+_epsilonOptimizationBBL*gainLossOptions::_epsilonForReRootFactor))// only for substantial improvement the tree will be re-rooted
		{
			_tr = tempTree;
			LOGnOUT(4,<<"tree re-rooted at "<<_tr.getRoot()->name()<<" with lL="<<newL<<endl);
			LOGnOUT(4,<<"sons of root are "<<_tr.getRoot()->getSon(0)->name()<<" , "<<_tr.getRoot()->getSon(1)->name()<<" , "<<_tr.getRoot()->getSon(2)->name()<<endl);
		}
		else{
			if(_unObservableData_p)	_unObservableData_p->setLforMissingData(_tr,_sp);	// go back to original tree value
		}

	}
	LOGnOUT(4,<<"*** After optimizeRoot=    "<<max(newL,oldL)<<endl);
	time(&t2);
	LOGnOUT(4,<<"optimizeRoot RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return newL;
}


/********************************************************************************************
optimizeBranchLengths
*********************************************************************************************/
MDOUBLE gainLossOptimizer::optimizeBranchLengths(const int outerIter)
{
	time_t t1;
	time(&t1);
	time_t t2;
	MDOUBLE tollForPairwiseDist=0.01; // the BBL default, epsilon per branch (brent's value)
	MDOUBLE bblEMfactor = 10;

	int numberOfBranchs = _tr.getNodesNum();
	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/5; // (is 1.5) for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationBBLIter = _epsilonOptimizationBBL*epsilonOptimizationIterFactor/bblEMfactor;	// The next iteration epsilon, multiply per-branch value

	MDOUBLE oldL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
	MDOUBLE newLnoUnObservableDataCorrection = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	tree oldTree = _tr;
	bool isFixedRoot = !_isReversible && !gainLossOptions::_isRootFreqEQstationary;

	MDOUBLE minLikeImprovmentForNoLS = 5.0;
	MDOUBLE minLikeImprovmentForNoSkip = 2.0;

	if(gainLossOptions::_isBblLS){
		if (gainLossOptions::_isBblEMbeforeLSWithMissSpecifiedModel) {
			// start with BBL-EM and additional iteration of Line-Search optimizeBranches
			MDOUBLE oldLnoUnObservableDataCorrection = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,NULL);
			LOGnOUT(4,<<"*** Prior BBL-EM... followed by addtional iteration of Line-Search "<<"\t"<<oldL<<endl);			
			if(isFixedRoot){
				LOGnOUT(4,<<"*** Start Fix Root BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				bblEMfixRoot bblEM1(_tr, _sc, *_sp, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_tr,_sp);
				newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			}
			else{
				LOGnOUT(4,<<"*** Start BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				// Note: likelihood does not improve with iterations compared to the likelihood under correction for UnObs, hence NULL
				bblEM bblEM1(_tr, _sc, *_sp, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_tr,_sp);
				newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			}
			
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			newL = bbl.optimizeBranches(_tr,_sp,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL, 1 ,bestL);
			if(newL<oldL){
				_tr = oldTree;
				LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
			}
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}
		// only Line-Search optimizeBranches
		else{
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			newL = bbl.optimizeBranches(_tr,_sp,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL,_maxNumOfIterationsBBL,bestL);
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}
	}
	else{
		if(!_isSkipBblEM){
			MDOUBLE oldLnoUnObservableDataCorrection = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,NULL);
			if(isFixedRoot){
				LOGnOUT(4,<<"*** Start Fix Root BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				bblEMfixRoot bblEM1(_tr, _sc, *_sp, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_tr,_sp);
				newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			}
			else{
				LOGnOUT(4,<<"*** Start BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				// Note: likelihood does not improve with iterations compared to the likelihood under correction for UnObs, hence NULL
				bblEM bblEM1(_tr, _sc, *_sp, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_tr,_sp);
				newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			}
		}
		// if include LS when BBL-EM fails
		if((newL-oldL < max(minLikeImprovmentForNoLS, abs(newL/10000)) ) && _isbblLSWhenbblEMdontImprove){ // Do LS if less than 5 likelihood points were gained
			LOGnOUT(4,<<" Only "<< newL-oldL<<" improvement with BBL-EM -> Perform BBL-LS one iteration"<<endl);
			if(gainLossOptions::_isSkipBblEMWhenbblEMdontImprove && (newL-oldL < minLikeImprovmentForNoSkip)){
				LOGnOUT(4,<<"Since no improvement (less than "<<minLikeImprovmentForNoSkip<<"), BBL-EM will be skipped next iteration, go directly to LS "<<endl);
				_isSkipBblEM = true; // once BBL-EM is not improving Like, next time - skip
			}
			_tr = oldTree;
			if(_unObservableData_p)
				_unObservableData_p->setLforMissingData(_tr,_sp);
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
			newL = bbl.optimizeBranches(_tr,_sp,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL, 1 ,bestL);
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}	
	}

	LOGnOUT(4,<<"*** After Branch Lengths Opt. returned Likelihood="<<"\t"<<newL<<endl);
	if((newL<oldL)){
		_tr = oldTree;
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_tr,_sp);
		newL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
		LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
	}	
	
	MDOUBLE postL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weightsUniqPatterns,_unObservableData_p);
	if( !DEQUAL(newL,postL) ){
		LOGnOUT(3,<<"***ERROR***: Diff returned L, and re-calculated L"<<" "<<newL<<" "<<postL<<" "<<postL-newL<<endl);
	}
	LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
	time(&t2);
	LOGnOUT(4,<<"Branch Lengths RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return newL;
}



/********************************************************************************************
optimizeBranchLengthsvv - new
*********************************************************************************************/
MDOUBLE gainLossOptimizer::optimizeBranchLengthsvv(const int outerIter)
{
	time_t t1;
	time(&t1);
	time_t t2;
	MDOUBLE tollForPairwiseDist=0.01; // the BBL default, epsilon per branch (brent's value)
	MDOUBLE bblEMfactor = 10;

	int numberOfBranchs = _tr.getNodesNum();
	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/5; // (is 1.5) for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationBBLIter = _epsilonOptimizationBBL*epsilonOptimizationIterFactor/bblEMfactor;	// The next iteration epsilon, multiply per-branch value
	MDOUBLE oldL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	
	MDOUBLE newLnoUnObservableDataCorrection = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	tree oldTree = _tr;
	//bool isFixedRoot = !_isReversible && !gainLossOptions::_isRootFreqEQstationary;

	MDOUBLE minLikeImprovmentForNoLS = 5.0;
	MDOUBLE minLikeImprovmentForNoSkip = 2.0;
	
	if(gainLossOptions::_isBblLS){
		if (gainLossOptions::_isBblEMbeforeLSWithMissSpecifiedModel) {
			// start with BBL-EM and additional iteration of Line-Search optimizeBranches
			MDOUBLE oldLnoUnObservableDataCorrection = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,NULL);
			LOGnOUT(4,<<"*** Prior BBL-EM... followed by addtional iteration of Line-Search "<<"\t"<<oldL<<endl);			
			//if(isFixedRoot){
			//	LOGnOUT(4,<<"*** spSimple - Fix Root BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
			//	LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
			//	bblEMfixRoot bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
			//		,NULL,&oldLnoUnObservableDataCorrection);
			//	newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
			//	if(_unObservableData_p)	
			//		_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
			//	newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
			//}
			//else{
				LOGnOUT(4,<<"*** _spSimple BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				// Note: likelihood does not improve with iterations compared to the likelihood under correction for UnObs, hence NULL
				bblEM bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)	
					_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
				newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
			//}
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);	// can be sent
			newL = bbl.optimizeBranches(_tr,_spVVec,_gainDist,_lossDist,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL,_maxNumOfIterationsBBL,bestL);
			if(newL<oldL){
				_tr = oldTree;
				LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
			}
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}
		// only Line-Search optimizeBranches
		else{
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);	// can be sent
			newL = bbl.optimizeBranches(_tr,_spVVec,_gainDist,_lossDist,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL,_maxNumOfIterationsBBL,bestL);
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}
	}
	else{
		if(!_isSkipBblEM){
			MDOUBLE oldLnoUnObservableDataCorrection = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,NULL);
			//if(isFixedRoot){
			//	LOGnOUT(4,<<"*** _spSimple Fix Root BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
			//	LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
			//	bblEMfixRoot bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
			//		,NULL,&oldLnoUnObservableDataCorrection);
			//	newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
			//	if(_unObservableData_p)	
			//		_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
			//	newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
			//}
			//else{
				LOGnOUT(4,<<"*** _spSimple BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
				LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
				// Note: likelihood does not improve with iterations compared to the likelihood under correction for UnObs, hence NULL
				bblEM bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist
					,NULL,NULL); // optional - likelihood &oldLnoUnObservableDataCorrection
				LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
				if(_unObservableData_p)	
					_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
				newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
			//}
		}
		// if include LS when BBL-EM fails

		if((newL-oldL < max(minLikeImprovmentForNoLS, abs(newL/10000)) ) && _isbblLSWhenbblEMdontImprove){ // Do LS if less than 5 likelihood points were gained
			LOGnOUT(4,<<" Only "<< newL-oldL<<" improvement with BBL-EM -> Perform BBL-LS one iteration"<<endl);
			if(gainLossOptions::_isSkipBblEMWhenbblEMdontImprove && (newL-oldL < minLikeImprovmentForNoSkip)){
				LOGnOUT(4,<<"Since no improvement (less than "<<minLikeImprovmentForNoSkip<<"), BBL-EM will be skipped next iteration, go directly to LS "<<endl);
				_isSkipBblEM = true; // once BBL-EM is not improving Like, next time - skip
			}
			_tr = oldTree;
			if(_unObservableData_p)	
				_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
			bblLS bbl;
			MDOUBLE bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);	// can be sent
			newL = bbl.optimizeBranches(_tr,_spVVec,_gainDist,_lossDist,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL,1,bestL);
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
		}	
	}
	LOGnOUT(4,<<"*** After Branch Lengths Opt. returned Likelihood="<<"\t"<<newL<<endl);
	if((newL<oldL)){
		_tr = oldTree;
		if(_unObservableData_p)	
			_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
		newL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
		LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
	}
	MDOUBLE postL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	if( !DEQUAL(newL,postL) ){
		LOGnOUT(3,<<"***ERROR***: Diff returned L, and re-calculated L"<<" "<<newL<<" "<<postL<<" "<<postL-newL<<endl);
	}
	LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
	time(&t2);
	LOGnOUT(4,<<"Branch Lengths RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return newL;
}



/********************************************************************************************
optimizeBranchLengthsvv
*********************************************************************************************/
//MDOUBLE gainLossOptimizer::optimizeBranchLengthsvv(const int outerIter)
//{
//	time_t t1;
//	time(&t1);
//	time_t t2;
//	MDOUBLE tollForPairwiseDist=0.01; // the BBL default, epsilon per branch (brent's value)
//	MDOUBLE bblEMfactor = 10;
//
//	int numberOfBranchs = _tr.getNodesNum();
//	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/5; // (is 1.5) for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
//	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
//	MDOUBLE epsilonOptimizationBBLIter = _epsilonOptimizationBBL*epsilonOptimizationIterFactor/bblEMfactor;	// The next iteration epsilon, multiply per-branch value
//	MDOUBLE oldL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
//
//	MDOUBLE newL = VERYSMALL;
//	tree oldTree = _tr;
//	bool isFixedRoot = !_isReversible && !gainLossOptions::_isRootFreqEQstationary;
//
//	if(gainLossOptions::_isBblLS){
//		if (gainLossOptions::_isBblEMbeforeLSWithMissSpecifiedModel) // not implemented to work well 
//		{
//			LOGnOUT(4,<<"*** Prior BBL-EM under stationarity assumption and one Sp\n WARN: this is not working well..."<<"\t"<<oldL<<endl);
//			LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
//			if(isFixedRoot){
//				LOGnOUT(4,<<"*** Start Fix Root BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
//				bblEMfixRoot bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor) , epsilonOptimizationBBLIter,tollForPairwiseDist,_unObservableData_p,&oldL);
//				newL = bblEM1.getTreeLikelihood();
//			}
//			else{
//				LOGnOUT(4,<<"*** Start BBL-EM Optimization with Likelihood="<<"\t"<<oldL<<endl);
//				bblEM bblEM1(_tr, _sc, *_spSimple, NULL, (int)(_maxNumOfIterationsBBL*bblEMfactor), epsilonOptimizationBBLIter,tollForPairwiseDist,_unObservableData_p,&oldL);
//				newL = bblEM1.getTreeLikelihood();
//			}
//			LOGnOUT(4,<<"*** After BBL-EM Likelihood (with wrong model)"<<"\t"<<newL<<endl);
//			if(newL<oldL){
//				_tr = oldTree;
//				LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
//			}
//		}
//		bblLS bbl;
//		MDOUBLE bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);	// can be sent
//		newL = bbl.optimizeBranches(_tr,_spVVec,_gainDist,_lossDist,_sc,_weightsUniqPatterns,_unObservableData_p,outerIter,_epsilonOptimizationBBL,_maxNumOfIterationsBBL,bestL);
//	}
//	else{
//		LOGnOUT(4,<<"!!!WARNING!!!: BBL-EM is not implemented for spVVec and is not performed."<<endl);
//	}
//
//	LOGnOUT(4,<<"*** After Branch Lengths Opt. returned Likelihood="<<"\t"<<newL<<endl);
//	if((newL<oldL)){
//		_tr = oldTree;
//		LOGnOUT(4,<<"NOTE: No improvment-> Retain previous tree"<<endl);
//	}
//	if(_unObservableData_p)	
//		_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
//	MDOUBLE postL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
//
//	if( !DEQUAL(newL,postL) ){
//		LOGnOUT(3,<<"***ERROR***: Diff returned L, and re-calculated L"<<" "<<newL<<" "<<postL<<" "<<postL-newL<<endl);
//	}
//
//	time(&t2);
//	LOGnOUT(4,<<"Branch Lengths RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
//	return newL;
//}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptimizer::printMixtureParams() 
{	
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(_sp->distr());
	for (int k = 0; k < pMixture->getComponentsNum(); ++k)
	{
		LOGnOUT(4, << "comp="<<k<<" Alp/Beta= "<<pMixture->getAlpha(k)/pMixture->getBeta(k)<<" alpha= "<<pMixture->getAlpha(k) << " beta= " <<pMixture->getBeta(k)<<" Prob= "<<pMixture->getComponentProb(k)<<endl);  
	}
}
/********************************************************************************************
*********************************************************************************************/
gainLossOptions::distributionType getRateDistributionType(distribution* dist)
{
	gainLossOptions::distributionType res;
	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)){
		res = gainLossOptions::GENERAL_GAMMA_PLUS_INV;
	}
	else if (dynamic_cast<generalGammaDistributionFixedCategories*>(dist)){
		res = gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES;
	}
	else if (dynamic_cast<generalGammaDistribution*>(dist)){
		res = gainLossOptions::GENERAL_GAMMA;
	}
	else{
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType");
	}
	return res;
}
