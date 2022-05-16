// $Id: bestAlpha.cpp 10046 2011-12-09 15:35:00Z rubi $

#include <iostream>
using namespace std;

#include "bestAlpha.h"
#include "bblEM.h"
#include "bblEMProportionalEB.h"
#include "bblLSProportionalEB.h"
#include "numRec.h"
#include "logFile.h"
#include "errorMsg.h"

#ifndef VERBOS
#define VERBOS
#endif
//void bestAlpha::checkAllocation() {
//	if (_pi->stocProcessFromLabel(0)->getPijAccelerator() == NULL) {
//		errorMsg::reportError(" error in function findBestAlpha");
//	}
//}
//
// @@@@ The method works with oldL,oldA,bestA and newL,newA.
// Only when it's about to end, the members _bestAlpha and _bestL are filled.

bestAlphaAndBBL::bestAlphaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE initAlpha,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBBL,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
//	LOG(5,<<"find Best Alpha and best BBL"<<endl);
//	LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	const MDOUBLE bx=initAlpha;
	const MDOUBLE ax=0;
	const MDOUBLE cx=upperBoundOnAlpha;
//
	MDOUBLE bestA=0;
	MDOUBLE oldA=0;
	int i=0;
	for (i=0; i < maxTotalIterations; ++i) {
		newL = -brent(ax,bx,cx,
		C_evalAlpha(et,sc,sp,weights),
		epsilonLoglikelihoodForAlphaOptimization,
		&bestA);
 
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBBL::bestAlphaAndBBL iteration " << i <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Alpha = " << bestA << endl);
#endif
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
		    oldL = newL;
		    oldA = bestA;
		} else {
		    oldL = newL;
		    oldA = bestA;

		    
		    _bestL = oldL;
		    _bestAlpha= oldA;
		    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
		    break;
		}

		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBBL::bestAlphaAndBBL iteration " << i <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));
#endif

		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			oldL = newL;
		}
		else {
		    oldL=newL;
		    _bestL = oldL;	
			_bestAlpha= oldA; 
			(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
			break;
		}
	}
	if (i==maxTotalIterations) {
		_bestL = newL;
		_bestAlpha= bestA;
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
	}
}

bestAlphasAndBBLProportional::bestAlphasAndBBLProportional(tree& et, //find Best Alphas (per gene - local and proportional factors - global) and best BBL
					   vector<sequenceContainer>& sc,
					   multipleStochasticProcess* msp,
					   gammaDistribution* pProportionDist,
					   Vdouble initLocalRateAlphas,
					   const MDOUBLE upperBoundOnLocalRateAlpha,
					   const MDOUBLE initGlobalRateAlpha,
					   const MDOUBLE upperBoundOnGlobalRateAlpha,
					   const int maxBBLIterations,
					   const int maxTotalIterations,
					   const bool optimizeSelectedBranches,
					   const bool optimizeTree,
					   const string branchLengthOptimizationMethod,
					   const bool optimizeLocalAlpha,
					   const bool optimizeGlobalAlpha,
					   const Vdouble * weights,
					   const MDOUBLE epsilonLoglikelihoodForLocalRateAlphaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForGlobalRateAlphaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBBL){
//	LOG(5,<<"find Best Alpha and best BBL"<<endl);
//	LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	
	
    if(initLocalRateAlphas.size() != sc.size()){
		MDOUBLE val = initLocalRateAlphas[0]; 
		initLocalRateAlphas.resize(sc.size(),val);
	}
	int spIndex;
	_bestGlobalAlpha = initGlobalRateAlpha;
	pProportionDist->setAlpha(_bestGlobalAlpha);
	_bestLocalAlphaVec = initLocalRateAlphas;
	for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
		(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
	}
	//First compute the likelihood
	_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
	if((!optimizeTree) && (!optimizeLocalAlpha) && (!optimizeGlobalAlpha)) return;
	MDOUBLE currentGlobalAlpha;
	currentGlobalAlpha = initGlobalRateAlpha;
	Vdouble currentLocalAlphaVec;
	Vdouble newLvec;
	//doubleRep newL;//DR
	MDOUBLE newL;
	//doubleRep oldL(VERYSMALL);//DR
	MDOUBLE oldL = VERYSMALL;
	currentLocalAlphaVec = initLocalRateAlphas;
	newLvec.resize(msp->getSPVecSize());
	//doubleRep epsilonLoglikelihoodForGlobalRateAlphaOptimizationDR(epsilonLoglikelihoodForGlobalRateAlphaOptimization);//DR
	string alphas;
	//doubleRep minusOne(-1.0);//DR
	int i;

	MDOUBLE a_localAlpha_x = 0.0;
	MDOUBLE c_localAlpha_x = upperBoundOnLocalRateAlpha;
	for(i=0; i < maxTotalIterations; ++i) {
		//Find best local alphas
		if(optimizeLocalAlpha){
			for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
				MDOUBLE b_localAlpha_x = _bestLocalAlphaVec[spIndex];
				newLvec[spIndex] = -brent(a_localAlpha_x,b_localAlpha_x,c_localAlpha_x,
				C_evalLocalAlpha(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
				epsilonLoglikelihoodForLocalRateAlphaOptimization,
				&currentLocalAlphaVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) {   
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestLocalAlphaVec[spIndex] = currentLocalAlphaVec[spIndex];
				}
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing local alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
			    (static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
			}
			LOGnOUT(2,<<"Done with local alpha optimization"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Local Alphas:");
			for(spIndex = 0;spIndex < _bestLocalAlphaVec.size();++spIndex){
				LOGnOUT(2,<<_bestLocalAlphaVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
		}
		//Find best global alpha
		if(optimizeGlobalAlpha){
			//doubleRep b_globalAlpha_x(_bestGlobalAlpha);//DR
            //doubleRep a_globalAlpha_x(0.0);//DR
			//doubleRep c_globalAlpha_x(upperBoundOnGlobalRateAlpha);//DR
			MDOUBLE b_globalAlpha_x = _bestGlobalAlpha;
            MDOUBLE a_globalAlpha_x = 0.0;
			MDOUBLE c_globalAlpha_x = upperBoundOnGlobalRateAlpha;

			//newL = minusOne*brentDoubleRep(a_globalAlpha_x,b_globalAlpha_x,c_globalAlpha_x,
			//C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
			//epsilonLoglikelihoodForGlobalRateAlphaOptimizationDR,
			//&_bestGlobalAlpha);//DR

			newL = -brent(a_globalAlpha_x,b_globalAlpha_x,c_globalAlpha_x,
			C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
			epsilonLoglikelihoodForGlobalRateAlphaOptimization,
			&currentGlobalAlpha);

			if (newL >= sumVdouble(_bestLvec)) { //converged   
				_bestGlobalAlpha = currentGlobalAlpha;
			}
			else
			{//likelihood went down!
				LOG(2,<<"likelihood went down in optimizing global alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
			}
            pProportionDist->setAlpha(_bestGlobalAlpha);
			//whether or not likelihood has improved we need to update _bestLvec 
			_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
			LOGnOUT(2,<<"Done with global alpha optimization"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Global Alpha:"<<_bestGlobalAlpha<<endl);
		}
		
		if(optimizeTree){
			if(branchLengthOptimizationMethod == "bblLS"){
				bblLSProportionalEB bblLSPEB1(et,sc,msp,pProportionDist,_bestLvec,optimizeSelectedBranches,maxBBLIterations,epsilonLoglikelihoodForBBL);
				_bestLvec = bblLSPEB1.getTreeLikelihoodVec();
				LOGnOUT(2,<<"Done with bblLS"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			else if(branchLengthOptimizationMethod == "bblEM"){
				bblEMProportionalEB bblEMPEB1(et,sc,msp,pProportionDist,optimizeSelectedBranches,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
				_bestLvec = bblEMPEB1.getTreeLikelihood();
				LOGnOUT(2,<<"Done with bblEM"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			LOGnOUT(2,<<et.stringTreeInPhylipTreeFormat()<<endl);
		}
		if (sumVdouble(_bestLvec) > oldL+epsilonLoglikelihoodForBBL) {
			//global and local alpha have already been updated individually
			oldL = sumVdouble(_bestLvec);
		}
		else {
			break;
		}
		LOGnOUT(2,<<"Done with optimization iteration "<<i<<". LL: "<<sumVdouble(_bestLvec)<<endl);
	}
}

bestBetaAndBBL::bestBetaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE initBeta,
					   const MDOUBLE upperBoundOnBeta,
					   const MDOUBLE epsilonLoglikelihoodForBetaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBBL,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
//	LOG(5,<<"find Best Beta and best BBL"<<endl);
//	LOG(5,<<" 1. bestBetaa::findBestBeta"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	const MDOUBLE bx=initBeta;
	const MDOUBLE ax=0;
	const MDOUBLE cx=upperBoundOnBeta;
//
	MDOUBLE bestB=0;
	MDOUBLE oldB=0;
	int i=0;
	for (i=0; i < maxTotalIterations; ++i) {
		newL = -brent(ax,bx,cx,
		C_evalBeta(et,sc,sp,weights),
		epsilonLoglikelihoodForBetaOptimization,
		&bestB);
 
#ifdef VERBOS
		LOG(5,<<"# bestBetaAndBBL::bestBetaAndBBL iteration " << i <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Beta = " << bestB << endl);
#endif
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
		    oldL = newL;
		    oldB = bestB;
		} else {
		    oldL = newL;
		    oldB = bestB;

		    
		    _bestL = oldL;
		    _bestBeta= oldB;
		    (static_cast<gammaDistribution*>(sp.distr()))->setBeta(bestB);
		    break;
		}

		(static_cast<gammaDistribution*>(sp.distr()))->setBeta(bestB);
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();
#ifdef VERBOS
		LOG(5,<<"# bestBetaAndBBL::bestBetaAndBBL iteration " << i <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));
#endif

		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			oldL = newL;
		}
		else {
		    oldL=newL;
		    _bestL = oldL;	
			_bestBeta= oldB; 
			(static_cast<gammaDistribution*>(sp.distr()))->setBeta(bestB);
			break;
		}
	}
	if (i==maxTotalIterations) {
		_bestL = newL;
		_bestBeta= bestB;
		(static_cast<gammaDistribution*>(sp.distr()))->setBeta(bestB);
	}
}

bestAlphaFixedTree::bestAlphaFixedTree(const tree& et, //findBestAlphaFixedTree
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization){
	//LOG(5,<<"findBestAlphaFixedTree"<<endl);
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
	const MDOUBLE bx=static_cast<gammaDistribution*>(sp.distr())->getAlpha();
	const MDOUBLE ax=0.0;

	
	_bestL = -brent(ax,bx,cx,
		C_evalAlpha(et,sc,sp,weights),
		epsilonLoglikelihoodForAlphaOptimization,
		&bestA);
	(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
	_bestAlpha= bestA;
}



bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE initAlpha,
					   const MDOUBLE initBeta,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE upperBoundOnBeta,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBetaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBBL,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
//	LOG(5,<<"find Best Alpha and Beta and best BBL"<<endl);
//	LOG(5,<<" 1. bestAlphaAndBetaAndBBL::findBestAlphaAndBeta"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	MDOUBLE bx=initAlpha;
	const MDOUBLE ax=0;
	const MDOUBLE cx=upperBoundOnAlpha;
	MDOUBLE ex=initBeta;
	const MDOUBLE dx=0;
	const MDOUBLE fx=upperBoundOnBeta;
	bool optimize = false;

//
	MDOUBLE bestA=0;
	MDOUBLE oldA=0;
	MDOUBLE bestB=0;
	MDOUBLE oldB=0;
	int i=0;
	for (i=0; i < maxTotalIterations; ++i) {
//optimize alpha
		newL = -brent(ax,bx,cx,
		C_evalAlpha(et,sc,sp,weights),
		epsilonLoglikelihoodForAlphaOptimization,
		&bestA);
		bx = bestA;
 
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL iteration " << i <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Alpha = " << bestA << endl);
#endif
		if(newL < oldL)
			errorMsg::reportError("likelihood decreased in alhpa optimization step in bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL");
		oldL = newL;
		oldA = bestA;
		_bestL = newL;
		_bestAlpha= bestA;
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			optimize = true;
		}
		(static_cast<generalGammaDistribution*>(sp.distr()))->setAlpha(bestA);

//optimize beta
		newL = -brent(dx,ex,fx,
		C_evalBeta(et,sc,sp,weights),
		epsilonLoglikelihoodForBetaOptimization,
		&bestB);
		ex = bestB;
 
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL iteration " << i <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Beta = " << bestB << endl);
#endif
		if(newL < oldL)
			errorMsg::reportError("likelihood decreased in beta optimization step in bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL");
	    oldL = newL;
	    oldB = bestB;
	    _bestL = oldL;
	    _bestBeta= oldB;
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			optimize = true;
		} 
		(static_cast<generalGammaDistribution*>(sp.distr()))->setBeta(bestB);

//bblEM
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL iteration " << i <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));
#endif
		if(newL < oldL)
			errorMsg::reportError("likelihood decreased in bbl optimization step in bestAlphaAndBetaAndBBL::bestAlphaAndBetaAndBBL");
		oldL = newL;
		_bestL = newL;
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			optimize = true;
		}
		if (!optimize)
			break;
	}
}

