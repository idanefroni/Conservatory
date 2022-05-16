// $Id: bestTamura92param.cpp 962 2006-11-07 15:13:34Z privmane $

#include "bestTamura92param.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "bblEMProportionalEB.h"
#include "bblLSProportionalEB.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestTamura92ParamFixedTree::bestTamura92ParamFixedTree(const tree& et, // find best TrTv and theta
													   const sequenceContainer& sc,
													   stochasticProcess& sp,
													   const Vdouble * weights,
													   const int maxTotalIterations,
													   const MDOUBLE epsilonLikelihoodImprovment,
													   const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
													   const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
													   const MDOUBLE upperBoundOnTrTv) {
	LOG(5,<<"Starting bestTamura92ParamFixedTree: find Best TrTv and theta"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = upperBoundOnTrTv*0.3;
	MDOUBLE prevTheta = 0.5;

	for (int i=0; i < maxTotalIterations; ++i) {
		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			oldL = newL;
			_bestL = newL;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
			}
			break;
		}
	}
}

bestTamura92ParamAndBBL::bestTamura92ParamAndBBL(tree& et, //find best TrTv, theta and best BBL
												 const sequenceContainer& sc,
												 stochasticProcess& sp,
												 const Vdouble * weights,
												 const int maxTotalIterations,
												 const MDOUBLE epsilonLikelihoodImprovment,
												 const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
												 const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
												 const MDOUBLE epsilonLoglikelihoodForBBL,
												 const MDOUBLE upperBoundOnTrTv,
												 const int maxBBLIterations){
	LOG(5,<<"Starting bestTamura92ParamAndBBL: find best TrTv, theta and BBL"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = upperBoundOnTrTv*0.3;
	MDOUBLE prevTheta = 0.5;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {
		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTheta(_bestTheta);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			oldL = newL;
			_bestL = newL;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
				et = prevTree;
			}
			break;
		}
	}
}

bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL( //find best TrTv, theta, Alpha and best branch lengths
	tree& et,
	const sequenceContainer& sc,
	stochasticProcess& sp,
	const Vdouble * weights,
	const int maxTotalIterations,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonLoglikelihoodForTrTvOptimization,
	const MDOUBLE epsilonLoglikelihoodForThetaOptimization,
	const MDOUBLE epsilonLoglikelihoodForAlphaOptimization,
	const MDOUBLE epsilonLoglikelihoodForBBL,
	const MDOUBLE upperBoundOnTrTv,
	const int maxBBLIterations,
	const MDOUBLE initAlpha,
	const MDOUBLE upperBoundOnAlpha)

{
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevTrTv = static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel())->getTrTv();
	MDOUBLE prevTheta = static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel())->getTheta();
	MDOUBLE prevAlpha = initAlpha;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {

		// optimize TrTv
		newL = -brent(0.0, prevTrTv, upperBoundOnTrTv,
					  C_evalTrTvParam(et,sc,sp,weights),
					  epsilonLoglikelihoodForTrTvOptimization,
					  &_bestTrTv);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTv);

		// optimize Theta
		newL = -brent(0.0, prevTheta, 1.0,
					  C_evalTheta(et,sc,sp,weights),
					  epsilonLoglikelihoodForThetaOptimization,
					  &_bestTheta);
		(static_cast<tamura92*>(sp.getPijAccelerator()->getReplacementModel()))->changeTheta(_bestTheta);

		// optimize Alpha
		newL = -brent(0.0, prevAlpha, upperBoundOnAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonLoglikelihoodForAlphaOptimization,
					  &_bestAlpha);
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha);
 
		LOG(5,<<"# bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL iteration " << i << ": after param optimization:" <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Alpha = " << _bestAlpha << endl);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();

		LOG(5,<<"# bestTamura92ParamAlphaAndBBL::bestTamura92ParamAlphaAndBBL iteration " << i << ": after branch lengths optimization:" <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
		    oldL = newL;
			_bestL = newL;
			prevTrTv = _bestTrTv;
			prevTheta = _bestTheta;
			prevAlpha = _bestAlpha;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestTrTv = prevTrTv;
				_bestTheta = prevTheta;
				et = prevTree;
			}
		    break;
		}
	}
}

bestTamura92ParamAlphaAndBBLProportional::bestTamura92ParamAlphaAndBBLProportional( //find best TrTv, theta, loca Alpha for each gene, global Alpha and best branch lengths
	tree& et,
	vector<sequenceContainer>& sc,
	multipleStochasticProcess* msp,
 	gammaDistribution* pProportionDist,
	Vdouble initLocalAlphas,
	Vdouble initLocalKappas,
	Vdouble initLocalThetas,
	const MDOUBLE upperBoundOnLocalAlpha,
	const MDOUBLE initGlobalAlpha,
	const MDOUBLE upperBoundOnGlobalAlpha,
	const MDOUBLE upperBoundOnTrTv,
	const int maxTotalIterations,
	const int maxBBLIterations,
	const bool optimizeSelectedBranches,
	const bool optimizeTree,
	const string branchLengthOptimizationMethod,
	const bool optimizeLocalParams,
	const bool optimizeGlobalAlpha,
	const Vdouble * weights,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonLoglikelihoodForLocalTrTvOptimization,
	const MDOUBLE epsilonLoglikelihoodForLocalThetaOptimization,
	const MDOUBLE epsilonLoglikelihoodForLocalAlphaOptimization,
	const MDOUBLE epsilonLoglikelihoodForGlobalAlphaOptimization,
	const MDOUBLE epsilonLoglikelihoodForBBL)

{
	LOG(5,<<"Starting bestTamura92ParamAlphaAndBBLProportional"<<endl);
	Vdouble currentTrTvVec,currentThetaVec,currentLocalAlphaVec;
	MDOUBLE currentGlobalAlpha = initGlobalAlpha;
	currentTrTvVec = initLocalKappas;
	currentThetaVec = initLocalThetas;
	currentLocalAlphaVec = initLocalAlphas;

	Vdouble newLvec;
	newLvec.resize(msp->getSPVecSize());
	//doubleRep oldL(VERYSMALL);//DR
	//doubleRep newL;//DR
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL;
	//doubleRep epsilonLoglikelihoodForGlobalAlphaOptimizationDR(epsilonLoglikelihoodForGlobalAlphaOptimization);//DR
	_bestLvec.resize(msp->getSPVecSize(),0.0);
	_bestLocalAlphaVec = initLocalAlphas;
	_bestGlobalAlpha = initGlobalAlpha;
	int spIndex;
	_bestTrTvVec = currentTrTvVec;
	_bestThetaVec = currentThetaVec;
	pProportionDist->setAlpha(_bestGlobalAlpha);
	for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
		(static_cast<tamura92*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTheta(_bestThetaVec[spIndex]);//safety
		(static_cast<tamura92*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTvVec[spIndex]);
		(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
	}	
	//first compute the likelihood;
	_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);	

	MDOUBLE ax_local = 0.0;
	MDOUBLE c_TrTv_x = upperBoundOnTrTv;
	MDOUBLE c_theta_x = 1.0;
	MDOUBLE c_localAlpha_x = upperBoundOnLocalAlpha;
	for (int i=0; i < maxTotalIterations; ++i) {
		if(optimizeLocalParams){
			for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
				//optimize Theta
				MDOUBLE theta_x(_bestThetaVec[spIndex]);
				newLvec[spIndex] = -brent(ax_local,theta_x,c_theta_x,
					  C_evalLocalTheta(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonLoglikelihoodForLocalThetaOptimization,
					  &currentThetaVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestThetaVec[spIndex] = currentThetaVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing TrTv param"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}				
				(static_cast<tamura92*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTheta(_bestThetaVec[spIndex]);//safety
				
				//optimize TrTv
				MDOUBLE TrTv_x(_bestTrTvVec[spIndex]);
				newLvec[spIndex] = -brent(ax_local,TrTv_x,c_TrTv_x,
					  C_evalLocalTrTvParam(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonLoglikelihoodForLocalTrTvOptimization,
					  &currentTrTvVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestTrTvVec[spIndex] = currentTrTvVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing TrTv param"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<tamura92*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestTrTvVec[spIndex]);//safety

				//optimize local alpha
				MDOUBLE localAlpha_x(_bestLocalAlphaVec[spIndex]);
				newLvec[spIndex] = -brent(ax_local,localAlpha_x, c_localAlpha_x,
					  C_evalLocalAlpha(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonLoglikelihoodForLocalAlphaOptimization,
					  &currentLocalAlphaVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestLocalAlphaVec[spIndex] = currentLocalAlphaVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing local alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]); //safety
			}
			LOGnOUT(2,<<"Done with Tamura92 local params optimization. LL: "<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Local Params:"<<endl);
			LOGnOUT(2,<<"TrTv:");
			for(spIndex = 0;spIndex < _bestTrTvVec.size();++spIndex){
				LOGnOUT(2,<<_bestTrTvVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"Theta:");
			for(spIndex = 0;spIndex < _bestThetaVec.size();++spIndex){
				LOGnOUT(2,<<_bestThetaVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"local alpha:");
			for(spIndex = 0;spIndex < _bestLocalAlphaVec.size();++spIndex){
				LOGnOUT(2,<<_bestLocalAlphaVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
		}
		if(optimizeGlobalAlpha){
			//doubleRep ax_global(0.0);//DR
			//doubleRep c_globalAlpha_x(upperBoundOnGlobalAlpha);//DR
			//doubleRep minusOne(-1.0);//DR
			MDOUBLE ax_global = 0.0;
			MDOUBLE c_globalAlpha_x = upperBoundOnGlobalAlpha;

			//optimize global alpha
			//doubleRep globalAlpha_x(prevGlobalAlpha);//DR
			MDOUBLE globalAlpha_x = _bestGlobalAlpha;
			//newL = minusOne*brentDoubleRep(ax_global,globalAlpha_x,c_globalAlpha_x,
			//		C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
			//		epsilonLoglikelihoodForGlobalAlphaOptimizationDR,
			//		&_bestGlobalAlpha);//DR
			newL = -brent(ax_global,globalAlpha_x,c_globalAlpha_x,
					C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
					epsilonLoglikelihoodForGlobalAlphaOptimization,
					&currentGlobalAlpha);
			if (newL >= sumVdouble(_bestLvec))
			{
				_bestGlobalAlpha = currentGlobalAlpha;
			} 
			else
			{//likelihood went down!
				LOG(2,<<"likelihood went down in optimizing global alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
			}
			pProportionDist->setAlpha(_bestGlobalAlpha); //safety
			//whether or not likelihood has improved we need to update _bestLvec 
			_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
			LOGnOUT(2,<<"Done with global alpha optimization"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Global Alpha:"<<_bestGlobalAlpha<<endl);
		}

		if(optimizeTree)
		{
			if(branchLengthOptimizationMethod == "bblLS"){
				bblLSProportionalEB bblLSPEB1(et,sc,msp,pProportionDist,_bestLvec,optimizeSelectedBranches,maxBBLIterations,epsilonLoglikelihoodForBBL);
				_bestLvec = bblLSPEB1.getTreeLikelihoodVec();
				LOGnOUT(2,<<"Done with bblLS"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			else if(branchLengthOptimizationMethod == "bblEM"){
				bblEMProportionalEB bblEMPEB1(et,sc,msp,pProportionDist,optimizeSelectedBranches,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);
				_bestLvec = bblEMPEB1.getTreeLikelihood();
				LOGnOUT(2,<<"Done with bblEM. LL: "<<sumVdouble(_bestLvec)<<endl);
			}
			LOGnOUT(2,<<et.stringTreeInPhylipTreeFormat()<<endl);
		}
		// check for improvement in the likelihood
		if (sumVdouble(_bestLvec) > oldL+epsilonLikelihoodImprovment) {
			//all params have already been updated
			oldL = sumVdouble(_bestLvec);
		} else {
			break;
		}
		LOGnOUT(4,<<"Done with optimization iteration "<<i<<". LL: "<<sumVdouble(_bestLvec)<<endl);
	}
}

