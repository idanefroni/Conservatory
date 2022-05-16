// $Id: bestHKYparam.cpp 10004 2011-11-13 04:40:13Z rubi $

#include "bestHKYparam.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "bblEMProportionalEB.h"
#include "bblLSProportionalEB.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestHkyParamFixedTree::bestHkyParamFixedTree(const tree& et, //findBestHkyParamFixedTree
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnHkyParam,
					   const MDOUBLE epsilonHkyParamOptimization){
	LOG(5,<<"findBestHkyParamFixedTree"<<endl);
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnHkyParam;// left, midle, right limit on HkyParam
	const MDOUBLE bx=cx*0.3;
	const MDOUBLE ax=0;

	
	_bestL = -brent(ax,bx,cx,
		C_evalHkyParam(et,sc,sp,weights),
		epsilonHkyParamOptimization,
		&bestA);
	_bestHkyParam= bestA;
	(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
}

bestHkyParamAndBBL::bestHkyParamAndBBL(tree& et, //find Best HkyParam and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnHkyParam,
					   const MDOUBLE epsilonHkyParamOptimization,
					   const MDOUBLE epsilonLikelihoodImprovment,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
	LOG(5,<<"find Best HkyParam and best BBL"<<endl);
//	LOG(5,<<" 1. bestHkyParam::findBestHkyParam"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	MDOUBLE oldL = VERYSMALL;
	_bestL = VERYSMALL;
	const MDOUBLE bx=upperBoundOnHkyParam*0.3;
	const MDOUBLE ax=0.01;
	const MDOUBLE cx=upperBoundOnHkyParam;
	MDOUBLE bestA=0;
	for (int i=0; i < maxTotalIterations; ++i) {
		_bestL = -brent(ax,bx,cx,
		C_evalHkyParam(et,sc,sp,weights),
		epsilonHkyParamOptimization,
		&bestA);

		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
		} 
		else {//LL converged
			if (_bestL > oldL)
				(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
			else
                _bestL = oldL;
            break;
		}
		_bestHkyParam = bestA;
		(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(bestA);
		LOG(5,<<"bestHkyParamAndBBL: trtv = "<<_bestHkyParam<<endl);
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLikelihoodImprovment);//maxIterations=1000
		_bestL =bblEM1.getTreeLikelihood();
		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
		}
		else {
			_bestL = oldL;
			break;
		}
	}
}

bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL( //find best TrTv (=HkyParam), Alpha and best branch lengths
	tree& et,
	const sequenceContainer& sc,
	stochasticProcess& sp,
	const Vdouble * weights,
	const int maxTotalIterations,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonHkyParamOptimization,
	const MDOUBLE epsilonAlphaOptimization,
	const MDOUBLE epsilonBBL,
	const MDOUBLE upperBoundOnHkyParam,
	const int maxBBLIterations,
	const MDOUBLE initAlpha,
	const MDOUBLE upperBoundOnAlpha)

{
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;

	// first guess for the parameters
	MDOUBLE prevHkyParam = static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel())->getTrTv();
	MDOUBLE prevAlpha = initAlpha;
	tree prevTree;

	for (int i=0; i < maxTotalIterations; ++i) {

		// optimize HkyParam
		newL = -brent(0.0, prevHkyParam, upperBoundOnHkyParam,
					  C_evalHkyParam(et,sc,sp,weights),
					  epsilonHkyParamOptimization,
					  &_bestHkyParam);
		(static_cast<hky*>(sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestHkyParam);
		LOG(5,<<"bestHkyParamAlphaAndBBL: trtv = "<<_bestHkyParam<<endl);
		// optimize Alpha
		newL = -brent(0.0, prevAlpha, upperBoundOnAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonAlphaOptimization,
					  &_bestAlpha);
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha);
		
		LOG(5,<<"# bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL iteration " << i << ": after param optimization:" <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
			  <<"# new hkyParam = " << _bestHkyParam << endl
		      <<"# new Alpha = " << _bestAlpha << endl);

		// optimize branch lengths
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonBBL);
		newL =bblEM1.getTreeLikelihood();

		LOG(5,<<"# bestHkyParamAlphaAndBBL::bestHkyParamAlphaAndBBL iteration " << i << ": after branch lengths optimization:" <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));

		// check for improvement in the likelihood
		if (newL > oldL+epsilonLikelihoodImprovment) {
		    oldL = newL;
			_bestL = newL;
			prevHkyParam = _bestHkyParam;
			prevAlpha = _bestAlpha;
			prevTree = et;
		} else {
			if (newL>oldL) {
				_bestL = newL;
			} else {
				_bestL = oldL;
				_bestHkyParam = prevHkyParam;
				et = prevTree;
			}
		    break;
		}
	}
}

bestHkyParamAlphaAndBBLProportional::bestHkyParamAlphaAndBBLProportional( //find best TrTv (=HkyParam), global Alpha, local Alpha, and best branch lengths
	tree& et,
	vector<sequenceContainer>& sc,
	multipleStochasticProcess* msp,
	gammaDistribution* pProportionDist,
	Vdouble initLocalAlphas,
	Vdouble initLocalKappas,
	const MDOUBLE upperBoundOnLocalAlpha,
	const MDOUBLE initGlobalAlpha,
	const MDOUBLE upperBoundOnGlobalAlpha,
	const MDOUBLE upperBoundOnHkyParam,
	const int maxTotalIterations,
	const int maxBBLIterations,
	const bool optimizeSelectedBranches,
	const bool optimizeTree,
	const string branchLengthOptimizationMethod,
	const bool optimizeLocalParams,
	const bool optimizeGlobalAlpha,
	const Vdouble * weights,
	const MDOUBLE epsilonLikelihoodImprovment,
	const MDOUBLE epsilonHkyParamOptimization,
	const MDOUBLE epsilonLocalAlphaOptimization,
	const MDOUBLE epsilonGlobalAlphaOptimization,
	const MDOUBLE epsilonBBL)

{
	LOG(5,<<"Starting bestHkyParamAlphaAndBBLProportional"<<endl);
	Vdouble current_HkyParamVec,currentLocalAlphaVec;
	MDOUBLE currentGlobalAlpha = initGlobalAlpha;
	current_HkyParamVec = initLocalKappas;
	currentLocalAlphaVec = initLocalAlphas;
	//doubleRep epsilonGlobalAlphaOptimizationDR(epsilonGlobalAlphaOptimization);//DR
	Vdouble newLvec;
	newLvec.resize(msp->getSPVecSize());
	//doubleRep oldL(VERYSMALL);//DR
	//doubleRep newL;
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL;
	_bestLvec.resize(msp->getSPVecSize(),0.0);
	_bestLocalAlphaVec = initLocalAlphas;
	_bestGlobalAlpha = initGlobalAlpha;
	int spIndex;
	//initial HKY params
	_bestHkyParamVec = initLocalKappas;
	pProportionDist->setAlpha(_bestGlobalAlpha);
	for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
		(static_cast<hky*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestHkyParamVec[spIndex]);
		(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
	}
	//first compute the likelihood;
	_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);	

	MDOUBLE ax_local = 0.0;
	MDOUBLE c_HKYParam_x = upperBoundOnHkyParam;
	MDOUBLE c_localAlpha_x = upperBoundOnLocalAlpha;
	for (int i=0; i < maxTotalIterations; ++i) {
		if(optimizeLocalParams){
			for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
				//optimize hky
				MDOUBLE hky_x(_bestHkyParamVec[spIndex]);
				newLvec[spIndex] = -brent(ax_local,hky_x,c_HKYParam_x,
					  C_evalLocalHkyParam(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonHkyParamOptimization,
					  &current_HkyParamVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestHkyParamVec[spIndex] = current_HkyParamVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing hky param"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<hky*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->changeTrTv(_bestHkyParamVec[spIndex]);//safety

				//optimize local alpha
				MDOUBLE localAlpha_x(_bestLocalAlphaVec[spIndex]);
				newLvec[spIndex] = -brent(ax_local,localAlpha_x,c_localAlpha_x,
					  C_evalLocalAlpha(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonLocalAlphaOptimization,
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
				(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
			}
			LOGnOUT(2,<<"Done with HKY local params optimization. LL: "<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Local Params:"<<endl);
			LOGnOUT(2,<<"HHY:");
			for(spIndex = 0;spIndex < _bestHkyParamVec.size();++spIndex){
				LOGnOUT(2,<<_bestHkyParamVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"local alpha:");
			for(spIndex = 0;spIndex < _bestLocalAlphaVec.size();++spIndex){
				LOGnOUT(2,<<_bestLocalAlphaVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
			LOGnOUT(2,<<"LL*: "<<sumVdouble(_bestLvec)<<endl);

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
			//		epsilonGlobalAlphaOptimizationDR,
			//		&_bestGlobalAlpha);//DR
			newL = -brent(ax_global,globalAlpha_x,c_globalAlpha_x,
					C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
					epsilonGlobalAlphaOptimization,
					&currentGlobalAlpha);
			if (newL >= sumVdouble(_bestLvec))
			{
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
		
		if(optimizeTree)
		{
			if(branchLengthOptimizationMethod == "bblLS"){
				bblLSProportionalEB bblLSPEB1(et,sc,msp,pProportionDist,_bestLvec,optimizeSelectedBranches,maxBBLIterations,epsilonBBL);
				_bestLvec = bblLSPEB1.getTreeLikelihoodVec();
				LOGnOUT(2,<<"Done with bblLS"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			else if(branchLengthOptimizationMethod == "bblEM"){
				bblEMProportionalEB bblEMPEB1(et,sc,msp,pProportionDist,optimizeSelectedBranches,NULL,maxBBLIterations,epsilonBBL);
				_bestLvec = bblEMPEB1.getTreeLikelihood();
				LOGnOUT(2,<<"Done with bblEM"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
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

