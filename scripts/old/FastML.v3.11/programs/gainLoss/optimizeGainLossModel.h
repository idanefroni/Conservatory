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


#ifndef ___OPTIMIZE_GLM
#define ___OPTIMIZE_GLM

#include "bblEM.h"
#include "bestAlpha.h"
#include "computePijComponent.h"
#include "computeUpAlg.h"
#include "definitions.h"
#include "gainLossModel.h"
#include "gammaDistribution.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include "numRec.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "talRandom.h"
#include "gainLossUtils.h"
#include "gainLossOptions.h"
#include "unObservableData.h"
#include "GamMixtureOptimizer.h"
#include "gammaDistributionFixedCategories.h"
#include "generalGammaDistributionPlusInvariant.h"
#include "mixtureDistribution.h"
#include "gammaUtilities.h"
#include "gainLossOptions.h"



class optimizeGainLossModel {
public:
	//explicit optimizeGainLossModel(const tree& tr, stochasticProcess& sp, const sequenceContainer &sc,
	//	const bool isReversible =false, /*const bool evalTheta =true,*/ 
	//	const MDOUBLE epsilonOptimization =0.1, const int numIterations =10,
	//	MDOUBLE* logLforMissingData =NULL,
	//	ostream& out=cout);
	explicit optimizeGainLossModel(const tree& tr, stochasticProcess& sp, const sequenceContainer &sc,
		const bool isReversible =false, /*const bool evalTheta =true,*/ 
		MDOUBLE epsilonOptimization =0.1, const int numIterations =10,
		Vdouble*  weights = NULL,
		unObservableData* unObservableData_p=NULL);


	//bool isUpdateGain(const MDOUBLE currBestL, MDOUBLE& currM1, const MDOUBLE lossLikelihoodImprovmet);
	MDOUBLE getBestMu1() {return _bestMu1;}
	MDOUBLE getBestMu2() {return _bestMu2;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestBeta() {return _bestBeta;}
	MDOUBLE getBestRateProbInvariant() {return _bestRateProbInvariant;}
	MDOUBLE getBestL() {return _bestL;}
	//void initMissingDataInfo();

	//MDOUBLE*   startingBestAlphaFixedTree(tree& tr,sequenceContainer& sc,stochasticProcess& sp);

private:
	MDOUBLE _bestMu1;
	MDOUBLE _bestMu2; // for non-reversible model only
	MDOUBLE _bestGainLossRatio;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestBeta;
	MDOUBLE _bestTheta;
	MDOUBLE _bestRateProbInvariant;
	MDOUBLE _bestL;
	////MDOUBLE _logLforMissingData;
	//MDOUBLE* _plogLforMissingData;
	//Vdouble* _pLforMissingDataPerCat;
	unObservableData* _unObservableData_p;
	Vdouble* _weightsUniqPatterns;
};

/********************************************************************************************
*********************************************************************************************/
/********************************************************************************************
*********************************************************************************************/
class C_evalParam{
public:
	C_evalParam(const tree& tr,
		const stochasticProcess& sp, const sequenceContainer &sc, int which_mu, bool isReversible,Vdouble*  weights, const unObservableData* unObservableData_p)
		: _tr(tr),/*_sp(sp),*/ _sc(sc),_which_param(which_mu),_isReversible(isReversible),_weights(weights)
	{
		
		_sp = sp.clone();	// the original sp is not effected
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;
		//unObservableData currUnObs(*unObservableData_p);
		
		//_weights = gainLossOptions::_weights;
		
		//if(gainLossOptions::_accountForMissingData){ // plogLforMissingData is not sent but it is needed (the change is local)
		//	_plogLforMissingData = &_logLforMissingData;
		//}
		//else{
		//	_plogLforMissingData = NULL;
		//}
		if ((_which_param>6) || (_which_param<0))
			errorMsg::reportError("Error in C_evalParam, error at _which_param");
	};
	virtual ~C_evalParam(){
		if(_sp)	delete _sp;
		if(_unObservableData_p)	delete _unObservableData_p;
	}

private:
	const tree& _tr;
	stochasticProcess* _sp;
	const sequenceContainer &_sc;
	int _which_param;
	bool _isReversible;
	unObservableData* _unObservableData_p;
	Vdouble* _weights;

public:
	enum paramName {gain,loss,rateAlpha,rateBeta,theta,rateProbInvariant,gainLossRatio};

	MDOUBLE operator() (MDOUBLE param) {
		MDOUBLE sumPijQij = 1.0;
		switch (_which_param) {
			case (C_evalParam::gain) : static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->setMu1(param,_isReversible); break;
			case (C_evalParam::loss) : static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->setMu2(param); break;
			case (C_evalParam::rateAlpha) : setRateAlpha(_sp->distr(),param);  break;
			case (C_evalParam::rateBeta) : setRateBeta(_sp->distr(),param);  break;
			case (C_evalParam::theta) : (static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel()))->setTheta(param); break;
			case (C_evalParam::rateProbInvariant) : static_cast<generalGammaDistributionPlusInvariant*>(_sp->distr())->setInvProb(param); break;
			case (C_evalParam::gainLossRatio) :
				if(gainLossOptions::_isOptimizeParamsWithLogMinMax) param = pow(10,param);
				static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->setMu1(sqrt(param),_isReversible);
				static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->setMu2( sqrt(1.0/param) );
				//norm_factor = normalizeQ(_sp);
				break;
		}
		sumPijQij = normalizeQ(_sp);
		if(_unObservableData_p){ 	_unObservableData_p->setLforMissingData(_tr,_sp); 	}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*_sp,_weights,_unObservableData_p);
		(static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel()))->norm( sumPijQij ); // reverse the normalization after likelihood computation.
		LOG(5,<<"for _which_param "<<_which_param<<" with val = "<<param<<" logL = "<<res<<endl);
		return -res;
	}
};



#endif
