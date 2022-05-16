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


#ifndef ___OPTIMIZE_GLM_VV
#define ___OPTIMIZE_GLM_VV

#include "bblEM.h"
#include "bestAlpha.h"
#include "computePijComponent.h"
#include "computeUpAlg.h"
#include "definitions.h"
#include "gainLossModel.h"
#include "gammaDistribution.h"
#include "generalGammaDistribution.h"
#include "generalGammaDistributionPlusInvariant.h"
#include "distributionPlusInvariant.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include "numRec.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "talRandom.h"
#include "gainLossUtils.h"
#include "unObservableData.h"


class optimizeGainLossModelVV {
public:
	explicit optimizeGainLossModelVV(const tree& tr, 
		vector<vector<stochasticProcess*> >& spVVec, const sequenceContainer &sc,
		distribution * gainDist, distribution * lossDist,
		const bool isReversible,  
		MDOUBLE epsilonOptimization, const int numIterations,
		Vdouble*  weights,
		unObservableData* unObservableData_p);

	MDOUBLE getBestGainAlpha() {return _bestGainAlpha;}
	MDOUBLE getBestGainBeta() {return _bestGainBeta;}
	MDOUBLE getBestGainProbInvariant() {return _bestGainProbInvariant;}

	MDOUBLE getBestLossAlpha() {return _bestLossAlpha;}
	MDOUBLE getBestLossBeta() {return _bestLossBeta;}
	MDOUBLE getBestLossProbInvariant() {return _bestLossProbInvariant;}

	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestRateAlpha() {return _bestRateAlpha;}
	MDOUBLE getBestRateProbInvariant() {return _bestRateProbInvariant;}

	MDOUBLE getBestL() {return _bestL;}

private:
	MDOUBLE _bestGainAlpha;
	MDOUBLE _bestGainBeta;
	MDOUBLE _bestGainProbInvariant;

	MDOUBLE _bestLossAlpha; // for non-reversible model only
	MDOUBLE _bestLossBeta; 
	MDOUBLE _bestLossProbInvariant;

	MDOUBLE _bestRateAlpha;
	MDOUBLE _bestRateProbInvariant;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
	MDOUBLE _bestGainLossRatio;
	unObservableData* _unObservableData_p;
	Vdouble* _weightsUniqPatterns;
};

/********************************************************************************************
*********************************************************************************************/
/********************************************************************************************
*********************************************************************************************/
class C_evalParamVV {
public:
	C_evalParamVV(const tree& tr,
		const vector<vector<stochasticProcess*> >& spVVec, const sequenceContainer &sc, int which_mu, 
		const distribution* gainDist, const distribution* lossDist,
		bool isReversible,Vdouble*  weights , const unObservableData* unObservableData_p)
		: _tr(tr),_sc(sc),_which_param(which_mu),_isReversible(isReversible),_weights(weights)
	{
		_gainDist=gainDist->clone();
		_lossDist=lossDist->clone();	
		_spVVec.resize(_gainDist->categories());
		for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
			_spVVec[gainCategor].resize(_lossDist->categories());
			for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){		
				_spVVec[gainCategor][lossCategor] = spVVec[gainCategor][lossCategor]->clone();
			}
		}
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;

	};
	virtual ~C_evalParamVV(){
		if(_spVVec[0][0]){
			for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
				for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){		
					delete _spVVec[gainCategor][lossCategor];
				}
			}
		}
		if(_gainDist)
			delete _gainDist;
		if(_lossDist)
			delete _lossDist;
		if(_unObservableData_p)
			delete _unObservableData_p;

	}
private:
	const tree& _tr;
	vector<vector<stochasticProcess*> > _spVVec;
	distribution* _gainDist; 
	distribution* _lossDist;
	const sequenceContainer &_sc;
	int _which_param;
	bool _isReversible;
	unObservableData* _unObservableData_p;
	Vdouble* _weights;

public:
	enum paramName {gainAlpha,gainBeta,gainProbInvariant,lossAlpha,lossBeta,lossProbInvariant,rateAlpha,rateProbInvariant,theta,gainLossRatio};

	MDOUBLE operator() (MDOUBLE param) {
		MDOUBLE gainLossRatioToCompleteByBeta = 1;
		MDOUBLE sumPijQij = 1;
		MDOUBLE previousAlpha = 1;
		MDOUBLE increaseToGainLossRatioInducedByAlphaModification = 1;

		switch (_which_param) {
			case (C_evalParamVV::gainAlpha) :
				if(1){ // keep gainLossRatio
					previousAlpha = getRateAlpha(_gainDist);
					increaseToGainLossRatioInducedByAlphaModification = param/previousAlpha;
					updateGainBeta(getRateBeta(_gainDist) * increaseToGainLossRatioInducedByAlphaModification,_spVVec,_gainDist,_lossDist);
				}
				updateGainAlpha(param,_spVVec,_gainDist,_lossDist);
				break;
			case (C_evalParamVV::gainBeta) : updateGainBeta(param,_spVVec,_gainDist,_lossDist); break;
			case (C_evalParamVV::gainProbInvariant) : updateGainProbInvariant(param,_gainDist); break;

			case (C_evalParamVV::lossAlpha) :
				if(1){ // keep gainLossRatio
					previousAlpha = getRateAlpha(_lossDist);
					increaseToGainLossRatioInducedByAlphaModification = param/previousAlpha;
					updateLossBeta(getRateBeta(_lossDist) * increaseToGainLossRatioInducedByAlphaModification,_spVVec,_gainDist,_lossDist);
				}
				updateLossAlpha(param,_spVVec,_gainDist,_lossDist);
				break;
			case (C_evalParamVV::lossBeta) : updateLossBeta(param,_spVVec,_gainDist,_lossDist); break;
			case (C_evalParamVV::lossProbInvariant) : updateLossProbInvariant(param,_lossDist); break;

			case (C_evalParamVV::gainLossRatio) :
				if(gainLossOptions::_isOptimizeParamsWithLogMinMax) param = pow(10,param);
				gainLossRatioToCompleteByBeta = param * (getRateAlpha(_lossDist)/getRateAlpha(_gainDist));
				if(gainLossOptions::_isUpdateOnlyGainBetaForRatio)
					updateGainBeta(getRateBeta(_lossDist)/gainLossRatioToCompleteByBeta,_spVVec,_gainDist,_lossDist);
				else{
					updateGainBeta(sqrt(1.0/gainLossRatioToCompleteByBeta),_spVVec,_gainDist,_lossDist);
					updateLossBeta(sqrt(gainLossRatioToCompleteByBeta),_spVVec,_gainDist,_lossDist);
				}				
				//norm_factor = normalizeQ(_spVVec, _gainDist, _lossDist);
				break;
			case (C_evalParamVV::rateAlpha) : updateRateAlpha(param,_spVVec,_gainDist,_lossDist); break;
			case (C_evalParamVV::rateProbInvariant) : updateRateProbInvariant(param,_spVVec,_gainDist,_lossDist); break;
			case (C_evalParamVV::theta) : updateTheta(param,_spVVec,_gainDist,_lossDist); break;	
		}
		sumPijQij = normalizeQ(_spVVec, _gainDist, _lossDist);

		if(_unObservableData_p)	_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
		MDOUBLE res = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,_spVVec,_gainDist,_lossDist,_weights,_unObservableData_p);
		normVec(sumPijQij,_spVVec, _gainDist, _lossDist); // reverse the normalization after likelihood computation.
		LOG(5,<<"with val= "<<param<<" which_param:: "<<_which_param<<" L="<<res<<endl);
		return -res;
	}
};


#endif
