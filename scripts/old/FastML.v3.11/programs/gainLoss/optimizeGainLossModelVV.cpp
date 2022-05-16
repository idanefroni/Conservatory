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
#include "optimizeGainLossModelVV.h"
#include "gainLossUtils.h"
#include "gainLossOptions.h"
#include "Parameters.h"

/********************************************************************************************
optimizeGainLossModel - for gain,Loss ~ Gamma(Alpha,Beta)
*********************************************************************************************/
optimizeGainLossModelVV::optimizeGainLossModelVV(const tree& tr, 
											 vector<vector<stochasticProcess*> >& spVVec, const sequenceContainer &sc,
											 distribution * gainDist, distribution * lossDist,
											 const bool isReversible, 
											 MDOUBLE epsilonOptimization, const int numIterations,
											 Vdouble*  weights,
											 unObservableData* unObservableData_p):
_weightsUniqPatterns(weights),_unObservableData_p(unObservableData_p)
{
	MDOUBLE MINIMUM_ALPHA_PARAM;
	if(gainLossOptions::_isAlphaLimit){
		MINIMUM_ALPHA_PARAM = 0.1;
	}
	else{
		MINIMUM_ALPHA_PARAM = ::MINIMUM_ALPHA_PARAM;
	}
	bool isAllowHigherAlpha = false;	// for distribution more 'gaussian' and Eq, need higher alpha, else 10.0
	MDOUBLE MAXIMUM_ALPHA_PARAM;
	if(isAllowHigherAlpha){
		MAXIMUM_ALPHA_PARAM = 100;
	}
	else{
		MAXIMUM_ALPHA_PARAM = ::MAXIMUM_ALPHA_PARAM;
	}
	MDOUBLE MINMUM_GAIN_LOSS_RATIO_PARAM;
	MDOUBLE MAXIMUM_GAIN_LOSS_RATIO_PARAM;
	if(gainLossOptions::_isOptimizeParamsWithLogMinMax){
		MINMUM_GAIN_LOSS_RATIO_PARAM = log10(::MINMUM_GAIN_LOSS_RATIO_PARAM);
		MAXIMUM_GAIN_LOSS_RATIO_PARAM = log10(::MAXIMUM_GAIN_LOSS_RATIO_PARAM);
	}else{
		MINMUM_GAIN_LOSS_RATIO_PARAM = ::MINMUM_GAIN_LOSS_RATIO_PARAM;
		MAXIMUM_GAIN_LOSS_RATIO_PARAM = ::MAXIMUM_GAIN_LOSS_RATIO_PARAM;
	}


	stochasticProcess sp = *spVVec[0][0];	

	bool optimizeBetaGain = isBetaOptimization(gainDist);
	bool optimizeBetaLoss = isBetaOptimization(lossDist);
	bool optimizeAlphasGainLoss = true;
	if(gainLossOptions::_optimizationLevel<=2){ // Vlow and below
		optimizeAlphasGainLoss = false;
		LOGnOUT(4,<<"No optimization of rate shape (Alphas) in low optimization level"<<endl);
	}
	bool optimizeGLProbInvariant = isInvariantOptimization(gainDist); // for both gain and loss
	
	bool optimizeRateAlpha = isAlphaOptimization((sp.distr()));
	bool optimizeRateProbInvariant = isInvariantOptimization((sp.distr())); 
	bool evalTheta = isThetaOptimization();
	
	MDOUBLE currBestL,currGainAlpha,currGainBeta,currGainProbInvariant,currLossAlpha,currLossBeta,currLossProbInvariant,currRateAlpha,currRateProbInvariant,currTheta,previousL,currGainLossRatio;
	MDOUBLE sumPijQij;
	//distribution* gainDistPrev=gainDist->clone();
	//distribution* lossDistPrev=lossDist->clone();
	//vector<vector<stochasticProcess*> > spVVecPrev;
	//spVVecPrev.resize(_gainDist->categories());
	//for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
	//	_spVVec[gainCategor].resize(_lossDist->categories());
	//	for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){		
	//		spVVecPrev[gainCategor][lossCategor] = spVVec[gainCategor][lossCategor]->clone();
	//	}
	//}
	//unObservableData* unObservableData_pPrev;
	//if(unObservableData_p)
	//	unObservableData_pPrev = unObservableData_p->clone();
	//else
	//	unObservableData_pPrev = NULL;

//Random Starts
	//unObservableData* currUnObservableData_p;
	//if(gainLossOptions::_accountForMissingData){
	//	currUnObservableData_p = new unObservableData(sc, &sp, gainLossAlphabet(),gainLossOptions::_minNumOfOnes);
	//	currUnObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
	//}
	//else{
	//	currUnObservableData_p = NULL;
	//}
	if(gainLossOptions::_initParamsAtRandPointsInOptimization){
		currGainAlpha =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_ALPHA_PARAM, MAXIMUM_ALPHA_PARAM);
		currGainBeta=talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_BETA_PARAM, MAXIMUM_BETA_PARAM);
		currGainProbInvariant = talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM);
		currLossAlpha =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_ALPHA_PARAM, MAXIMUM_ALPHA_PARAM);
		currLossBeta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_BETA_PARAM, MAXIMUM_BETA_PARAM); 
		currLossProbInvariant =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM);

        currRateProbInvariant =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM);		
		currRateAlpha =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_ALPHA_PARAM, MAXIMUM_ALPHA_PARAM);
		currTheta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_FREQ_PARAM, MAXIMUM_FREQ_PARAM);
	}
	else{
		currBestL=VERYSMALL;
		currGainAlpha=1;	//Gain
		currGainBeta=1;
		currGainProbInvariant = 0.1;
		currLossAlpha=1; // Loss (for non-reversible model only)
		currLossBeta=1; 
		currLossProbInvariant = 0.1;

		currRateAlpha=1;	//Rate
		currRateProbInvariant = 0.1;
		currTheta = 0.5;
		currGainLossRatio = 1;
	}

	int numberOfParameters = 1;
// initialize
	// Gain
	_bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
	if(optimizeGLProbInvariant) {
		_bestGainProbInvariant = static_cast<generalGammaDistributionPlusInvariant*>(gainDist)->getInvProb();
		++numberOfParameters;
	}
	//_bestGainAlpha = static_cast<generalGammaDistribution*>(gainDist)->getAlpha();
	_bestGainAlpha = getRateAlpha(gainDist);

	//if(optimizeBetaGain) _bestGainBeta = static_cast<generalGammaDistribution*>(gainDist)->getBeta();
	if(optimizeBetaGain) {
		_bestGainBeta = getRateBeta(gainDist);
		++numberOfParameters;
	}	
	// Loss
	if (!isReversible){
		if(optimizeGLProbInvariant) {
			_bestLossProbInvariant = static_cast<generalGammaDistributionPlusInvariant*>(lossDist)->getInvProb();
			++numberOfParameters;
		}
		//_bestLossAlpha = static_cast<generalGammaDistribution*>(lossDist)->getAlpha();
		//if(optimizeBetaLoss) _bestLossBeta = static_cast<generalGammaDistribution*>(lossDist)->getBeta();
		_bestLossAlpha = getRateAlpha(lossDist);
		if(optimizeBetaLoss){
			_bestLossBeta = getRateBeta(lossDist);
			++numberOfParameters;
		}
	}
	// overall rate
	if(optimizeRateAlpha){
		_bestRateAlpha = getRateAlpha(static_cast<gammaDistribution*>(sp.distr()));
		++numberOfParameters;
	}
	if(optimizeRateProbInvariant){
		_bestRateProbInvariant = static_cast<generalGammaDistributionPlusInvariant*>((sp.distr()))->getInvProb();
		++numberOfParameters;	}

	if(evalTheta){
		++numberOfParameters;
	}
	_bestTheta = static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->getTheta(); // taken either way
	_bestGainLossRatio = computeExpOfGainByExpOfLossRatio(gainDist, lossDist); //(_bestGainAlpha/_bestGainBeta)/(_bestLossAlpha/_bestLossBeta);
	MDOUBLE epsilonOptimizationIterFactor = numberOfParameters;  
	epsilonOptimizationIterFactor = max(3.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationIter = epsilonOptimization*epsilonOptimizationIterFactor;	// for e=0.1 next iteration only for ~0.5 logL points

// optimize
	LOGnOUT(3,<<"### "<<"optimization starting- epsilonOptParam="<<epsilonOptimization<<" epsilonOptIter= "<<epsilonOptimizationIter<<", MaxNumIterations="<<numIterations<<endl);
	LOGnOUT(3,<<"start optimization with:" <<endl<<" L= "<<_bestL<<endl
		<<" gainLossRatio= "<<_bestGainLossRatio<<endl
		<<" GainAlpha= "<<_bestGainAlpha<<" GainBeta= "<<_bestGainBeta<<endl
		<<" LossAlpha= "<<_bestLossAlpha<<" LossBeta= "<<_bestLossBeta<<endl);
	if(evalTheta)	LOGnOUT(3,<<" Theta= "<<_bestTheta);	
	if(optimizeGLProbInvariant) LOGnOUT(3,<<" GainProbInvariant= "<<_bestGainProbInvariant<<" LossProbInvariant= "<<_bestLossProbInvariant<<endl);
	if(optimizeRateAlpha) LOGnOUT(3,<<" RateAlpha= "<<_bestRateAlpha<<endl);
	if(optimizeRateProbInvariant) LOGnOUT(3,<<" RateProbInvariant= "<<_bestRateProbInvariant<<endl);
	
	int iter;
	for (iter=1;iter<=numIterations;iter++)
	{
		previousL = _bestL;	// before loop likelihood
		LOGnOUT(4,<<"\n---- iter="<<iter<<endl);		
		//bool isOptimizeModelParametersInRandomOrderNoReturns = true;
		//int numOfParameters = 9;
		//Vint paramsAlreadyOptimizedV;
		//int curParam;
		//for(int parInd=1; parInd<=numOfParameters; ++parInd){
		//	bool isParamAlreadyOptimized = true;
		//	if(isOptimizeModelParametersInRandomOrderNoReturns){
		//		while (isParamAlreadyOptimized) {
		//			curParam = floor(talRandom::giveRandomNumberBetweenTwoPoints(1,numOfParameters+1));
		//			int::iterator begin = paramsAlreadyOptimizedV.begin();
		//			int::iterator end = paramsAlreadyOptimizedV.end();


		//			if(! paramsAlreadyOptimizedV.find(begin,end,curParam)){
		//				isParamAlreadyOptimized = false;
		//			}
		//		}

		//	}
		//	else if (gainLossOptions::_isStartWithTheta) {
		//		curParam = C_evalParamVV::theta;

		//	}
		//	else{
		//		curParam = parInd;
		//	}
		//}

// optimization - Freq (Theta)
		if (gainLossOptions::_isStartWithTheta && evalTheta && !gainLossOptions::_isRootFreqEQstationary){
			currBestL = -brent(MINIMUM_FREQ_PARAM,_bestTheta,MAXIMUM_FREQ_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::theta, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),
					epsilonOptimization*gainLossOptions::_epsilonOptimizationThetaFactor,&currTheta);
			if (currBestL>_bestL) {
				updateTheta(currTheta,spVVec,gainDist,lossDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST				
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tTheta="<<currTheta<<endl);
				_bestTheta=currTheta;
				_bestL=currBestL;
			}
		}
		// gainLoss ratio
		if(gainLossOptions::_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately && !Parameters::getInt("_keepUserGainLossRatio")){
			currBestL = -brent(MINMUM_GAIN_LOSS_RATIO_PARAM,_bestGainLossRatio,MAXIMUM_GAIN_LOSS_RATIO_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::gainLossRatio, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currGainLossRatio);		
			if(gainLossOptions::_isOptimizeParamsWithLogMinMax) currGainLossRatio = pow(10,currGainLossRatio);
			if (currBestL>_bestL) {
				MDOUBLE gainLossRatioToCompleteByBeta = currGainLossRatio * (getRateAlpha(lossDist)/getRateAlpha(gainDist));
				if(gainLossOptions::_isUpdateOnlyGainBetaForRatio){
					currGainBeta = (getRateBeta(lossDist)/gainLossRatioToCompleteByBeta);
					updateGainBeta(currGainBeta,spVVec,gainDist,lossDist);
				}else{
					currGainBeta = sqrt(1.0/gainLossRatioToCompleteByBeta);
					currLossBeta = sqrt(gainLossRatioToCompleteByBeta);
					updateGainBeta(currGainBeta,spVVec,gainDist,lossDist);
					updateLossBeta(currLossBeta,spVVec,gainDist,lossDist);
				}
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGainLossRatio="<<currGainLossRatio<<endl);
				_bestGainLossRatio=currGainLossRatio;
				_bestGainBeta=currGainBeta;
				_bestLossBeta=currLossBeta;
				_bestL=currBestL;
				MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}

			}
		}else{
			// optimization - GainBeta
			if(optimizeBetaGain && !Parameters::getInt("_keepUserGainLossRatio")){
				currBestL = -brent(MINIMUM_BETA_PARAM,_bestGainBeta,MAXIMUM_BETA_PARAM,
					C_evalParamVV(tr,spVVec,sc,C_evalParamVV::gainBeta, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currGainBeta);		
				if (currBestL>_bestL) {
					updateGainBeta(currGainBeta,spVVec,gainDist,lossDist);
					sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
					if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
					LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGainBeta="<<currGainBeta<<endl);
					_bestGainBeta=currGainBeta;
					_bestL=currBestL;
					//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
					//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
					//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
					//}
				}
			}
	// optimization - LossBeta
			if(optimizeBetaLoss && !Parameters::getInt("_keepUserGainLossRatio")){
				currBestL = -brent(MINIMUM_BETA_PARAM,_bestLossBeta,MAXIMUM_BETA_PARAM,
					C_evalParamVV(tr,spVVec,sc,C_evalParamVV::lossBeta, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currLossBeta);				
				if (currBestL>_bestL) {
					updateLossBeta(currLossBeta,spVVec,gainDist,lossDist);
					sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
					if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
					LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLossBeta="<<currLossBeta<<endl);
					_bestLossBeta=currLossBeta;
					_bestL=currBestL;
					//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
					//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
					//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
					//}
				}
			}
		}

// optimization - GainAlpha
		if(optimizeAlphasGainLoss){
			currBestL = -brent(MINIMUM_ALPHA_PARAM,_bestGainAlpha,MAXIMUM_ALPHA_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::gainAlpha, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currGainAlpha);		
			if (currBestL>_bestL) {
				if(1){ // keep gainLossRatio
					MDOUBLE previousAlpha = getRateAlpha(gainDist);
					MDOUBLE increaseToGainLossRatioInducedByAlphaModification = currGainAlpha/previousAlpha;
					currGainBeta = getRateBeta(gainDist)*increaseToGainLossRatioInducedByAlphaModification;
					updateGainBeta( currGainBeta, spVVec,gainDist,lossDist);
					_bestGainBeta = currGainBeta;
				}
				updateGainAlpha(currGainAlpha,spVVec,gainDist,lossDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGainAlpha="<<currGainAlpha<<endl);
				_bestGainAlpha=currGainAlpha;
				_bestL=currBestL;
							//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
			//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
			//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
			//}
			}
		}
// optimization - GainProbInvariant
		if(optimizeGLProbInvariant && !Parameters::getInt("_keepUserGainLossRatio")){
			currBestL = -brent(MINIMUM_PROB_PARAM,_bestGainProbInvariant,MAXIMUM_PROB_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::gainProbInvariant, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currGainProbInvariant);		
			if (currBestL>_bestL) {
				updateGainProbInvariant(currGainProbInvariant,gainDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGainProbInvariant="<<currGainProbInvariant<<endl);
				_bestGainProbInvariant=currGainProbInvariant;
				_bestL=currBestL;
				//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}
			}
		}
// optimization - LossAlpha
		if (!isReversible ){
			if(optimizeAlphasGainLoss){
				currBestL = -brent(MINIMUM_ALPHA_PARAM,_bestLossAlpha,MAXIMUM_ALPHA_PARAM,
					C_evalParamVV(tr,spVVec,sc,C_evalParamVV::lossAlpha, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currLossAlpha);
				if (currBestL>_bestL) {
					if(1){ // keep gainLossRatio
						MDOUBLE previousAlpha = getRateAlpha(lossDist);
						MDOUBLE increaseToGainLossRatioInducedByAlphaModification = currLossAlpha/previousAlpha;
						currLossBeta = getRateBeta(lossDist)*increaseToGainLossRatioInducedByAlphaModification;
						updateLossBeta(  currLossBeta, spVVec,gainDist,lossDist);
						_bestLossBeta = currLossBeta;				
					}
					updateLossAlpha(currLossAlpha,spVVec,gainDist,lossDist);
					sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
					if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
					LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLossAlpha="<<currLossAlpha<<endl);
					_bestLossAlpha=currLossAlpha;
					_bestL=currBestL;
					//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
					//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
					//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
					//}
				}
			}
// optimization - LossprobInvariant
			if(optimizeGLProbInvariant && !Parameters::getInt("_keepUserGainLossRatio"))
			{
				currBestL = -brent(MINIMUM_PROB_PARAM,_bestLossProbInvariant,MAXIMUM_PROB_PARAM,
					C_evalParamVV(tr,spVVec,sc,C_evalParamVV::lossProbInvariant, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currLossProbInvariant);		
				if (currBestL>_bestL) {
					updateLossProbInvariant(currLossProbInvariant,lossDist);
					sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
					if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
					LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLossProbInvariant="<<currLossProbInvariant<<endl);
					_bestLossProbInvariant=currLossProbInvariant;
					_bestL=currBestL;
					//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
					//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
					//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
					//}
				}
			}
		}
// optimize rateAlpha - additionally (inner sp)...
		if(optimizeRateAlpha){
			currBestL = -brent(MINIMUM_ALPHA_PARAM,_bestRateAlpha,MAXIMUM_ALPHA_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::rateAlpha, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currRateAlpha);
			if (currBestL>_bestL) {
				updateRateAlpha(currRateAlpha,spVVec,gainDist,lossDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tRate Alpha="<<currRateAlpha<<endl);
				_bestRateAlpha=currRateAlpha;
				_bestL=currBestL;
				//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}
			}
		}
// optimization - RateprobInvariant
		if(optimizeRateProbInvariant){
			currBestL = -brent(MINIMUM_PROB_PARAM,_bestRateProbInvariant,MAXIMUM_PROB_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::rateProbInvariant, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currRateProbInvariant);		
			if (currBestL>_bestL) {
				updateRateProbInvariant(currRateProbInvariant,spVVec,gainDist,lossDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLossProbInvariant="<<currRateProbInvariant<<endl);
				_bestRateProbInvariant=currRateProbInvariant;
				_bestL=currBestL;
				//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}
			}
		}		

// optimization - Freq (Theta)
		if (!gainLossOptions::_isStartWithTheta && evalTheta && !gainLossOptions::_isRootFreqEQstationary){
			currBestL = -brent(MINIMUM_FREQ_PARAM,_bestTheta,MAXIMUM_FREQ_PARAM,
				C_evalParamVV(tr,spVVec,sc,C_evalParamVV::theta, gainDist,lossDist,isReversible,_weightsUniqPatterns,_unObservableData_p),
					epsilonOptimization*gainLossOptions::_epsilonOptimizationThetaFactor,&currTheta);
			if (currBestL>_bestL) {
				updateTheta(currTheta,spVVec,gainDist,lossDist);
				sumPijQij = normalizeQ(spVVec, gainDist, lossDist);	// TEST	
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tTheta="<<currTheta<<endl);
				_bestTheta=currTheta;
				_bestL=currBestL;
				//MDOUBLE currentlogL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}
			}
		}
		if (!(_bestL>previousL+epsilonOptimizationIter)){	// previousL is before loop likelihood - if no epsilon improvment => break
			if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist); //not clear needed...
			LOGnOUT(3,<<" model optimization converged. Iter= "<<iter<<" Likelihood="<<_bestL<<endl);
			_bestL=max(_bestL,currBestL);	// not to reduce likelihood. currBestL, returning from brent may be lower
			//if(!DEQUAL(currBestL,_bestL)){ //DEQUAL(currentlogL,bestL)
			//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currBestL-_bestL <<"\n");
			//}
			break;
		}
		if(gainLossOptions::_simulatedAnnealing)
			epsilonOptimization = max(epsilonOptimization*gainLossOptions::_simulatedAnnealingCoolingFactor,0.1*gainLossOptions::_simulatedAnnealingMinEpsilonFactor);  //simulated annealing
	}
	if (iter>=numIterations){
		_bestL=max(_bestL,currBestL);	// not to reduce likelihood. currBestL, returning from brent may be lower
		LOGnOUT(3,<<" Too many iterations in optimizeGainLossModelVV. Iter= "<<iter<< " Last optimized parameters are used."<<endl);
	}
	//if(currUnObservableData_p) delete currUnObservableData_p;
}




