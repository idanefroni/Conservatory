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
#include "optimizeGainLossModel.h"
#include "Parameters.h"

optimizeGainLossModel::optimizeGainLossModel(const tree& tr, stochasticProcess& sp, const sequenceContainer &sc,
											 const bool isReversible, /*const bool evalTheta,*/ 
											 MDOUBLE epsilonOptimization, const int numIterations,
											 Vdouble*  weights,
											 unObservableData* unObservableData_p):
_weightsUniqPatterns(weights), _unObservableData_p(unObservableData_p)
{
	//_weights = gainLossOptions::_weights;			// since - no weights are used over positions

	MDOUBLE MINIMUM_ALPHA_PARAM;
	if(gainLossOptions::_isAlphaLimit){
		MINIMUM_ALPHA_PARAM = 0.1;
	}
	else{
		MINIMUM_ALPHA_PARAM = ::MINIMUM_ALPHA_PARAM;
	}
	MDOUBLE MINIMUM_GAIN_PARAM;
	if(gainLossOptions::_isGainLimit){
		MINIMUM_GAIN_PARAM = 0.1;
	}
	else{
		MINIMUM_GAIN_PARAM = ::MINIMUM_GAIN_PARAM;
	}
	MDOUBLE MAXIMUM_GAIN_PARAM;
	if(gainLossOptions::_gainLossRateAreFreq){
		MAXIMUM_GAIN_PARAM = 0.9999;
	}
	else{
		MAXIMUM_GAIN_PARAM = ::MAXIMUM_GAIN_PARAM;
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

	bool isAllowHigherAlpha = true;	// for distribution more 'gaussian' and Eq, need higher alpha, else 10.0
	MDOUBLE MAXIMUM_ALPHA_PARAM;
	if(isAllowHigherAlpha){
		MAXIMUM_ALPHA_PARAM = 100;
	}
	else{
		MAXIMUM_ALPHA_PARAM = ::MAXIMUM_ALPHA_PARAM;
	}

	bool optimizeAlpha = isAlphaOptimization(sp.distr());
	bool optimizeBeta = isBetaOptimization(sp.distr());
	bool optimizeMixture = isMixOptimization(sp.distr());
	bool probInvariant = isInvariantOptimization(sp.distr());
	bool evalTheta = isThetaOptimization();

	MDOUBLE previousL;
	MDOUBLE currBestL=VERYSMALL;
	MDOUBLE currM1=0.1;
	MDOUBLE currM2=1; // for non-reversible model only
	MDOUBLE currAlpha=1;
	MDOUBLE currBeta=1;
	MDOUBLE currTheta = 0.5;
	MDOUBLE currRateProbInvariant = 0.05;
	MDOUBLE currGainLossRatio = 1;
	MDOUBLE incrementFactorForGain = gainLossOptions::_slopeFactorForGain;	// forces slow climb for gain param
	MDOUBLE sumPijQij;
// MissingData
	//unObservableData* currUnObservableData_p;
	//if(gainLossOptions::_accountForMissingData){
	//	currUnObservableData_p = new unObservableData(sc, &sp, gainLossAlphabet(),gainLossOptions::_minNumOfOnes);
	//	currUnObservableData_p->setLforMissingData(tr,&sp);
	//}
	//else{
	//	currUnObservableData_p = NULL;
	//}

// currSeeds
	if(gainLossOptions::_initParamsAtRandPointsInOptimization){
		currM1 =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_GAIN_PARAM, MAXIMUM_GAIN_PARAM);
		currM2=talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_LOSS_PARAM, MAXIMUM_LOSS_PARAM);
		currAlpha = talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_ALPHA_PARAM, MAXIMUM_ALPHA_PARAM);
		currBeta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_BETA_PARAM, MAXIMUM_BETA_PARAM);
		currTheta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_FREQ_PARAM, MINIMUM_FREQ_PARAM); 
		currRateProbInvariant =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM);
	}

// initialize - best
	int numberOfParameters = 1;
	_bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,sp,_weightsUniqPatterns,_unObservableData_p); //PerCat
	_bestMu1 = static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->getMu1();
	if (!isReversible){
		_bestMu2 = static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->getMu2();
		++numberOfParameters;
	}
	if(optimizeAlpha){
		_bestAlpha = getRateAlpha(sp.distr());
		++numberOfParameters;
	}
	if(optimizeBeta){
		_bestBeta = getRateBeta(sp.distr());
		++numberOfParameters;
	}
	if(evalTheta)
		++numberOfParameters;
	_bestTheta = static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->getTheta(); // take eiter way	
	if(probInvariant){
		_bestRateProbInvariant = static_cast<generalGammaDistributionPlusInvariant*>(sp.distr())->getInvProb();
		++numberOfParameters;
	}
	_bestGainLossRatio = _bestMu1/_bestMu2;
	MDOUBLE epsilonOptimizationIterFactor = numberOfParameters;  
	epsilonOptimizationIterFactor = max(3.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationIter = epsilonOptimization*epsilonOptimizationIterFactor;	// for e=0.1 next iteration only for ~0.5 logL points

// optimize
	LOGnOUT(3,<<"### "<<"optimization starting- epsilonOptParam="<<epsilonOptimization<<" epsilonOptIter= "<<epsilonOptimizationIter<<", MaxNumIterations="<<numIterations<<endl);
	LOGnOUT(3,<<"start optimization with:"<<endl<<" L= "<<_bestL<<" gainLossRatio= "<<_bestGainLossRatio<<" gain= "<<_bestMu1);
	if(!isReversible) LOGnOUT(3,<<" loss= "<<_bestMu2);
	if(optimizeAlpha) LOGnOUT(3,<<" Alpha= "<<_bestAlpha);
	if(optimizeBeta) LOGnOUT(3,<<" Beta= "<<_bestBeta);
	if(optimizeMixture) LOGnOUT(3,<<" ");
	if(evalTheta)	LOGnOUT(3,<<" Theta= "<<_bestTheta);
	if(probInvariant)	LOGnOUT(3,<<" RateProbInvariant= "<<_bestRateProbInvariant<<"\n");
	if(optimizeMixture)  printMixtureParams(&sp);
	LOGnOUT(3,<<endl);

	int iter;
	for (iter=1;iter<=numIterations;iter++){
		previousL = _bestL;		// breaking out of loop when no (>epsilon) improvement is made by comparing to previousL	
		LOGnOUT(4,<<"\n---- iter="<<iter<<endl);
// optimization - Freq (Theta)
		if (gainLossOptions::_isStartWithTheta && evalTheta && !gainLossOptions::_isRootFreqEQstationary){
			currBestL = -brent(MINIMUM_FREQ_PARAM,_bestTheta,MAXIMUM_FREQ_PARAM,C_evalParam(tr,sp,sc,C_evalParam::theta,isReversible,_weightsUniqPatterns,_unObservableData_p),
				epsilonOptimization*gainLossOptions::_epsilonOptimizationThetaFactor,&currTheta);
			if (currBestL>_bestL) {
				static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setTheta(currTheta);
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tTheta= "<<currTheta<<endl);
				_bestTheta=currTheta;
				_bestL=currBestL;			
			}
		}
// optimize mixture
		if(optimizeMixture){
			GamMixtureOptimizer optGamma(&sp, sc, tr, _unObservableData_p);
			int maxIterations = 1;
			if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::EM) 
				currBestL  = optGamma.findBestParam(GamMixtureOptimizer::EM, maxIterations, epsilonOptimization, gainLossOptions::_weights);
			else if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::ONE_DIM) 
				currBestL  = optGamma.findBestParam(GamMixtureOptimizer::ONE_DIM, maxIterations, epsilonOptimization, gainLossOptions::_weights);
			else errorMsg::reportError("unknown type in gammmaMixtureOptimizerAlgType");
			if (currBestL>_bestL) {
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\timprovment in optimize gammaMixture params"<<endl);
				_bestL=currBestL;
			}
		}
		// gainLoss ratio
		if(gainLossOptions::_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately && !Parameters::getInt("_keepUserGainLossRatio")){
			currBestL = -brent(MINMUM_GAIN_LOSS_RATIO_PARAM , _bestGainLossRatio, MAXIMUM_GAIN_LOSS_RATIO_PARAM, C_evalParam(tr,sp,sc,C_evalParam::gainLossRatio,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currGainLossRatio);
			if(gainLossOptions::_isOptimizeParamsWithLogMinMax) currGainLossRatio = pow(10,currGainLossRatio);
			if (currBestL>_bestL) {
				_bestMu1=sqrt(currGainLossRatio);
				_bestMu2=sqrt(1.0/currGainLossRatio);
				static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setMu1(_bestMu1,isReversible);
				static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->setMu2(_bestMu2);
				sumPijQij = normalizeQ(&sp); //TEST				
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGainLossRatio= "<<currGainLossRatio<<endl);
				_bestGainLossRatio=currGainLossRatio;
				_bestL=currBestL;			
			}
		}else if(!Parameters::getInt("_keepUserGainLossRatio")){
		// optimization - gain
			if(!gainLossOptions::_isSkipGainOptimization){
				if(gainLossOptions::_incrementFactorForGain) 
					currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,min((_bestMu1*incrementFactorForGain),MAXIMUM_GAIN_PARAM),C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currM1);		
				else if(gainLossOptions::_lossBiggerGainLimit) 
					currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,min(_bestMu2,MAXIMUM_GAIN_PARAM),C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currM1);		
				else 
					currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,MAXIMUM_GAIN_PARAM,C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currM1);	
			}
			if (currBestL>_bestL) {
				static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setMu1(currM1,isReversible);
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGain= "<<currM1<<endl);
				_bestMu1=currM1;
				_bestL=currBestL;			
			}
		// optimization - loss
			if (!isReversible & !gainLossOptions::_gainLossRateAreFreq){
				if(gainLossOptions::_lossBiggerGainLimit) currBestL = -brent(max(_bestMu1,MINIMUM_LOSS_PARAM),_bestMu2,MAXIMUM_LOSS_PARAM,C_evalParam(tr,sp,sc,C_evalParam::loss,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currM2);
				else currBestL = -brent(MINIMUM_LOSS_PARAM,_bestMu2,MAXIMUM_LOSS_PARAM,C_evalParam(tr,sp,sc,C_evalParam::loss,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currM2);
				if (currBestL>_bestL) {
					static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->setMu2(currM2);
					sumPijQij = normalizeQ(&sp); //TEST
					if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
					LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLoss= "<<currM2<<endl);
					_bestMu2=currM2;
					_bestL=currBestL;
				}
			}
		}
// optimize Beta 
		if(optimizeBeta && !gainLossOptions::_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately){
			currBestL = -brent(MINIMUM_BETA_PARAM,_bestBeta,MAXIMUM_BETA_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateBeta,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currBeta);			
			if (currBestL>_bestL) {
				setRateBeta(sp.distr(),currBeta);
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tBeta= "<<currBeta<<endl);
				_bestBeta=currBeta;
				_bestL=currBestL;			
			}
		}
// optimize Alpha - 3 options  (all results with same values)
		if(optimizeAlpha){
			currBestL = -brent(MINIMUM_ALPHA_PARAM,_bestAlpha,MAXIMUM_ALPHA_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateAlpha,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currAlpha);			
			if (currBestL>_bestL) {
				setRateAlpha(sp.distr(),currAlpha);
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tAlpha= "<<currAlpha<<endl);
				_bestAlpha=currAlpha;
				_bestL=currBestL;			
			}
		}
// optimization - probInvariant
		if (probInvariant){
			currBestL = -brent(MINIMUM_PROB_PARAM,_bestRateProbInvariant,MAXIMUM_PROB_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateProbInvariant,isReversible,_weightsUniqPatterns,_unObservableData_p),epsilonOptimization,&currRateProbInvariant);
			if (currBestL>_bestL) {
				static_cast<generalGammaDistributionPlusInvariant*>(sp.distr())->setInvProb(currRateProbInvariant);				
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tRateProbInvariant= "<<currRateProbInvariant<<endl);
				_bestRateProbInvariant=currRateProbInvariant;
				_bestL=currBestL;			
			}
		}
// optimization - Freq (Theta)
		if (!gainLossOptions::_isStartWithTheta && evalTheta && !gainLossOptions::_isRootFreqEQstationary){
			currBestL = -brent(MINIMUM_FREQ_PARAM,_bestTheta,MAXIMUM_FREQ_PARAM,C_evalParam(tr,sp,sc,C_evalParam::theta,isReversible,_weightsUniqPatterns,_unObservableData_p),
				epsilonOptimization*gainLossOptions::_epsilonOptimizationThetaFactor,&currTheta);
			if (currBestL>_bestL) {
				static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setTheta(currTheta);
				sumPijQij = normalizeQ(&sp); //TEST
				if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp);
				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tTheta= "<<currTheta<<endl);
				_bestTheta=currTheta;
				_bestL=currBestL;
				//MDOUBLE currentlogL =likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,sp,_weightsUniqPatterns,_unObservableData_p);
				//if(!DEQUAL(currentlogL,_bestL)){ //DEQUAL(currentlogL,bestL)
				//	LOGnOUT(3,<<"!!! ERROR: different likelihood after optimizeGainLossModel,diff= "<<currentlogL-_bestL <<"\n");
				//}
			}
		}
		//else{
		//	_bestTheta = _bestMu1/(_bestMu1+_bestMu2);
		//}

		if (!(_bestL>previousL+epsilonOptimizationIter))	// no significant improvement -> break 
		{
			//if(_unObservableData_p) _unObservableData_p->setLforMissingData(tr,&sp); // Done after each update, not here
			//_bestL=max(_bestL,currBestL);	// not to reduce likelihood
			LOGnOUT(3,<<" model optimization converged. Iter= "<<iter<<" Likelihood="<<_bestL<<endl);
			break;
		}
		if(gainLossOptions::_simulatedAnnealing)
			epsilonOptimization = max(epsilonOptimization*gainLossOptions::_simulatedAnnealingCoolingFactor,0.1*gainLossOptions::_simulatedAnnealingMinEpsilonFactor);  //simulated annealing
	}
	if (iter>=numIterations){
		_bestL=max(_bestL,currBestL);	// not to reduce likelihood
		LOGnOUT(3,<<" Too many iterations in optimizeGainLossModel. Iter= "<<iter<<" Last optimized parameters are used. iter="<<iter<<endl);
	}
	//if(currUnObservableData_p) delete currUnObservableData_p;
}

/********************************************************************************************
*********************************************************************************************/
//bool optimizeGainLossModel::isUpdateGain(const MDOUBLE currBestL, MDOUBLE& currM1, const MDOUBLE lossLikelihoodImprovmet)
//{
//	bool isUpdateGain = false;
//	if((currBestL-_bestL)>lossLikelihoodImprovmet){
//		isUpdateGain =true;
//	}
//	return isUpdateGain;
//}

/********************************************************************************************
*********************************************************************************************/
//void optimizeGainLossModel::initMissingDataInfo()
//{
//	//if(gainLossOptions::_accountForMissingData && (_plogLforMissingData==NULL)){ // plogLforMissingData was't sent but it is needed
//	//	LOGnOUT(4,<<"----------plogLforMissingData was't sent but it is needed"<<endl);
//	//	_plogLforMissingData = &_logLforMissingData;
//	//}
//}



/********************************************************************************************
*********************************************************************************************/
//optimizeGainLossModel::optimizeGainLossModel(const tree& tr, stochasticProcess& sp, const sequenceContainer &sc,
//					   const bool isReversible, /*const bool evalTheta,*/ 
//					   const MDOUBLE epsilonOptimization, const int numIterations,
//					   MDOUBLE* plogLforMissingData,
//					   /*const MDOUBLE upperValueOfParam, const MDOUBLE lowerValueOfParam,*/
//					   ostream& out):
//_plogLforMissingData(plogLforMissingData)
//{
//	//initMissingDataInfo();
//	
//	MDOUBLE MINIMUM_ALPHA_PARAM;
//	if(gainLossOptions::_isAlphaLimit){
//		MINIMUM_ALPHA_PARAM = 0.3;
//	}
//	else{
//		MINIMUM_ALPHA_PARAM = ::MINIMUM_ALPHA_PARAM;
//	}
//	MDOUBLE MAXIMUM_GAIN_PARAM;
//	if(gainLossOptions::_gainLossRateAreFreq){
//		MAXIMUM_GAIN_PARAM = 0.9999;
//	}
//	else{
//		MAXIMUM_GAIN_PARAM = ::MAXIMUM_GAIN_PARAM;
//	}
//
//
//	bool optimizeAlpha = isAlphaOptimization(sp.distr());
//	bool optimizeBeta = isBetaOptimization(sp.distr());
//	bool optimizeMixture = isMixOptimization(sp.distr());
//	bool probInvariant = isInvariantOptimization(sp.distr());
//	bool evalTheta = isThetaOptimization();
//	
//	MDOUBLE previousL;
//	MDOUBLE currBestL=VERYSMALL;
//	MDOUBLE currM1=0.1;
//	MDOUBLE currM2=1; // for non-reversible model only
//	MDOUBLE currAlpha=1;
//	MDOUBLE currBeta=1;
//	MDOUBLE currTheta = 0.5;
//	MDOUBLE currRateProbInvariant = 0.05;
//	MDOUBLE lossLikelihoodImprovmet = 0;
//	MDOUBLE incrementFactorForGain = gainLossOptions::_slopeFactorForGain;	// forces slow climb for gain param
//	MDOUBLE currLogLforMissingData;
//	MDOUBLE* currpLogLforMissingData;
//	if(gainLossOptions::_accountForMissingData){
//		currpLogLforMissingData = &currLogLforMissingData;
//		*currpLogLforMissingData = *_plogLforMissingData;
//	}
//	else
//		currpLogLforMissingData = NULL;	
//
//
//	
//	if(gainLossOptions::_initParamsAtRandPointsInOptimization){
//		currM1 =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_GAIN_PARAM, MAXIMUM_GAIN_PARAM);
//		currM2=talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_LOSS_PARAM, MAXIMUM_LOSS_PARAM);
//		currAlpha = talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_ALPHA_PARAM, MAXIMUM_ALPHA_PARAM);
//		currBeta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_BETA_PARAM, MAXIMUM_BETA_PARAM);
//		currTheta =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM); 
//		currRateProbInvariant =talRandom::giveRandomNumberBetweenTwoPoints(MINIMUM_PROB_PARAM, MAXIMUM_PROB_PARAM);
//	}
//
//// initialize
//	_bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,sp,0,currpLogLforMissingData);
//	_bestMu1 = static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->getMu1();
//	if (!isReversible){
//		_bestMu2 = static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->getMu2();	}
//	if(optimizeAlpha){
//		_bestAlpha = getRateAlpha(sp.distr());	}
//	if(optimizeBeta){
//		_bestBeta = getRateBeta(sp.distr());	}
//	if(evalTheta){
//		_bestTheta = static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->getTheta();	}
//	if(probInvariant){
//		_bestRateProbInvariant = static_cast<generalGammaDistributionPlusInvariant*>(sp.distr())->getInvProb();	}
//
//// optimize
//	LOGnOUT(3,<<"### "<<"optimization starting- 'epsilonOptimization'="<<epsilonOptimization<<" 'numIterations'="<<numIterations<<endl);
//	LOGnOUT(3,<<"start optimization with:"<<endl<<" L= "<<_bestL<<" gain= "<<_bestMu1);
//	if(!isReversible) LOGnOUT(3,<<" loss= "<<_bestMu2);
//	if(optimizeAlpha) LOGnOUT(3,<<" Alpha= "<<_bestAlpha);
//	if(optimizeBeta) LOGnOUT(3,<<" Beta= "<<_bestBeta);
//	if(optimizeMixture) LOGnOUT(3,<<" ");
//	if(evalTheta)	LOGnOUT(3,<<" Theta= "<<_bestTheta);
//	if(probInvariant)	LOGnOUT(3,<<" RateProbInvariant= "<<_bestRateProbInvariant);
//	LOGnOUT(3,<<endl);
//	
//	int iter;
//	for (iter=1;iter<=numIterations;iter++){
//		previousL = _bestL;
//		//bool changed=false;
//		LOGnOUT(4,<<"iter="<<iter<<endl);
//// optimization - gain
//		if(gainLossOptions::_incrementFactorForGain) currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,min((_bestMu1*incrementFactorForGain),MAXIMUM_GAIN_PARAM),C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,currpLogLforMissingData),epsilonOptimization,&currM1);		
//		if(gainLossOptions::_lossBiggerGainLimit) currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,min(_bestMu2,MAXIMUM_GAIN_PARAM),C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,currpLogLforMissingData),epsilonOptimization,&currM1);		
//		else currBestL = -brent(MINIMUM_GAIN_PARAM,_bestMu1,MAXIMUM_GAIN_PARAM,C_evalParam(tr,sp,sc,C_evalParam::gain,isReversible,currpLogLforMissingData),epsilonOptimization,&currM1);		
//		if (currBestL>_bestL) {
//			//lossLikelihoodImprovmet *= requiredPresentOflastLikelihoodImprovmet;
//			//if (isUpdateGain(currBestL,currM1,lossLikelihoodImprovmet)) {
//			static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setMu1(currM1,isReversible);
//			*_plogLforMissingData = *currpLogLforMissingData;
//			LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tGain= "<<currM1<<endl);
//			_bestMu1=currM1;
//			_bestL=currBestL;			
//		}
//		//if (currBestL>_bestL+epsilonOptimization) {
//		////if (isUpdateGain(currBestL,currM1,lossLikelihoodImprovmet)) {
//		//	changed=true;
//		//	_bestL=currBestL;			
//		//}
//
//// optimization - loss
//		if (!isReversible & !gainLossOptions::_gainLossRateAreFreq){
//			if(gainLossOptions::_lossBiggerGainLimit) currBestL = -brent(max(_bestMu1,MINIMUM_LOSS_PARAM),_bestMu2,MAXIMUM_LOSS_PARAM,C_evalParam(tr,sp,sc,C_evalParam::loss,isReversible,currpLogLforMissingData),epsilonOptimization,&currM2);
//			else currBestL = -brent(MINIMUM_LOSS_PARAM,_bestMu2,MAXIMUM_LOSS_PARAM,C_evalParam(tr,sp,sc,C_evalParam::loss,isReversible,currpLogLforMissingData),epsilonOptimization,&currM2);
//			if (currBestL>_bestL) {
//				lossLikelihoodImprovmet = currBestL-_bestL;
//				static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->setMu2(currM2);
//				*_plogLforMissingData = *currpLogLforMissingData;
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tLoss= "<<currM2<<endl);
//				_bestMu2=currM2;
//				_bestL=currBestL;
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//
//// optimize Alpha - 3 options  (all results with same values)
//		if(optimizeAlpha){
//			currBestL = -brent(MINIMUM_ALPHA_PARAM,_bestAlpha,MAXIMUM_ALPHA_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateAlpha,isReversible,currpLogLforMissingData),epsilonOptimization,&currAlpha);			
//			if (currBestL>_bestL) {
//				//static_cast<gammaDistribution*>(sp.distr())->setAlpha(currAlpha);
//				setRateAlpha(sp.distr(),currAlpha);
//				*_plogLforMissingData = *currpLogLforMissingData;
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tAlpha= "<<currAlpha<<endl);
//				_bestAlpha=currAlpha;
//				_bestL=currBestL;			
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//// optimize Beta 
//		if(optimizeBeta){
//			currBestL = -brent(MINIMUM_BETA_PARAM,_bestBeta,MAXIMUM_BETA_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateBeta,isReversible,currpLogLforMissingData),epsilonOptimization,&currBeta);			
//			if (currBestL>_bestL) {
//				setRateBeta(sp.distr(),currBeta);
//				*_plogLforMissingData = *currpLogLforMissingData;
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tBeta= "<<currBeta<<endl);
//				_bestBeta=currBeta;
//				_bestL=currBestL;			
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//// optimize mixture
//		if(optimizeMixture){
//			//Vint pointsNum(1, 1);
//			//Vint iterNum(1, 1), const vector<GamMixtureOptimizer::OptimAlg> optAlgs, const Vdouble tols
//			
//			GamMixtureOptimizer optGamma(&sp, sc, tr);
//			if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::EM) currBestL  = optGamma.findBestParam(GamMixtureOptimizer::EM, 1, epsilonOptimization, NULL);
//			else if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::ONE_DIM) currBestL  = optGamma.findBestParam(GamMixtureOptimizer::ONE_DIM, 1, epsilonOptimization, NULL);
//			else errorMsg::reportError("unknown type in gammmaMixtureOptimizerAlgType");
//
//			if (currBestL>_bestL) {
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\timprovment in optimize gammaMixture params"<<endl);
//				_bestL=currBestL;
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//// optimization - Freq (Theta)
//		if (evalTheta){
//			currBestL = -brent(MINIMUM_PROB_PARAM,_bestTheta,MAXIMUM_PROB_PARAM,C_evalParam(tr,sp,sc,C_evalParam::theta,isReversible,currpLogLforMissingData),epsilonOptimization,&currTheta);
//			if (currBestL>_bestL) {
//				static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setTheta(currTheta);
//				*_plogLforMissingData = *currpLogLforMissingData;
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tTheta= "<<currTheta<<endl);
//				_bestTheta=currTheta;
//				_bestL=currBestL;			
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//// optimization - Prob
//		if (probInvariant){
//			currBestL = -brent(MINIMUM_PROB_PARAM,_bestRateProbInvariant,MAXIMUM_PROB_PARAM,C_evalParam(tr,sp,sc,C_evalParam::rateProbInvariant,isReversible,currpLogLforMissingData),epsilonOptimization,&currRateProbInvariant);
//			if (currBestL>_bestL) {
//				static_cast<generalGammaDistributionPlusInvariant*>(sp.distr())->setInvProb(currRateProbInvariant);
//				*_plogLforMissingData = *currpLogLforMissingData;
//				LOGnOUT(4,<<"currBestL= "<<currBestL<<"\tRateProbInvariant= "<<currRateProbInvariant<<endl);
//				_bestRateProbInvariant=currRateProbInvariant;
//				_bestL=currBestL;			
//			}
//			//if (currBestL>_bestL+epsilonOptimization) {
//			//	changed=true;
//			//	_bestL=currBestL;			
//			//}
//		}
//		if (!(_bestL>previousL+epsilonOptimization))	// no significant improvement -> break 
//		{
//			_bestL=max(_bestL,currBestL);	// not to reduce likelihood
//			break;
//		}
//	}
//	if (iter>=numIterations){
//		_bestL=max(_bestL,currBestL);	// not to reduce likelihood
//		LOGnOUT(3,<<"WARNING: Too many iterations in optimizeGainLossModel. Last optimized parameters are used. iter="<<iter<<endl);
//	}
//}







