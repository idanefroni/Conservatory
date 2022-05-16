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


#ifndef ___GAIN_LOSS_OPTIMIZER
#define ___GAIN_LOSS_OPTIMIZER

#include "definitions.h"
#include "replacementModel.h"
#include "gainLoss.h"
#include "gainLossOptions.h"
#include "mixtureDistribution.h"
#include "unObservableData.h"
/********************************************************************************************
The optimization flow:
*note: the gainLossOptimizer changes "in place" (byRef)
*note: the C_evalParam makes a copy

gainLoss (-> startOptimizations -> optimizationsManyStarts[/optimizationsVVManyStarts] )
	gainLossOptimizer//overloaded for spVVec (-> optimizations[/optimizationsSPvv] 
													-> optimizeBranchLengths[/optimizeBranchLengthsSpvv]
													-> optimizeParameters [/optimizeParametersSPvv])
			optimizeGainLossModel (->brent)
				C_evalParam (->setParam)
					likelihoodComputation

*********************************************************************************************/

/********************************************************************************************
gainLossOptimizer
*********************************************************************************************/
class gainLossOptimizer
{

public:
	explicit gainLossOptimizer(tree& tr, stochasticProcess* sp, const sequenceContainer &sc, 
		const MDOUBLE epsilonOptimization, const int numIterations,
		const MDOUBLE epsilonOptimizationModel, const int numIterationsModel,
		const MDOUBLE epsilonOptimizationBBL, const int numIterationsBBL,
		Vdouble * weights,
		unObservableData* unObservableData_p, bool performOptimizationsBBL, bool isbblLSWhenbblEMdontImprove);

	explicit gainLossOptimizer(tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist, 
		const sequenceContainer &sc, 
		const MDOUBLE epsilonOptimization, const int numIterations,
		const MDOUBLE epsilonOptimizationModel, const int numIterationsModel,
		const MDOUBLE epsilonOptimizationBBL, const int numIterationsBBL,
		Vdouble * weights,
		unObservableData* _unObservableData_p, bool performOptimizationsBBL, bool isbblLSWhenbblEMdontImprove); 	


	virtual ~gainLossOptimizer(){;}
	
	MDOUBLE getBestL(){return _bestL;}
	tree getOptTree(){return _tr;};
	gainLossOptions::distributionType getRateDistributionType(distribution* dist);


protected:
	//func
	//void initMissingDataInfo();
	void optimizations();
	void optimizationsSPvv();

	MDOUBLE  optimizeParameters();
	MDOUBLE  optimizeParametersSPvv();

	MDOUBLE optimizeBranchLengths(const int outerIter);
	MDOUBLE optimizeBranchLengthsvv(const int outerIter);

	MDOUBLE optimizeRoot();
	MDOUBLE optimizeRootSPvv();

	void printMixtureParams();

protected:
	//members
	MDOUBLE _bestL;
	MDOUBLE _epsilonOptimization;
	int _maxNumOfIterations;
	MDOUBLE _epsilonOptimizationModel;
	int _maxNumOfIterationsModel;
	MDOUBLE _epsilonOptimizationBBL;
	int _maxNumOfIterationsBBL;
	////MDOUBLE _logLforMissingData;
	//MDOUBLE* _plogLforMissingData;		
	//Vdouble* _pLforMissingDataPerCat;	// used foreach rate category
	unObservableData*  _unObservableData_p;

	Vdouble* _weightsUniqPatterns;

	bool _performOptimizationsBBL;
	bool _isbblLSWhenbblEMdontImprove;

	stochasticProcess *_sp;
	MDOUBLE _bestGain;
	MDOUBLE _bestLoss;
	MDOUBLE _bestAlphaRate;
	MDOUBLE _bestBetaRate;
	MDOUBLE _bestRateProbInvariant;
	stochasticProcess *_spSimple;
	MDOUBLE _bestTheta;
	Vdouble _freq;
	MDOUBLE _bestGainAlpha;
	MDOUBLE _bestGainBeta;
	MDOUBLE _bestGainProbInvariant;

	MDOUBLE _bestLossAlpha;
	MDOUBLE _bestLossBeta;
	MDOUBLE _bestLossProbInvariant;

	MDOUBLE _gainExp;
	MDOUBLE _lossExp;

	MDOUBLE _gainSTD;
	MDOUBLE _lossSTD;

	bool _isReversible;
	bool _isSkipBblEM;
	tree _tr;
	sequenceContainer _sc;

	vector<vector<stochasticProcess*> > _spVVec; //save stochasticProcess for each category
	distribution* _gainDist;
	distribution* _lossDist;
	gainLossOptions::distributionType _rateDistributionType;

};


#endif
