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
/********************************************************************************************
gainLossOptions - a class that contains all the parameters for the gainLossProjest as static
	use the 'Parameters' class to read info from txt file.
	initDefault. (+Parameters::addParameter)
	getParamsFromFile. ->with alterations of defults for consistancy
	verifyConsistParams.
*********************************************************************************************/
#include "gainLossOptions.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "Parameters.h"
#include <iostream>
#include <cmath>

using namespace std;

// recognize all the static members defined at .h
int gainLossOptions::_alphabet_size;
string gainLossOptions::_seqFile;
string gainLossOptions::_treeFile;
string gainLossOptions::_treeFileOrig;	// used for branchDiff calc. functionality

string gainLossOptions::_rootAt;
string gainLossOptions::_logFile;
int gainLossOptions::_logValue;
string gainLossOptions::_referenceSeq; 
string gainLossOptions::_outDir;
string gainLossOptions::_treeOutFile;
//string gainLossOptions::_outFile;
//string gainLossOptions::_outFileNotNormalize;
//string gainLossOptions::_outFileGain4Site;
//string gainLossOptions::_outFileLoss4Site;
//string gainLossOptions::_outFileLikeofPos;
//string gainLossOptions::_outFilePosteriorExpectationOfChange;
//gainLossOptions::discretizationType gainLossOptions::_discretizationType;
gainLossOptions::treeSearchAlgType gainLossOptions::_treeSearchAlg;
gainLossOptions::gammmaMixtureOptimizerAlgType gainLossOptions::_gammmaMixtureOptimizerAlg;
gainLossOptions::distributionType gainLossOptions::_gainDistributionType; 
gainLossOptions::distributionType gainLossOptions::_lossDistributionType;
gainLossOptions::distributionType gainLossOptions::_rateDistributionType;
gainLossOptions::rateEstimationMethodType gainLossOptions::_rateEstimationMethod;
gainLossOptions::characterFreqEvalType gainLossOptions::_characterFreqEval;
gainLossOptions::discretizationType gainLossOptions::_rateDiscretizationType;
MDOUBLE gainLossOptions::_userGainLossRatio;
bool gainLossOptions::_keepUserGainLossRatio;
MDOUBLE gainLossOptions::_userAlphaGain;
MDOUBLE gainLossOptions::_userBetaGain;
MDOUBLE gainLossOptions::_userProbInvariantGain;
MDOUBLE gainLossOptions::_userAlphaLoss;
MDOUBLE gainLossOptions::_userBetaLoss;
MDOUBLE gainLossOptions::_userProbInvariantLoss;
MDOUBLE gainLossOptions::_userProbInvariantRate;
MDOUBLE gainLossOptions::_userRateInvariantVal;
MDOUBLE gainLossOptions::_userAlphaRate;
MDOUBLE gainLossOptions::_userBetaRate;
MDOUBLE gainLossOptions::_userGain;
MDOUBLE gainLossOptions::_userLoss;
MDOUBLE gainLossOptions::_userTheta;
MDOUBLE gainLossOptions::_userAlphaGainMax;
MDOUBLE gainLossOptions::_userBetaGainMax;
MDOUBLE gainLossOptions::_userProbInvariantGainMax;
MDOUBLE gainLossOptions::_userAlphaLossMax;
MDOUBLE gainLossOptions::_userBetaLossMax;
MDOUBLE gainLossOptions::_userProbInvariantLossMax;
MDOUBLE gainLossOptions::_userProbInvariantRateMax;
MDOUBLE gainLossOptions::_userAlphaRateMax;
MDOUBLE gainLossOptions::_userBetaRateMax;
MDOUBLE gainLossOptions::_userGainMax;
MDOUBLE gainLossOptions::_userLossMax;
MDOUBLE gainLossOptions::_userThetaMax;
MDOUBLE gainLossOptions::_userAlphaGainMin;
MDOUBLE gainLossOptions::_userBetaGainMin;
MDOUBLE gainLossOptions::_userProbInvariantGainMin;
MDOUBLE gainLossOptions::_userAlphaLossMin;
MDOUBLE gainLossOptions::_userBetaLossMin;
MDOUBLE gainLossOptions::_userProbInvariantLossMin;
MDOUBLE gainLossOptions::_userProbInvariantRateMin;
MDOUBLE gainLossOptions::_userAlphaRateMin;
MDOUBLE gainLossOptions::_userBetaRateMin;
MDOUBLE gainLossOptions::_userGainMin;
MDOUBLE gainLossOptions::_userLossMin;
MDOUBLE gainLossOptions::_userThetaMin;

MDOUBLE gainLossOptions::_probCutOffPrintEvent;
bool gainLossOptions::_isFewCutOffCounts;
MDOUBLE gainLossOptions::_probCutOffCounts;

int gainLossOptions::_numberOfGainCategories;
int gainLossOptions::_numberOfLossCategories;
int gainLossOptions::_numberOfRateCategories;
int gainLossOptions::_numberOfRateComponents;
int gainLossOptions::_maxNumOfIterations;
int gainLossOptions::_maxNumOfIterationsModel; 
int gainLossOptions::_maxNumOfIterationsBBL;
int gainLossOptions::_maxNumOfIterationsManyStarts;
int gainLossOptions::_numberOfRandPointsInOptimization;
int gainLossOptions::_numberOfRandStartPoints;

int gainLossOptions::_numOfSimulationsForPotExp;

gainLossOptions::optimizationLevel gainLossOptions::_optimizationLevel;
MDOUBLE gainLossOptions::_epsilonOptimizationModel;
MDOUBLE gainLossOptions::_epsilonOptimizationBBL;

MDOUBLE gainLossOptions::_epsilonOptimizationIterationCycleManyStarts;
MDOUBLE gainLossOptions::_epsilonFactor_Model;
MDOUBLE gainLossOptions::_epsilonFactor_BBL; 
MDOUBLE gainLossOptions::_numIterationsFactor_Model; 
MDOUBLE gainLossOptions::_numIterationsFactor_BBL; 

MDOUBLE gainLossOptions::_epsilonOptimizationIterationCycle;
bool gainLossOptions::_gainLossDist;
bool gainLossOptions::_calculateRate4site;
bool gainLossOptions::_calculeGainLoss4site;
MDOUBLE gainLossOptions::_likelihoodLandscapeIncrement;
bool gainLossOptions::_printLikelihoodLandscape;
bool gainLossOptions::_printLikelihoodLandscapeAlphaRate;
bool gainLossOptions::_printLikelihoodLandscapeGainLoss;
bool gainLossOptions::_printLikelihoodLandscapeTheta;
bool gainLossOptions::_optAlphaInIteration;
bool gainLossOptions::_optBBL_LS_InIteration;
bool gainLossOptions::_optBBL_EM_InIteration;
bool gainLossOptions::_printP11forgain;
bool gainLossOptions::_printTree;
bool gainLossOptions::_printSeq;
bool gainLossOptions::_printPij_t;
bool gainLossOptions::_printLofPos;
bool gainLossOptions::_printLofPosBothModels;
bool gainLossOptions::_performOptimizations;
bool gainLossOptions::_correctOptimizationEpsilon;
bool gainLossOptions::_performOptimizationsBBL;
bool gainLossOptions::_performOptimizationsBBLOnlyOnce;

bool gainLossOptions::_isBblLS;
bool gainLossOptions::_isbblLSWhenbblEMdontImprove;
bool gainLossOptions::_isSkipBblEMWhenbblEMdontImprove;

bool gainLossOptions::_isInitGainLossByEmpiricalFreq;
bool gainLossOptions::_isBBLEMwithSimpleSpBeforeFullOptimization;
bool gainLossOptions::_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately;
bool gainLossOptions::_isOptimizeInvariantCategoryProb;

bool gainLossOptions::_isUpdateOnlyGainBetaForRatio;
bool gainLossOptions::_isComputeLikelihoodDuringInit;

bool gainLossOptions::_isBblEMbeforeLSWithMissSpecifiedModel;
bool gainLossOptions::_isBblForceFactorCorrection;
MDOUBLE gainLossOptions::_BblFactorCorrection;

bool gainLossOptions::_isSkipFirstParamsOptimization;
bool gainLossOptions::_isOptimizeParamsWithLogMinMax;

bool gainLossOptions::_isMultipleAllBranchesByFactorAtStart;
bool gainLossOptions::_isNormalizeAtStart;

bool gainLossOptions::_performOptimizationsROOT;
bool gainLossOptions::_performOptimizationsBBLManyStarts;
bool gainLossOptions::_simulatedAnnealing;
MDOUBLE gainLossOptions::_simulatedAnnealingMinEpsilonFactor;
MDOUBLE gainLossOptions::_simulatedAnnealingCoolingFactor;
bool gainLossOptions::_performOptimizationsManyStarts;
bool gainLossOptions::_gainLossDistPlusInvariant;
bool gainLossOptions::_isHGT_normal_Pij;
bool gainLossOptions::_isHGT_with_Q;
bool gainLossOptions::_initParamsAtRandPoints;
bool gainLossOptions::_initParamsAtRandPointsInOptimization;
bool gainLossOptions::_calculePosteriorExpectationOfChange;
bool gainLossOptions::_simulatePosteriorExpectationOfChange;
bool gainLossOptions::_isOnlySimulateSeq;

bool gainLossOptions::_modelOptimizationSimPostExp;
bool gainLossOptions::_BBLOptimizationSimPostExp;
bool gainLossOptions::_initParamsAtRandPointsInSimPostExp;
bool gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos;
bool gainLossOptions::_isFlatTreeBeforOpt;
bool gainLossOptions::_isbBLEMwithSimpleSpSimulatePostExp;
MDOUBLE gainLossOptions::_noiseLevelInGammaSimulation;
bool gainLossOptions::_isTheataFromObservedFreq;
bool gainLossOptions::_isRootFreqEQstationaryInSimulations;
bool gainLossOptions::_isMatrixGainLossFromRatioInSimulations;
bool gainLossOptions::_isFlatSpBeforeOpt;
MDOUBLE gainLossOptions::_epsilonOptForPostExpSimFactor;
MDOUBLE gainLossOptions::_numOfIterationsOptForPostExpSimFactor;
MDOUBLE gainLossOptions::_loss2gainRatioToSim;

bool gainLossOptions::_printAncestralReconstructPosterior;
bool gainLossOptions::_saveProbChanges_PosNodeXY;
bool gainLossOptions::_isComputeDistanceFromRootForRecent;

bool gainLossOptions::_printTreesWithProbabilityValuesAsBP;
bool gainLossOptions::_printTreesWithExpectationValuesAsBP;
bool gainLossOptions::_calculateAncestralReconstruct;
bool gainLossOptions::_printTreesWithAncestralReconstructAsBP;
bool gainLossOptions::_printAncestralReconstructFullData;

bool gainLossOptions::_printDEBUGinfo;
bool gainLossOptions::_printPropExpOfChangeFullData;
bool gainLossOptions::_printExpPerPosPerBranchMatrix;
bool gainLossOptions::_printComputedCorrelations;
bool gainLossOptions::_performParametricBootstapCorrelation;
bool gainLossOptions::_usePosSpecificSimulations;
bool gainLossOptions::_isConsiderNegativeCorrelations;
bool gainLossOptions::_isDivideBinsByRange;
bool gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound;
bool gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound;
MDOUBLE gainLossOptions::_relativeSizeOfOverLappedBins;

bool gainLossOptions::_isPrintpairWiseCorrelationsAndNmin;
bool gainLossOptions::_isPrintCorrelationsOfAllPairs_Corr;
bool gainLossOptions::_isPrintCorrelationsOfAllPairs_pVal;
bool gainLossOptions::_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH;
bool gainLossOptions::_isAllCorrTypeReqruiedToBeSignificant;
bool gainLossOptions::_isNminBasedOnCountBranchesOverCutOff;

int gainLossOptions::_numOfBinsInParametricBootstrapSimulations;
bool gainLossOptions::_isAddSimulationsWithLowRate;
bool gainLossOptions::_isFDRcorrectionForPValInCorrelation;
bool gainLossOptions::_isComputeQVals;
MDOUBLE gainLossOptions::_pValueCutOffForBootStrap;
MDOUBLE gainLossOptions::_minExpThresholdForPValComputationForCorrelatingPair;
bool gainLossOptions::_isUpdateMinExpThresholdGivenSimulaitonsQuantile;
bool gainLossOptions::_isUpdateMinExpThresholdGivenRealDataQuantile;
MDOUBLE gainLossOptions::_updateMinExpThresholdGivenRealDataQuantileVal;

bool gainLossOptions::_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel;
bool gainLossOptions::_isCompExtremeValDistribution;

MDOUBLE gainLossOptions::_minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair;

bool gainLossOptions::_isCorrelateWithPearson;
bool gainLossOptions::_isCorrelateWithSpearman;
bool gainLossOptions::_isCorrelationsBasedOnMaxParsimonyMapping;

bool gainLossOptions::_isAlsoCorrelateWithLoss;
bool gainLossOptions::_isAlsoCorrelateWithBoth;
bool gainLossOptions::_isOnlyCorrelateWithBoth;
bool gainLossOptions::_isUseRateForSiteAsNminForCorrelations;

bool gainLossOptions::_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur;
bool gainLossOptions::_isRemoveSimulatedPositionsBasedOnMP;
MDOUBLE gainLossOptions::_minNumOfMPEvent2RemoveSimulatedPositions;
bool gainLossOptions::_isUpdateminNumOfMPEvent2RemoveSimulatedPositions;


bool gainLossOptions::_printComputedCorrelationsAllSites;
bool gainLossOptions::_isIgnoreCorrelationAmongSelectedSites;
bool gainLossOptions::_isNormalizeForBranchExpInCorrCompute;
bool gainLossOptions::_isNormalizeByExpectationPerBranch;
string gainLossOptions::_selectedSitesForCorrelation;
bool gainLossOptions::_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation;
int gainLossOptions::_checkCoEvolWithUnionPAP_against_pos;


bool gainLossOptions::_isReversible;
bool gainLossOptions::_isRootFreqEQstationary;
bool gainLossOptions::_initRandomGammaMixuteParam;
bool gainLossOptions::_incrementFactorForGain;
bool gainLossOptions::_lossBiggerGainLimit;
MDOUBLE gainLossOptions::_slopeFactorForGain;

bool gainLossOptions::_isStartWithTheta;
bool gainLossOptions::_isSkipGainOptimization;
MDOUBLE gainLossOptions::_epsilonOptimizationThetaFactor;
bool gainLossOptions::_isAlphaLimit;
bool gainLossOptions::_isGainLimit;
//MDOUBLE gainLossOptions::_probCutOffSum;
MDOUBLE gainLossOptions::_maxRateForML;
MDOUBLE gainLossOptions::_minBranchLength;
MDOUBLE gainLossOptions::_maxBranchLength;
MDOUBLE gainLossOptions::_epsilonForReRootFactor;
MDOUBLE gainLossOptions::_percentOfImprovManySarts;
MDOUBLE gainLossOptions::_percentOfImprov;

bool gainLossOptions::_calculeBranchLegthDiffFactor;

gainLossOptions::simulationType gainLossOptions::_simulationType;
bool gainLossOptions::_isMPratio;
bool gainLossOptions::_isInitGainLossByEmpiricalFreqSimulatePostExp;
bool gainLossOptions::_is3states;
MDOUBLE gainLossOptions::_3statesGain;
MDOUBLE gainLossOptions::_3statesMore;
MDOUBLE gainLossOptions::_3statesLess;
MDOUBLE gainLossOptions::_3statesLoss;
MDOUBLE gainLossOptions::_3states0;
MDOUBLE gainLossOptions::_3states1;

bool gainLossOptions::_simulateSequences;
bool gainLossOptions::_isReversibleSim;
bool gainLossOptions::_useTheSameSpForSim;
int gainLossOptions::_numberOfSequences2simulate;
int gainLossOptions::_numberOfPositions2simulate;
int gainLossOptions::_numberOfIterations2simulate;
int gainLossOptions::_numberOfIterationsForPrintResults;
MDOUBLE gainLossOptions::_percentileOfNminWithCorr1RequiredForLastIteration;


gainLossOptions::distributionType gainLossOptions::_rateDistributionTypeSim;
bool gainLossOptions::_gainEQlossSim;
bool gainLossOptions::_calculateRate4siteSim;
bool gainLossOptions::_writeSeqSim;
bool gainLossOptions::_accountForMissingData;
bool gainLossOptions::_gainEQloss;
bool gainLossOptions::_gainLossRateAreFreq;
bool gainLossOptions::_findCoEvolvingSitesOldNotWorking;					// for the co evolving project
int gainLossOptions::_numberOfSequences2simulateForCoEvol;	// for the co evolving project
Vdouble* gainLossOptions::_weights;
int gainLossOptions::_minNumOfOnes;
int gainLossOptions::_minNumOfZeros;

ostream* gainLossOptions::_outPtr;
bool gainLossOptions::_isAnaliticComputeJumps;
bool gainLossOptions::_isSequenceUniqPattern;
bool gainLossOptions::_isRemovePositionsWithHighPercentOfMissingData;
MDOUBLE gainLossOptions::_fractionOfMissingDataToRemove;

bool gainLossOptions::_isOnlyComputeLikelihood;
bool gainLossOptions::_isNormalizeQ;
bool gainLossOptions::_isNormalizeQinSpVVec;
bool gainLossOptions::_isNormalizeQandTreeafterOpt;
bool gainLossOptions::_isFlatUserParameters;
bool gainLossOptions::_isAlphaEqBetaManipulation;
bool gainLossOptions::_calculeBranchLegthDiffFactorFromInputTrees;
bool gainLossOptions::_intersectTreeAndSeq;

bool gainLossOptions::_isOnlyParsimony;
bool gainLossOptions::_calculeMaxParsimonyChange;
bool gainLossOptions::_calculeMaxParsimonyChangeSeveralGainLossRatios;
string gainLossOptions::_costMatrixfile;
gainLossOptions::costMatrixType  gainLossOptions::_costMatrixType;
MDOUBLE gainLossOptions::_costMatrixGainLossRatio;


//ofstream gainLossOptions::_out_f;
//string gainLossOptions::_mainType;


gainLossOptions::~gainLossOptions(){}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::initOptions(const string& paramFileName)
{	
	getOutDirFromFile(paramFileName);	// first set _outDir to be used next
	createDir("", gainLossOptions::_outDir);
	initDefault();
	getParamsFromFile(paramFileName);
	verifyConsistParams();
}

/********************************************************************************************
initDefault
*********************************************************************************************/
void gainLossOptions::initDefault()
{
// all the default values are stored in the gainLossOptions:: static members
//################### Basic parameters:
// input (general)
	_seqFile = "";				// essential - fasta file with presence(1)/absence(0) for each species over all gene families (positions)
	_treeFile = "";				// basic	 - if not given - calculated based on distanceTable
	_treeFileOrig = "";			// for brachLength Diff.
	_rootAt ="";				// name of node to be root (the tree must contain names of internal nodes)
	_referenceSeq = "non";		// the results are printed with this seq in each positions. (default - first)
	//static string _mainType;

// output
	//_outDir = "RESULTS";		// concatenated after current dir location 'pwd'
	_logFile = _outDir + "//" + "log.txt";		// print-outs of the running progress including the estimated parameters optimization
	_logValue = 5;				// verbosity level - ~4 - normal, >7 - load of info
	_treeOutFile = _outDir + "//" + "TheTree.ph";		// "TheTree.ph" - tree after BBL and other changes
// all of these files are still part of the output, but names are fixed
	//static string _outFile;		// Rate4Site results (normalized - Ave=0, Sd=1)
	//static string _outFileNotNormalize;		// Rate4Site results (original)
	//static string _outFileGain4Site;		// gain4Site results
	//static string _outFileLoss4Site;		// loss4Site results
	//static string _outFileLikeofPos;		// compare to model with gainRate=0
	//static string _outFilePosteriorExpectationOfChange;		// exp01, exp10 per gene

//################################################## Model params
	_alphabet_size =2;						// 2 - presence(1)/absence(0)
	_gainLossDist =false;					// GLM (mixture)
	_accountForMissingData =true;			// for phyletic patterns - must be true
	_minNumOfOnes = 1;						// for COG and EggNOG only patterns with 3 or more are observable
	_minNumOfZeros = 0;						// for indels, change to 1_isRemoveSimulatedPositionsBasedOnMP 

	_gainEQloss =false;						// M1 (the basic model - gain=loss)
	_isReversible =false;					// if(_isReversible==False) -> the root is fixed
	_isRootFreqEQstationary =true;			// same "-"
	_gainLossDistPlusInvariant =false;		// Automatically True if GENERAL_GAMMA_PLUS_INV or GAMMA_PLUS_INV
	_gainLossRateAreFreq =false;			// test parameter where gain+loss = 1, and the "r_Q" is external

//Each of the rates governing the stochastic process are assumed to be sampled from a prior distribution.	
	_rateDistributionType =GAMMA;
	_gainDistributionType =GENERAL_GAMMA; //(only for the mixture models - _gainLossDist 1) 
	_lossDistributionType =GENERAL_GAMMA; //(only for the mixture models - _gainLossDist 1)
	_numberOfGainCategories = 3;		//	gain 3-5 - the overall number of stochasticProcess 9-25
	_numberOfLossCategories = 3;		//	loss 3-5
	_numberOfRateCategories = 4;		//	discretization usually 4-16
	_numberOfRateComponents = 3;		//	gammaMix
	_rateDiscretizationType =QUANTILE;	// QUANTILE, LAGUERRE - only in use for gammaMix

//################################################## computations (What calculations are processed)
	_calculateRate4site =true;
	_rateEstimationMethod =ebExp;	// mlRate (only option for UNIFORM) or posteriorBayesianExpectation
	_calculeGainLoss4site =true;
	_calculePosteriorExpectationOfChange =true;
	_calculateAncestralReconstruct =true;
	_simulatePosteriorExpectationOfChange =false;		// simulate PostExp (To test to accuracy of the stochastic mapping)
	_isOnlySimulateSeq =false;							// no mapping or parsimony is done

	_simulateSequences =false;							// Test the rate4site computation
	_calculateRate4siteSim =false;						// Test the rate4site computation
	_calculeBranchLegthDiffFactor =true;				// if BBL is used for each branch - compare length before/after
	_findCoEvolvingSitesOldNotWorking =false;						// for the co evolving project
	_saveProbChanges_PosNodeXY =true;					// used for AnsetralReconstruc - posterior
	_isComputeDistanceFromRootForRecent =false;			// used to classify branches
	_printAncestralReconstructPosterior =true;			// huge file...
	_isOnlyParsimony = false;							// only parsimony computation and Return
	_calculeMaxParsimonyChange = true;
	_calculeMaxParsimonyChangeSeveralGainLossRatios = false;

//################################################## Prints
	_printTree =true;
	_printSeq =true;
	_printPij_t =true;
	_printLofPos =true;
	_printLofPosBothModels =false;
	_printTreesWithProbabilityValuesAsBP =false;
	_printTreesWithExpectationValuesAsBP =false;
	_printTreesWithAncestralReconstructAsBP =false;
	_printPropExpOfChangeFullData =false;				// Could be a huge file, if probCutOff is 0.0
	_printExpPerPosPerBranchMatrix =false;				// Used as input for COMAP	
	_printComputedCorrelations =false;				//
	_performParametricBootstapCorrelation =false;				
	_usePosSpecificSimulations =false;				
	_isConsiderNegativeCorrelations =false;				
	_isDivideBinsByRange =false;						// if true, each bin will get different number of samples, but the rate(Nmin) is eq-partitioned
	_isSortVectorOfCorrelationsBinsByLowerRateBound =false;
	_isSortVectorOfCorrelationsBinsByMidRateBound =true; // if true, the bins are overlapping
	_relativeSizeOfOverLappedBins = 0.25;					// if 0.25, 25% of samples per bin

	_isPrintpairWiseCorrelationsAndNmin =false;
	_isPrintCorrelationsOfAllPairs_Corr =false;	// huge files
	_isPrintCorrelationsOfAllPairs_pVal =false;	// huge files

	_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH =true; // if true, only pairs with PVal significant after BH will be printed
	_isAllCorrTypeReqruiedToBeSignificant =false; // if true, only pairs with PVal significant after BH will be printed
	_isNminBasedOnCountBranchesOverCutOff =false;	// it true, Nmin is an integer= the number of branches with probEvent>cuttoff 

	_numOfBinsInParametricBootstrapSimulations =5;				
	_isAddSimulationsWithLowRate =false;				// true seems problematics with Mixture (GL) models				
	_isFDRcorrectionForPValInCorrelation =true;
	_isComputeQVals =false;
	_pValueCutOffForBootStrap = 0.05;							// was 0.05
	_minExpThresholdForPValComputationForCorrelatingPair = 1.0;	// if 0, no Nmin is imposed, 2.0, 3.0 are possible values
	_isUpdateMinExpThresholdGivenSimulaitonsQuantile = false;	// 0.25 quantile (more "relevant" simulation)
	_isUpdateMinExpThresholdGivenRealDataQuantile = false;		// Given real data, minR is defined by the 0.1 percentile (updated only is higher) 
	_updateMinExpThresholdGivenRealDataQuantileVal = 0.1;		// if 0.2, Nmin is for sites above the 0.2 percentile rate

	_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel = false;	// elevate Nmin Threshold if: (A) freqOfHighCorr was too high (B) freqOfHighCorr is reduced consistently with higher Nmin (C) new Nmin is lower than medianNminOfRealData
	_isCompExtremeValDistribution = false;	// pValue is also estimated assuming EVD distribution

	_minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair = 1;	// e.g., if =1, with 500 species, minT = 5

	_isCorrelateWithPearson =true;			//o
	_isCorrelateWithSpearman =false;			//
	_isCorrelationsBasedOnMaxParsimonyMapping =false;			//

	_isAlsoCorrelateWithLoss =true;					// not fully functional !
	_isAlsoCorrelateWithBoth =true;					//
	_isOnlyCorrelateWithBoth =true;					// if true, only gain.concat.loss correlations are computed
	_isUseRateForSiteAsNminForCorrelations =false;	// 

	_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur =false;	// Remove simulated position with too low/high occur to save later computation time (quick and (VERY) dirty)
	_isRemoveSimulatedPositionsBasedOnMP =true;							// Remove simulated positions with less than 2 events based on max parsimony (quick and dirty)
	_minNumOfMPEvent2RemoveSimulatedPositions =1;							// If 1, then gain+loss events must be above 1 (at least one event). Must be higher for many genomes
	_isUpdateminNumOfMPEvent2RemoveSimulatedPositions =true;				// If true, add 0.2 events for every sqrt(num Of species)

	_printComputedCorrelationsAllSites =false;		//
	_isIgnoreCorrelationAmongSelectedSites =false;		// High correlation is due to shared branch length and topology
	_isNormalizeForBranchExpInCorrCompute =false;		// The values per-branch are normalized to remove branch-dependent signal
	_isNormalizeByExpectationPerBranch =true;			// else, by branch length

	_selectedSitesForCorrelation = "";				// in this file, for each position, the correlation with all other positions if computed.
	_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation = false; // the last position is a trait (with possible unknown). If true, (1) unknown removed, (2) correlation only against last
	_checkCoEvolWithUnionPAP_against_pos = 0; // if 0, not perforing union

	_printAncestralReconstructFullData =false;			// huge file...
	_printDEBUGinfo =false;								// huge file...
	_printLikelihoodLandscape =false;					// test purpose (Ad-hoc)
	_likelihoodLandscapeIncrement = 0.05;
	_printLikelihoodLandscapeAlphaRate =false;			// test purpose (Ad-hoc)
	_printLikelihoodLandscapeGainLoss =false;			// test purpose (Ad-hoc)
	_printLikelihoodLandscapeTheta =false;				// test purpose (Ad-hoc)
	_optAlphaInIteration =false;				
	_optBBL_LS_InIteration =false;				
	_optBBL_EM_InIteration =false;
	_printP11forgain =false;							// test purpose (Ad-hoc)
	
//################################################## optimizations
	_performOptimizations =true;			// model parameters are numerically estimated to maximize likelihood
	_performOptimizationsBBL = false;		// 
	_performOptimizationsBBLOnlyOnce = true;
	_isBblLS = false;						// possibly after BBL-EM, to make further improvement
	_isbblLSWhenbblEMdontImprove = false;	//If No improvement with BBL-EM -> Perform BBL-LS one iteration
	_isSkipBblEMWhenbblEMdontImprove = true; // Since no improvement, BBL-EM will be skipped next iteration, go directly to LS

	_isInitGainLossByEmpiricalFreq=true;				// the sp is initialized with the empirical 0 and 1 freq
	_isBBLEMwithSimpleSpBeforeFullOptimization=true;	// before optimization - BBL-EM is performed with simplified sp
	_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately=true;					// gain/loss is estimated (not separately gain, loss...)
	_isOptimizeInvariantCategoryProb=true;					

	_isUpdateOnlyGainBetaForRatio=false;				// work in progress...
	_isComputeLikelihoodDuringInit=true;				// true, unless fast/parsimony run is performed	

	_isBblEMbeforeLSWithMissSpecifiedModel = true;	// if both _isBblLS and this is true, after BBL-EM additional iteration is done
	_isBblForceFactorCorrection = true;
	_BblFactorCorrection = 2.0;

	_isSkipFirstParamsOptimization = false;
	_isOptimizeParamsWithLogMinMax = true;			// when the parameter is a positive and values are e.g., [0.01,100] brent works better for [-2,2]
	_isMultipleAllBranchesByFactorAtStart = true;
	_isNormalizeAtStart = true;

	_performOptimizationsROOT = false;
	_performOptimizationsManyStarts =false;			// several models are chosen are starting point for optimization
	_performOptimizationsBBLManyStarts = false;
	_correctOptimizationEpsilon =false;				// according to dataset size (was initial likelihood), abs(_logL) * gainLossOptions::_epsilonOptimizationIterationCycle *  gainLossOptions::_percentOfImprov
	_simulatedAnnealing =false;						// epsilon is lowered with iterations
	_simulatedAnnealingMinEpsilonFactor =0.2;		// lower normal epsilons (Model, BBL, Both). e.g., 0.1*0.2=0.02 - the new epsilon
	_simulatedAnnealingCoolingFactor =0.8;			// to lower epsilons each iteration

	_gammmaMixtureOptimizerAlg = ONE_DIM;	// ONE_DIM or EM (not fully functional)
	_characterFreqEval =optimizeOverTree;	// "-F option" the estimation of freq at root: FiftyFifty, LeavesAve, optimizeOverTree
	
	_isStartWithTheta =false;			// the optimization loop of the parameter will start with Theta
	_isSkipGainOptimization =false;		// 
	_epsilonOptimizationThetaFactor =1.0;			// allows for different optimization Theta

	_isAlphaLimit =true;				// 0.3 - for Alpha <<0.3, the following computations are erroneous [BUG?]
	_isGainLimit =false;					// 0.1 - for Gain  <<0.1, the following computations are erroneous [BUG?]
	_isHGT_normal_Pij =true;			// test parameter - 
	_isHGT_with_Q =true;				// test parameter - 
	_incrementFactorForGain =false;		// test parameter - 
	_lossBiggerGainLimit =false;		// test parameter - 
	_slopeFactorForGain =2.0;			// test parameter - limit growth in gain estimation
										// if the log-likelihood after optimization is lower than this threshold - then optimize again.
	_optimizationLevel = low;
	_epsilonOptimizationIterationCycle =1.0;			// 1 cycle(model+BBL) epsilon.
	_epsilonOptimizationModel =0.01;	// (was 0.05) Used by cEval for each parameter, the iteration epsilon is x3(or number of parameters)
	_epsilonOptimizationBBL =0.02;		// (was 0.1) Used by cEvel for each branch, the iteration epsilon is x5(or number of branches)
	//enum optimizationLevel {Vlow, low, mid, high, Vhigh};
	
	_epsilonOptimizationIterationCycleManyStarts = 2.0;	// 	epsilonOptimizationManyStarts = max(epsilonOptimization, abs(_logL)*gainLossOptions::_percentOfImprovManySarts);
	_percentOfImprovManySarts = 0.0001;		// epsilonOptimization = abs(logL)*_percentOfImprovManySarts
	_epsilonFactor_Model = 0.01;
	_epsilonFactor_BBL = 0.02;
	
	_maxNumOfIterationsManyStarts = 1;		// the basic number of manyStarts option (Model and BBL factors are used, 3 and 2, respectively)
	_numIterationsFactor_Model = 3;
	_numIterationsFactor_BBL = 2;

	_maxNumOfIterations = 3;				// 3
	_maxNumOfIterationsModel = 10;			// 30  
	_maxNumOfIterationsBBL = 5;				// 10

	_epsilonForReRootFactor =10;			// only for substantial improvement the tree will be re-rooted
	_percentOfImprov = 0.00001;			// for lL=-200,000 the epsilon is 0.2, epsilonOptimization = abs(logL)*_percentOfImprov*epsilonOptimization
	
	_initParamsAtRandPoints =false;
	_initParamsAtRandPointsInOptimization =true;
	_initRandomGammaMixuteParam =true;
	_numberOfRandPointsInOptimization = 10;	//10
	_numberOfRandStartPoints = 300;			//10, the loop will break before if L is improved

	//##################################################  all the model parameters can be given by the user
	_userGainLossRatio = VERYBIG;	// If given (< VERYBIG), all the related parameter are adapted
	_keepUserGainLossRatio = false;	// If given other than 1, all the related parameter are adapted
	_userGain = 0.2;		//
	_userLoss = 0.8;		//
	_userTheta =0.5;		// default 0.5 - otherwise, counting is done prior to optimization
	_userAlphaGain =1.0;	// smaller Alpha => wide distribution with divergent categories. Gain with narrower distribution.	
	_userBetaGain =2.0;		// the Alpha/Beta is the excpectation
	_userProbInvariantGain= 0.05;	// was
	_userAlphaLoss =0.5;			// loss had wider distribution (sites with no loss)	
	_userBetaLoss =0.25;			// Thus, gain:loss is 1:4
	_userProbInvariantLoss= 0.05;	// 
	_userAlphaRate =0.5;			// 
	_userBetaRate =0.5;
	_userProbInvariantRate = 0.05;	// 
	_userRateInvariantVal = 1e-6;	//
	_isFlatUserParameters = false;

// for initRand - Rand(x){min<x<max}
	_userGainMax =2.0;
	_userLossMax =5.0;
	_userThetaMax =0.9;
	_userAlphaGainMax =2.0;				//1	
	_userBetaGainMax =5.0;				//2
	_userProbInvariantGainMax= 0.1;
	_userAlphaLossMax =3.0;				//2
	_userBetaLossMax =2.0;
	_userProbInvariantLossMax= 0.1;
	_userProbInvariantRateMax = 0.1;
	_userAlphaRateMax =2.0;
	_userBetaRateMax =2.0;
	
	_userGainMin =0.1;
	_userLossMin =0.1;
	_userThetaMin =0.01;				//0.1
	_userAlphaGainMin =0.05;
	_userBetaGainMin =0.05;
	_userProbInvariantGainMin= 0.0;
	_userAlphaLossMin =0.05;
	_userBetaLossMin =0.05;
	_userProbInvariantLossMin= 0.0;
	_userProbInvariantRateMin = 0.0;
	_userAlphaRateMin =0.05;
	_userBetaRateMin =0.05;

//################################################## PostExp (Stochastic mapping based Counting)
	_numOfSimulationsForPotExp = 100000;	// the counting (expectation) is based on simulations - val: >1000 - accurate enough
	//_probCutOffSum =0.3;					// the cutOff to "ProbabilityPerPosPerBranch.txt"
	_probCutOffCounts = 0.3;				// the cutOff to estimate HGT count (0.45) "gainLossProbExpCountPerPos.txt"
	_isFewCutOffCounts = true;				// the cutOff to estimate HGT count - Few (0.1,...,0.9) "gainLossProbExpCountPerPos.txt"
	_probCutOffPrintEvent = 0.05;			// the cutOff for perPosperBranch (so that file is not too big) (0.05)

//################################################## simulate PostExp (To test to accuracy of the stochastic mapping)
	_simulationType = Gamma; // Uniform
	_isMPratio = false;
	_isInitGainLossByEmpiricalFreqSimulatePostExp = true;
	_is3states = false;
	_3statesGain = 0.66; //gain (0->1)
	_3statesMore=2.68; //more (1->more)
	_3statesLess=2.68; // less (more->1) 
	_3statesLoss=1.34; // loss (1->0)
	_3states0=0.5;
	_3states1=0.2;	//_3states2+= 1 - _3states0 + _3states1;

	_numberOfPositions2simulate =8000;	// The number of positions, seqLen, note the after Nmin filter, if there are X sites, X^2/2 pairs are computed 
	_numberOfIterations2simulate = 100; // max number of simulating iteration in parametric bootstrap, without convergence
	_numberOfIterationsForPrintResults = 5; // if =3, each 3 simulation iterations, results are updated (thus, temp results are available)

	_percentileOfNminWithCorr1RequiredForLastIteration = 10; // if 2, median Nmin wity Cor=1 is required for last simulation iteration, if 10, the ten-percentile is required for convergence

	_modelOptimizationSimPostExp =true;
	_BBLOptimizationSimPostExp =true;				// changed to tree, since the branch length are "erased"		
	_epsilonOptForPostExpSimFactor = 10;				// 1 is for normal accuracy
	_numOfIterationsOptForPostExpSimFactor = 0.1;		// 1 is for normal accuracy
	_loss2gainRatioToSim = 3;						// loss rate is 3 time that of gain
	
	_initParamsAtRandPointsInSimPostExp =true;		// these 3 options could be used as:  enum simulationType {GAMMA, UNI, MP};

	_noiseLevelInGammaSimulation =0.5;	
	_isMatrixGainLossFromRatioInSimulations =true;
	_initRootFreqAtRandPointsInSimPostExpEachPos =false;	// not required, in current settings
	_isTheataFromObservedFreq =true;						// The theta is taken from observed freq +random perturbation
	_isRootFreqEQstationaryInSimulations =true;						// 
	_isFlatSpBeforeOpt =true;								// need to change to T when performing initParamsFromTrueEstimation
	_isFlatTreeBeforOpt =true;								// In simulations - Flat the tree before Opt
	_isbBLEMwithSimpleSpSimulatePostExp =true;				// In simulations - Do BBL-EM simple

//################################################## CoEvolvingSites
	_numberOfSequences2simulate =100;
	_numberOfSequences2simulateForCoEvol = 100; // number of simulations used in the co-evoving computations
	_useTheSameSpForSim =true;
	_isReversibleSim =false;
	_rateDistributionTypeSim =GAMMA;
	_gainEQlossSim =false;
	_writeSeqSim =true;

//################################################## Misc.
	_maxRateForML =100.0;
	_minBranchLength =0.0000001;
	_maxBranchLength =10.0;
	_treeSearchAlg = njML;		// To construct tree from distanceTable (JC or others)
	_weights = NULL;			// positions are weighted (not in use)
	_isOnlyComputeLikelihood = false;
	_isSequenceUniqPattern = false;
	_isRemovePositionsWithHighPercentOfMissingData = false;
	_fractionOfMissingDataToRemove = 0.5;

	_isAnaliticComputeJumps = true;
	_isNormalizeQ = false;					// true, but it is required to change optimizeModel (such that the model is not copied, but reference is sent).
	_isNormalizeQinSpVVec = false;			// update of method is required, otherwise, global changes are made
	_isNormalizeQandTreeafterOpt = true;					// after bug fixed.
	_isAlphaEqBetaManipulation = false;		// This manipulation produces an un normalized Q matrices 
	_calculeBranchLegthDiffFactorFromInputTrees = false;	// input 2 trees - compute logL diff per branch length
	_intersectTreeAndSeq = false;	// input tree and seq (not the same taxa) - intersect, write seq and tree and return

	_outPtr =&cout;

	_costMatrixfile = "";
	_costMatrixType = gainLossCost;
	_costMatrixGainLossRatio = 2.001; // add 0.001 as tie breaker

// all the parameters are added to the static: ParamList paramList (vector<Parameter>);	
	//Parameters::addParameter("_mainType", _mainType);
	Parameters::addParameter("_alphabet_size", _alphabet_size);
	Parameters::addParameter("_treeFile", _treeFile);
	Parameters::addParameter("_treeFileOrig", _treeFileOrig);
	Parameters::addParameter("_seqFile", _seqFile);
	Parameters::addParameter("_logFile", _logFile);
	Parameters::addParameter("_numOfSimulationsForPotExp", _numOfSimulationsForPotExp);
	Parameters::addParameter("_logValue", _logValue);
	Parameters::addParameter("_referenceSeq", _referenceSeq);
	//Parameters::addParameter("_outFile", _outFile);
	//Parameters::addParameter("_outFileNotNormalize", _outFileNotNormalize);
	//Parameters::addParameter("_outFileGain4Site", _outFileGain4Site);
	//Parameters::addParameter("_outFileLoss4Site", _outFileLoss4Site);
	//Parameters::addParameter("_outFileLikeofPos", _outFileLikeofPos);
	Parameters::addParameter("_treeOutFile", _treeOutFile);
	Parameters::addParameter("_isOnlyComputeLikelihood", (_isOnlyComputeLikelihood == true) ? 1 : 0);
	Parameters::addParameter("_isSequenceUniqPattern", (_isSequenceUniqPattern == true) ? 1 : 0);
	Parameters::addParameter("_isRemovePositionsWithHighPercentOfMissingData", (_isRemovePositionsWithHighPercentOfMissingData == true) ? 1 : 0);
	Parameters::addParameter("_fractionOfMissingDataToRemove", _fractionOfMissingDataToRemove);

	Parameters::addParameter("_isAnaliticComputeJumps", (_isAnaliticComputeJumps == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeQ", (_isNormalizeQ == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeQinSpVVec", (_isNormalizeQinSpVVec == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeQandTreeafterOpt", (_isNormalizeQandTreeafterOpt == true) ? 1 : 0);
	Parameters::addParameter("_isFlatUserParameters", (_isFlatUserParameters == true) ? 1 : 0);
	Parameters::addParameter("_isAlphaEqBetaManipulation", (_isAlphaEqBetaManipulation == true) ? 1 : 0);
	Parameters::addParameter("_calculeBranchLegthDiffFactorFromInputTrees", (_calculeBranchLegthDiffFactorFromInputTrees == true) ? 1 : 0);
	Parameters::addParameter("_intersectTreeAndSeq", (_intersectTreeAndSeq == true) ? 1 : 0);
	
	//Parameters::addParameter("_discretizationType", _discretizationType);
	Parameters::addParameter("_gainDistributionType", getDistributionType(_gainDistributionType));
	Parameters::addParameter("_lossDistributionType", getDistributionType(_lossDistributionType));
	Parameters::addParameter("_rateDistributionType", getDistributionType(_rateDistributionType));	
	
	Parameters::addParameter("_userGainLossRatio", _userGainLossRatio);
	Parameters::addParameter("_keepUserGainLossRatio", _keepUserGainLossRatio);
	Parameters::addParameter("_userAlphaGain", _userAlphaGain);
	Parameters::addParameter("_userBetaGain", _userBetaGain);
	Parameters::addParameter("_userProbInvariantGain", _userProbInvariantGain);
	Parameters::addParameter("_userAlphaLoss", _userAlphaLoss);
	Parameters::addParameter("_userBetaLoss", _userBetaLoss);
	Parameters::addParameter("_userProbInvariantLoss", _userProbInvariantLoss);	
	Parameters::addParameter("_userProbInvariantRate", _userProbInvariantRate);
	Parameters::addParameter("_userRateInvariantVal", _userRateInvariantVal);
	Parameters::addParameter("_userAlphaRate", _userAlphaRate);
	Parameters::addParameter("_userBetaRate", _userBetaRate);
	Parameters::addParameter("_userGain", _userGain);
	Parameters::addParameter("_userLoss", _userLoss);
	Parameters::addParameter("_userTheta", _userTheta);

	Parameters::addParameter("_userAlphaGainMax", _userAlphaGainMax);
	Parameters::addParameter("_userBetaGainMax", _userBetaGainMax);
	Parameters::addParameter("_userProbInvariantGainMax", _userProbInvariantGainMax);
	Parameters::addParameter("_userAlphaLossMax", _userAlphaLossMax);
	Parameters::addParameter("_userBetaLossMax", _userBetaLossMax);
	Parameters::addParameter("_userProbInvariantLossMax", _userProbInvariantLossMax);	
	Parameters::addParameter("_userProbInvariantRateMax", _userProbInvariantRateMax);
	Parameters::addParameter("_userAlphaRateMax", _userAlphaRateMax);
	Parameters::addParameter("_userBetaRateMax", _userBetaRateMax);
	Parameters::addParameter("_userGainMax", _userGainMax);
	Parameters::addParameter("_userLossMax", _userLossMax);
	Parameters::addParameter("_userThetaMax", _userThetaMax);

	Parameters::addParameter("_userAlphaGainMin", _userAlphaGainMin);
	Parameters::addParameter("_userBetaGainMin", _userBetaGainMin);
	Parameters::addParameter("_userProbInvariantGainMin", _userProbInvariantGainMin);
	Parameters::addParameter("_userAlphaLossMin", _userAlphaLossMin);
	Parameters::addParameter("_userBetaLossMin", _userBetaLossMin);
	Parameters::addParameter("_userProbInvariantLossMin", _userProbInvariantLossMin);	
	Parameters::addParameter("_userProbInvariantRateMin", _userProbInvariantRateMin);
	Parameters::addParameter("_userAlphaRateMin", _userAlphaRateMin);
	Parameters::addParameter("_userBetaRateMin", _userBetaRateMin);
	Parameters::addParameter("_userGainMin", _userGainMin);
	Parameters::addParameter("_userLossMin", _userLossMin);
	Parameters::addParameter("_userThetaMin", _userThetaMin);
	Parameters::addParameter("_probCutOffPrintEvent", _probCutOffPrintEvent);
	Parameters::addParameter("_probCutOffCounts", _probCutOffCounts);
	Parameters::addParameter("_isFewCutOffCounts", _isFewCutOffCounts);

	
	Parameters::addParameter("_characterFreqEval", getCharacterFreqEvalType(_characterFreqEval));
	Parameters::addParameter("_treeSearchAlg", getTreeSearchAlgType(_treeSearchAlg));
	Parameters::addParameter("_gammmaMixtureOptimizerAlg", getGammmaMixtureOptimizerAlgType(_gammmaMixtureOptimizerAlg));
	//Parameters::addParameter("_optimizeBranchLengths", _optimizeBranchLengths);
	Parameters::addParameter("_rateEstimationMethod", getRateEstimationMethodType(_rateEstimationMethod));
	Parameters::addParameter("_rateDiscretizationType", getDiscretizationType(_rateDiscretizationType));

	Parameters::addParameter("_numberOfGainCategories", _numberOfGainCategories);
	Parameters::addParameter("_numberOfLossCategories", _numberOfLossCategories);
	Parameters::addParameter("_numberOfRateCategories", _numberOfRateCategories);
	Parameters::addParameter("_numberOfRateComponents", _numberOfRateComponents);
	
	Parameters::addParameter("_maxNumOfIterations", _maxNumOfIterations);
	Parameters::addParameter("_maxNumOfIterationsModel", _maxNumOfIterationsModel);
	Parameters::addParameter("_maxNumOfIterationsBBL", _maxNumOfIterationsBBL);
	Parameters::addParameter("_maxNumOfIterationsManyStarts", _maxNumOfIterationsManyStarts);
	Parameters::addParameter("_numberOfRandPointsInOptimization", _numberOfRandPointsInOptimization);
	Parameters::addParameter("_numberOfRandStartPoints", _numberOfRandStartPoints);

	Parameters::addParameter("_optimizationLevel", getOptimizationLevelType(_optimizationLevel));	
	Parameters::addParameter("_epsilonOptimizationIterationCycle", _epsilonOptimizationIterationCycle); 
	Parameters::addParameter("_epsilonOptimizationModel", _epsilonOptimizationModel); 
	Parameters::addParameter("_epsilonOptimizationBBL", _epsilonOptimizationBBL); 
	Parameters::addParameter("_epsilonOptimizationIterationCycleManyStarts", _epsilonOptimizationIterationCycleManyStarts);

	Parameters::addParameter("_epsilonFactor_Model", _epsilonFactor_Model);
	Parameters::addParameter("_epsilonFactor_BBL", _epsilonFactor_BBL);
	Parameters::addParameter("_numIterationsFactor_Model", _numIterationsFactor_Model);
	Parameters::addParameter("_numIterationsFactor_BBL", _numIterationsFactor_BBL);

	Parameters::addParameter("_epsilonOptForPostExpSimFactor", _epsilonOptForPostExpSimFactor);
	Parameters::addParameter("_numOfIterationsOptForPostExpSimFactor", _numOfIterationsOptForPostExpSimFactor);
	Parameters::addParameter("_loss2gainRatioToSim", _loss2gainRatioToSim);

	Parameters::addParameter("_gainLossDist", (_gainLossDist == true) ? 1 : 0);
	Parameters::addParameter("_calculateRate4site", (_calculateRate4site == true) ? 1 : 0);
	Parameters::addParameter("_calculeGainLoss4site", (_calculeGainLoss4site == true) ? 1 : 0);
	Parameters::addParameter("_printLikelihoodLandscape", (_printLikelihoodLandscape == true) ? 1 : 0);
	Parameters::addParameter("_likelihoodLandscapeIncrement", _likelihoodLandscapeIncrement);
	Parameters::addParameter("_printLikelihoodLandscapeAlphaRate", (_printLikelihoodLandscapeAlphaRate == true) ? 1 : 0);
	Parameters::addParameter("_printLikelihoodLandscapeGainLoss", (_printLikelihoodLandscapeGainLoss == true) ? 1 : 0);
	Parameters::addParameter("_printLikelihoodLandscapeTheta", (_printLikelihoodLandscapeTheta == true) ? 1 : 0);
	Parameters::addParameter("_optAlphaInIteration", (_optAlphaInIteration == true) ? 1 : 0);
	Parameters::addParameter("_optBBL_LS_InIteration", (_optBBL_LS_InIteration == true) ? 1 : 0);
	Parameters::addParameter("_optBBL_EM_InIteration", (_optBBL_EM_InIteration == true) ? 1 : 0);
	Parameters::addParameter("_printP11forgain", (_printP11forgain == true) ? 1 : 0);

	Parameters::addParameter("_printTree", (_printTree == true) ? 1 : 0);
	Parameters::addParameter("_printSeq", (_printSeq == true) ? 1 : 0);
	Parameters::addParameter("_printPij_t", (_printPij_t == true) ? 1 : 0);
	Parameters::addParameter("_printLofPos", (_printLofPos == true) ? 1 : 0);
	Parameters::addParameter("_printLofPosBothModels", (_printLofPosBothModels == true) ? 1 : 0);
	Parameters::addParameter("_performOptimizations", (_performOptimizations == true) ? 1 : 0);
	Parameters::addParameter("_correctOptimizationEpsilon", (_correctOptimizationEpsilon == true) ? 1 : 0);
	Parameters::addParameter("_performOptimizationsROOT", (_performOptimizationsROOT == true) ? 1 : 0);
	Parameters::addParameter("_performOptimizationsBBL", (_performOptimizationsBBL == true) ? 1 : 0);
	Parameters::addParameter("_performOptimizationsBBLOnlyOnce", (_performOptimizationsBBLOnlyOnce == true) ? 1 : 0);
	Parameters::addParameter("_isBblLS", (_isBblLS == true) ? 1 : 0);
	Parameters::addParameter("_isbblLSWhenbblEMdontImprove", (_isbblLSWhenbblEMdontImprove == true) ? 1 : 0);
	Parameters::addParameter("_isSkipBblEMWhenbblEMdontImprove", (_isSkipBblEMWhenbblEMdontImprove == true) ? 1 : 0);

	Parameters::addParameter("_isInitGainLossByEmpiricalFreq", (_isInitGainLossByEmpiricalFreq == true) ? 1 : 0);
	Parameters::addParameter("_isBBLEMwithSimpleSpBeforeFullOptimization", (_isBBLEMwithSimpleSpBeforeFullOptimization == true) ? 1 : 0);
	Parameters::addParameter("_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately", (_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately == true) ? 1 : 0);
	Parameters::addParameter("_isOptimizeInvariantCategoryProb", (_isOptimizeInvariantCategoryProb == true) ? 1 : 0);
	Parameters::addParameter("_isUpdateOnlyGainBetaForRatio", (_isUpdateOnlyGainBetaForRatio == true) ? 1 : 0);
	Parameters::addParameter("_isComputeLikelihoodDuringInit", (_isComputeLikelihoodDuringInit == true) ? 1 : 0);

	Parameters::addParameter("_isBblEMbeforeLSWithMissSpecifiedModel", (_isBblEMbeforeLSWithMissSpecifiedModel == true) ? 1 : 0);
	Parameters::addParameter("_isBblForceFactorCorrection", (_isBblForceFactorCorrection == true) ? 1 : 0);
	Parameters::addParameter("_BblFactorCorrection", _BblFactorCorrection);

	Parameters::addParameter("_isSkipFirstParamsOptimization", (_isSkipFirstParamsOptimization == true) ? 1 : 0);
	Parameters::addParameter("_isOptimizeParamsWithLogMinMax", (_isOptimizeParamsWithLogMinMax == true) ? 1 : 0);
	Parameters::addParameter("_isMultipleAllBranchesByFactorAtStart", (_isMultipleAllBranchesByFactorAtStart == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeAtStart", (_isNormalizeAtStart == true) ? 1 : 0);

	Parameters::addParameter("_performOptimizationsBBLManyStarts", (_performOptimizationsBBLManyStarts == true) ? 1 : 0);
	Parameters::addParameter("_simulatedAnnealing", (_simulatedAnnealing == true) ? 1 : 0);
	Parameters::addParameter("_simulatedAnnealingMinEpsilonFactor", _simulatedAnnealingMinEpsilonFactor);
	Parameters::addParameter("_simulatedAnnealingCoolingFactor", _simulatedAnnealingCoolingFactor);
	Parameters::addParameter("_performOptimizationsManyStarts", (_performOptimizationsManyStarts == true) ? 1 : 0);
	Parameters::addParameter("_gainLossDistPlusInvariant", (_gainLossDistPlusInvariant == true) ? 1 : 0);
	Parameters::addParameter("_isHGT_normal_Pij", (_isHGT_normal_Pij == true) ? 1 : 0);
	Parameters::addParameter("_isHGT_with_Q", (_isHGT_with_Q == true) ? 1 : 0);
	Parameters::addParameter("_initParamsAtRandPoints", (_initParamsAtRandPoints == true) ? 1 : 0);
	Parameters::addParameter("_initParamsAtRandPointsInOptimization", (_initParamsAtRandPointsInOptimization == true) ? 1 : 0);
	Parameters::addParameter("_calculePosteriorExpectationOfChange", (_calculePosteriorExpectationOfChange == true) ? 1 : 0);
	Parameters::addParameter("_simulatePosteriorExpectationOfChange", (_simulatePosteriorExpectationOfChange == true) ? 1 : 0);
	Parameters::addParameter("_isOnlySimulateSeq", (_isOnlySimulateSeq == true) ? 1 : 0);

	Parameters::addParameter("_modelOptimizationSimPostExp", (_modelOptimizationSimPostExp == true) ? 1 : 0);
	Parameters::addParameter("_BBLOptimizationSimPostExp", (_BBLOptimizationSimPostExp == true) ? 1 : 0);
	Parameters::addParameter("_initParamsAtRandPointsInSimPostExp", (_initParamsAtRandPointsInSimPostExp == true) ? 1 : 0);
	Parameters::addParameter("_initRootFreqAtRandPointsInSimPostExpEachPos", (_initRootFreqAtRandPointsInSimPostExpEachPos == true) ? 1 : 0);
	Parameters::addParameter("_isFlatTreeBeforOpt", (_isFlatTreeBeforOpt == true) ? 1 : 0);
	Parameters::addParameter("_isbBLEMwithSimpleSpSimulatePostExp", (_isbBLEMwithSimpleSpSimulatePostExp == true) ? 1 : 0);
	Parameters::addParameter("_noiseLevelInGammaSimulation", _noiseLevelInGammaSimulation);

	Parameters::addParameter("_isTheataFromObservedFreq", (_isTheataFromObservedFreq == true) ? 1 : 0);
	Parameters::addParameter("_isRootFreqEQstationaryInSimulations", (_isRootFreqEQstationaryInSimulations == true) ? 1 : 0);
	Parameters::addParameter("_isMatrixGainLossFromRatioInSimulations", (_isMatrixGainLossFromRatioInSimulations == true) ? 1 : 0);
	Parameters::addParameter("_isFlatSpBeforeOpt", (_isFlatSpBeforeOpt == true) ? 1 : 0);

	Parameters::addParameter("_printTreesWithProbabilityValuesAsBP", (_printTreesWithProbabilityValuesAsBP == true) ? 1 : 0);
	Parameters::addParameter("_printTreesWithExpectationValuesAsBP", (_printTreesWithExpectationValuesAsBP == true) ? 1 : 0);

	Parameters::addParameter("_printTreesWithAncestralReconstructAsBP", (_printTreesWithAncestralReconstructAsBP == true) ? 1 : 0);
	Parameters::addParameter("_printAncestralReconstructFullData", (_printAncestralReconstructFullData == true) ? 1 : 0);
	Parameters::addParameter("_printDEBUGinfo", (_printDEBUGinfo == true) ? 1 : 0);
	Parameters::addParameter("_printPropExpOfChangeFullData", (_printPropExpOfChangeFullData == true) ? 1 : 0);
	Parameters::addParameter("_printExpPerPosPerBranchMatrix", (_printExpPerPosPerBranchMatrix == true) ? 1 : 0);
	Parameters::addParameter("_printComputedCorrelations", (_printComputedCorrelations == true) ? 1 : 0);
	Parameters::addParameter("_performParametricBootstapCorrelation", (_performParametricBootstapCorrelation == true) ? 1 : 0);
	Parameters::addParameter("_usePosSpecificSimulations", (_usePosSpecificSimulations == true) ? 1 : 0);
	Parameters::addParameter("_isConsiderNegativeCorrelations", (_isConsiderNegativeCorrelations == true) ? 1 : 0);
	Parameters::addParameter("_isDivideBinsByRange", (_isDivideBinsByRange == true) ? 1 : 0);
	Parameters::addParameter("_isSortVectorOfCorrelationsBinsByLowerRateBound", (_isSortVectorOfCorrelationsBinsByLowerRateBound == true) ? 1 : 0);
	Parameters::addParameter("_isSortVectorOfCorrelationsBinsByMidRateBound", (_isSortVectorOfCorrelationsBinsByMidRateBound == true) ? 1 : 0);
	Parameters::addParameter("_relativeSizeOfOverLappedBins", _relativeSizeOfOverLappedBins);

	Parameters::addParameter("_isPrintpairWiseCorrelationsAndNmin", (_isPrintpairWiseCorrelationsAndNmin == true) ? 1 : 0);
	Parameters::addParameter("_isPrintCorrelationsOfAllPairs_Corr", (_isPrintCorrelationsOfAllPairs_Corr == true) ? 1 : 0);
	Parameters::addParameter("_isPrintCorrelationsOfAllPairs_pVal", (_isPrintCorrelationsOfAllPairs_pVal == true) ? 1 : 0);

	Parameters::addParameter("_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH", (_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH == true) ? 1 : 0);
	Parameters::addParameter("_isAllCorrTypeReqruiedToBeSignificant", (_isAllCorrTypeReqruiedToBeSignificant == true) ? 1 : 0);
	Parameters::addParameter("_isNminBasedOnCountBranchesOverCutOff", (_isNminBasedOnCountBranchesOverCutOff == true) ? 1 : 0);


	Parameters::addParameter("_numOfBinsInParametricBootstrapSimulations", _numOfBinsInParametricBootstrapSimulations);	
	Parameters::addParameter("_isAddSimulationsWithLowRate", (_isAddSimulationsWithLowRate == true) ? 1 : 0);
	Parameters::addParameter("_isFDRcorrectionForPValInCorrelation", (_isFDRcorrectionForPValInCorrelation == true) ? 1 : 0);
	Parameters::addParameter("_isComputeQVals", (_isComputeQVals == true) ? 1 : 0);
	Parameters::addParameter("_pValueCutOffForBootStrap", _pValueCutOffForBootStrap);
	Parameters::addParameter("_minExpThresholdForPValComputationForCorrelatingPair", _minExpThresholdForPValComputationForCorrelatingPair);
	Parameters::addParameter("_isUpdateMinExpThresholdGivenSimulaitonsQuantile", _isUpdateMinExpThresholdGivenSimulaitonsQuantile); // is Wrong AddParameter? Not the bool type
	Parameters::addParameter("_isUpdateMinExpThresholdGivenRealDataQuantile", _isUpdateMinExpThresholdGivenRealDataQuantile);
	Parameters::addParameter("_updateMinExpThresholdGivenRealDataQuantileVal", _updateMinExpThresholdGivenRealDataQuantileVal);
	Parameters::addParameter("_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel", _isUpdateMinExpThresholdGivenHighFractionOfHighCorrel);
	Parameters::addParameter("_isCompExtremeValDistribution", _isCompExtremeValDistribution);
	Parameters::addParameter("_minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair", _minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair);

	Parameters::addParameter("_isCorrelateWithPearson", (_isCorrelateWithPearson == true) ? 1 : 0);
	Parameters::addParameter("_isCorrelateWithSpearman", (_isCorrelateWithSpearman == true) ? 1 : 0);
	Parameters::addParameter("_isCorrelationsBasedOnMaxParsimonyMapping", (_isCorrelationsBasedOnMaxParsimonyMapping == true) ? 1 : 0);
	Parameters::addParameter("_isAlsoCorrelateWithLoss", (_isAlsoCorrelateWithLoss == true) ? 1 : 0);
	Parameters::addParameter("_isAlsoCorrelateWithBoth", (_isAlsoCorrelateWithBoth == true) ? 1 : 0);
	Parameters::addParameter("_isOnlyCorrelateWithBoth", (_isOnlyCorrelateWithBoth == true) ? 1 : 0);
	Parameters::addParameter("_isUseRateForSiteAsNminForCorrelations", (_isUseRateForSiteAsNminForCorrelations == true) ? 1 : 0);
	Parameters::addParameter("_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur", (_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur == true) ? 1 : 0);
	Parameters::addParameter("_isRemoveSimulatedPositionsBasedOnMP", (_isRemoveSimulatedPositionsBasedOnMP == true) ? 1 : 0);
	Parameters::addParameter("_minNumOfMPEvent2RemoveSimulatedPositions", _minNumOfMPEvent2RemoveSimulatedPositions);
	Parameters::addParameter("_isUpdateminNumOfMPEvent2RemoveSimulatedPositions", (_isUpdateminNumOfMPEvent2RemoveSimulatedPositions == true) ? 1 : 0);


	Parameters::addParameter("_printComputedCorrelationsAllSites", (_printComputedCorrelationsAllSites == true) ? 1 : 0);
	Parameters::addParameter("_isIgnoreCorrelationAmongSelectedSites", (_isIgnoreCorrelationAmongSelectedSites == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeForBranchExpInCorrCompute", (_isNormalizeForBranchExpInCorrCompute == true) ? 1 : 0);
	Parameters::addParameter("_isNormalizeByExpectationPerBranch", (_isNormalizeByExpectationPerBranch == true) ? 1 : 0);


	Parameters::addParameter("_selectedSitesForCorrelation", _selectedSitesForCorrelation);
	Parameters::addParameter("_calculateAncestralReconstruct", (_calculateAncestralReconstruct == true) ? 1 : 0);
	Parameters::addParameter("_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation", (_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation == true) ? 1 : 0);
	Parameters::addParameter("_checkCoEvolWithUnionPAP_against_pos", _checkCoEvolWithUnionPAP_against_pos);


	Parameters::addParameter("_isReversible", (_isReversible == true) ? 1 : 0);
	Parameters::addParameter("_isRootFreqEQstationary", (_isRootFreqEQstationary == true) ? 1 : 0);
	Parameters::addParameter("_initRandomGammaMixuteParam", (_initRandomGammaMixuteParam == true) ? 1 : 0);
	Parameters::addParameter("_incrementFactorForGain", (_incrementFactorForGain == true) ? 1 : 0);
	Parameters::addParameter("_lossBiggerGainLimit", (_lossBiggerGainLimit == true) ? 1 : 0);
	Parameters::addParameter("_slopeFactorForGain", _slopeFactorForGain);
	Parameters::addParameter("_isStartWithTheta", (_isStartWithTheta == true) ? 1 : 0);
	Parameters::addParameter("_isSkipGainOptimization", (_isSkipGainOptimization == true) ? 1 : 0);
	Parameters::addParameter("_epsilonOptimizationThetaFactor", _epsilonOptimizationThetaFactor);
	Parameters::addParameter("_isAlphaLimit", (_isAlphaLimit == true) ? 1 : 0);
	Parameters::addParameter("_isGainLimit", (_isGainLimit == true) ? 1 : 0);
	//Parameters::addParameter("_probCutOffSum", _probCutOffSum);
	Parameters::addParameter("_maxRateForML", _maxRateForML);
	Parameters::addParameter("_minBranchLength", _minBranchLength);
	Parameters::addParameter("_maxBranchLength", _maxBranchLength);
	Parameters::addParameter("_epsilonForReRootFactor", _epsilonForReRootFactor);
	Parameters::addParameter("_percentOfImprovManySarts", _percentOfImprovManySarts);
	Parameters::addParameter("_percentOfImprov", _percentOfImprov);
	Parameters::addParameter("_calculeBranchLegthDiffFactor", (_calculeBranchLegthDiffFactor == true) ? 1 : 0);
	
	Parameters::addParameter("_simulationType", getSimulationType(_simulationType));
	Parameters::addParameter("_isMPratio", (_isMPratio == true) ? 1 : 0);
	Parameters::addParameter("_isInitGainLossByEmpiricalFreqSimulatePostExp", (_isInitGainLossByEmpiricalFreqSimulatePostExp == true) ? 1 : 0);
	Parameters::addParameter("_is3states", (_is3states == true) ? 1 : 0);
	Parameters::addParameter("_3statesGain", _3statesGain);
	Parameters::addParameter("_3statesMore", _3statesMore);
	Parameters::addParameter("_3statesLess", _3statesLess);
	Parameters::addParameter("_3statesLoss", _3statesLoss);
	Parameters::addParameter("_3states0", _3states0);
	Parameters::addParameter("_3states1", _3states1);

	Parameters::addParameter("_simulateSequences", (_simulateSequences == true) ? 1 : 0);
	Parameters::addParameter("_numberOfSequences2simulate", _numberOfSequences2simulate);
	Parameters::addParameter("_numberOfPositions2simulate", _numberOfPositions2simulate);
	Parameters::addParameter("_numberOfIterations2simulate", _numberOfIterations2simulate);
	Parameters::addParameter("_numberOfIterationsForPrintResults", _numberOfIterationsForPrintResults);
	Parameters::addParameter("_percentileOfNminWithCorr1RequiredForLastIteration", _percentileOfNminWithCorr1RequiredForLastIteration);


	Parameters::addParameter("_useTheSameSpForSim", (_useTheSameSpForSim == true) ? 1 : 0);
	Parameters::addParameter("_isReversibleSim", (_isReversibleSim == true) ? 1 : 0);
	Parameters::addParameter("_rateDistributionTypeSim", getDistributionType(_rateDistributionTypeSim));
	Parameters::addParameter("_gainEQlossSim", (_gainEQlossSim == true) ? 1 : 0);
	Parameters::addParameter("_calculateRate4siteSim", (_calculateRate4siteSim == true) ? 1 : 0);
	Parameters::addParameter("_writeSeqSim", (_writeSeqSim == true) ? 1 : 0);

	Parameters::addParameter("_accountForMissingData", (_accountForMissingData == true) ? 1 : 0);
	Parameters::addParameter("_gainEQloss", (_gainEQloss == true) ? 1 : 0);
	Parameters::addParameter("_gainLossRateAreFreq", (_gainLossRateAreFreq == true) ? 1 : 0);

	Parameters::addParameter("_findCoEvolvingSitesOldNotWorking", (_findCoEvolvingSitesOldNotWorking == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_saveProbChanges_PosNodeXY", (_saveProbChanges_PosNodeXY == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_isComputeDistanceFromRootForRecent", (_isComputeDistanceFromRootForRecent == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_printAncestralReconstructPosterior", (_printAncestralReconstructPosterior == true) ? 1 : 0);
	Parameters::addParameter("_minNumOfOnes", _minNumOfOnes); // 1,3
	Parameters::addParameter("_minNumOfZeros", _minNumOfZeros); // 0,1

	Parameters::addParameter("_isOnlyParsimony", (_isOnlyParsimony == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_calculeMaxParsimonyChange", (_calculeMaxParsimonyChange == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_calculeMaxParsimonyChangeSeveralGainLossRatios", (_calculeMaxParsimonyChangeSeveralGainLossRatios == true) ? 1 : 0);// for the co evolving project
	Parameters::addParameter("_costMatrixType", getCostMatrixType(_costMatrixType));
	Parameters::addParameter("_costMatrixfile", _costMatrixfile);
	Parameters::addParameter("_costMatrixGainLossRatio", _costMatrixGainLossRatio);
}
/********************************************************************************************
getParamsFromFile
*********************************************************************************************/
void gainLossOptions::readParameters(const string& paramFileName)
{
	ifstream params(paramFileName.c_str());
	if(params.good())
		Parameters::readParameters(params);		// only place where params are read, updateParameter(paramName, param.c_str()) used
	params.close();
}
/********************************************************************************************
getParamsFromFile
*********************************************************************************************/
void gainLossOptions::getParamsFromFile(const string& paramFileName)
{	
	readParameters(paramFileName);
	readFromParameters2gainLossOptions();	
	updateDependencies();	
	readParameters(paramFileName); // if specifically asked for other value in paramFile, now without updated...
	updateParamsInRangeOverrideParamFile();
	readFromParameters2gainLossOptions(); 
}


/********************************************************************************************
Updates... Verify consistencies
*********************************************************************************************/
void gainLossOptions::readFromParameters2gainLossOptions(){
//_mainType = Parameters::getString("_mainType");
	_outDir = Parameters::getString("_outDir");
	_alphabet_size = Parameters::getInt("_alphabet_size");
	_minNumOfOnes = Parameters::getInt("_minNumOfOnes");
	_minNumOfZeros = Parameters::getInt("_minNumOfZeros");
	_numOfSimulationsForPotExp = Parameters::getInt("_numOfSimulationsForPotExp");

	_gainLossRateAreFreq = (Parameters::getInt("_gainLossRateAreFreq") == 1) ? true : false;
	_isOnlyComputeLikelihood = (Parameters::getInt("_isOnlyComputeLikelihood") == 1) ? true : false;
	_isSequenceUniqPattern = (Parameters::getInt("_isSequenceUniqPattern") == 1) ? true : false;
	_isRemovePositionsWithHighPercentOfMissingData = (Parameters::getInt("_isRemovePositionsWithHighPercentOfMissingData") == 1) ? true : false;
	_fractionOfMissingDataToRemove = Parameters::getFloat("_fractionOfMissingDataToRemove");

	_isAnaliticComputeJumps = (Parameters::getInt("_isAnaliticComputeJumps") == 1) ? true : false;
	_isNormalizeQ = (Parameters::getInt("_isNormalizeQ") == 1) ? true : false;
	_isNormalizeQinSpVVec = (Parameters::getInt("_isNormalizeQinSpVVec") == 1) ? true : false;
	_isNormalizeQandTreeafterOpt = (Parameters::getInt("_isNormalizeQandTreeafterOpt") == 1) ? true : false;
	_isFlatUserParameters = (Parameters::getInt("_isFlatUserParameters") == 1) ? true : false;
	_isAlphaEqBetaManipulation = (Parameters::getInt("_isAlphaEqBetaManipulation") == 1) ? true : false;
	_calculeBranchLegthDiffFactorFromInputTrees = (Parameters::getInt("_calculeBranchLegthDiffFactorFromInputTrees") == 1) ? true : false;
	_intersectTreeAndSeq = (Parameters::getInt("_intersectTreeAndSeq") == 1) ? true : false;

	_gainEQloss = (Parameters::getInt("_gainEQloss") == 1) ? true : false;
	_isRootFreqEQstationary = (Parameters::getInt("_isRootFreqEQstationary") == 1) ? true : false;
	_isReversible = (Parameters::getInt("_isReversible") == 1) ? true : false;
	_gainLossDist = (Parameters::getInt("_gainLossDist") == 1) ? true : false;
	
	
	_rateDistributionType = getDistributionType(Parameters::getString("_rateDistributionType"));
	if(_rateDistributionType == UNIFORM){
		_rateEstimationMethod = mlRate;
		Parameters::updateParameter("_rateEstimationMethod","mlRate");        
	}
	_gainDistributionType = getDistributionType(Parameters::getString("_gainDistributionType"));
	_lossDistributionType = getDistributionType(Parameters::getString("_lossDistributionType"));

	_lossBiggerGainLimit = (Parameters::getInt("_lossBiggerGainLimit") == 1) ? true : false;
	_userGainLossRatio = Parameters::getFloat("_userGainLossRatio");
	_keepUserGainLossRatio = (Parameters::getInt("_keepUserGainLossRatio") == 1) ? true : false;

	_userGain = Parameters::getFloat("_userGain");
	_userLoss = Parameters::getFloat("_userLoss");
	if((_lossBiggerGainLimit) && (_userLoss <= _userGain)){
		_userGain = 0.5;
		Parameters::updateParameter("_userGain","0.5");
		_userLoss = 1.5;
		Parameters::updateParameter("_userLoss","1.5");	
	}
	_performOptimizationsBBL = (Parameters::getInt("_performOptimizationsBBL") == 1) ? true : false;
	_performOptimizationsBBLOnlyOnce = (Parameters::getInt("_performOptimizationsBBLOnlyOnce") == 1) ? true : false;
	_isBblLS = (Parameters::getInt("_isBblLS") == 1) ? true : false;
	_isbblLSWhenbblEMdontImprove = (Parameters::getInt("_isbblLSWhenbblEMdontImprove") == 1) ? true : false;
	_isSkipBblEMWhenbblEMdontImprove = (Parameters::getInt("_isSkipBblEMWhenbblEMdontImprove") == 1) ? true : false;

	_isBblEMbeforeLSWithMissSpecifiedModel = (Parameters::getInt("_isBblEMbeforeLSWithMissSpecifiedModel") == 1) ? true : false;
	_isBblForceFactorCorrection = (Parameters::getInt("_isBblForceFactorCorrection") == 1) ? true : false;
	_BblFactorCorrection = Parameters::getFloat("_BblFactorCorrection");

	_isSkipFirstParamsOptimization = (Parameters::getInt("_isSkipFirstParamsOptimization") == 1) ? true : false;
	_isOptimizeParamsWithLogMinMax = (Parameters::getInt("_isOptimizeParamsWithLogMinMax") == 1) ? true : false;
	_isMultipleAllBranchesByFactorAtStart = (Parameters::getInt("_isMultipleAllBranchesByFactorAtStart") == 1) ? true : false;
	_isNormalizeAtStart = (Parameters::getInt("_isNormalizeAtStart") == 1) ? true : false;

	_performOptimizationsBBLManyStarts = (Parameters::getInt("_performOptimizationsBBLManyStarts") == 1) ? true : false;
	_simulatedAnnealing = (Parameters::getInt("_simulatedAnnealing") == 1) ? true : false;
	_simulatedAnnealingMinEpsilonFactor = Parameters::getFloat("_simulatedAnnealingMinEpsilonFactor");
	_simulatedAnnealingCoolingFactor = Parameters::getFloat("_simulatedAnnealingCoolingFactor");

	_performOptimizationsManyStarts = (Parameters::getInt("_performOptimizationsManyStarts") == 1) ? true : false;
	if(_performOptimizationsManyStarts == 1){
		_initParamsAtRandPointsInOptimization = true;
		Parameters::updateParameter("_initParamsAtRandPointsInOptimization","1");	
	}	
	_seqFile = Parameters::getString("_seqFile");
	_simulatePosteriorExpectationOfChange = (Parameters::getInt("_simulatePosteriorExpectationOfChange") == 1) ? true : false;
	_isOnlySimulateSeq = (Parameters::getInt("_isOnlySimulateSeq") == 1) ? true : false;

	if(_seqFile=="" && _simulatePosteriorExpectationOfChange==0) errorMsg::reportError("_seqFile is needed");
	_treeFile = Parameters::getString("_treeFile");
	_treeFileOrig = Parameters::getString("_treeFileOrig");

	_rootAt = Parameters::getString("_rootAt");
	_logFile= Parameters::getString("_logFile");
	_logValue = Parameters::getInt("_logValue");
	_referenceSeq = Parameters::getString("_referenceSeq");
	//_outFile = Parameters::getString("_outFile");
	_treeOutFile = Parameters::getString("_treeOutFile");
	
	//_discretizationType = Parameters::getString("_discretizationType");
	_treeSearchAlg = getTreeSearchAlgType(Parameters::getString("_treeSearchAlg"));
	_gammmaMixtureOptimizerAlg = getGammmaMixtureOptimizerAlgType(Parameters::getString("_gammmaMixtureOptimizerAlg"));
	//_optimizeBranchLengths = Parameters::getString("_optimizeBranchLengths");

	_characterFreqEval = getCharacterFreqEvalType(Parameters::getString("_characterFreqEval"));
	_rateEstimationMethod = getRateEstimationMethodType(Parameters::getString("_rateEstimationMethod"));
	_rateDiscretizationType = getDiscretizationType(Parameters::getString("_rateDiscretizationType"));

	_numberOfGainCategories = Parameters::getInt("_numberOfGainCategories");
	_numberOfLossCategories = Parameters::getInt("_numberOfLossCategories");
	_numberOfRateCategories = Parameters::getInt("_numberOfRateCategories");
	_numberOfRateComponents = Parameters::getInt("_numberOfRateComponents");
	
	_maxNumOfIterations = Parameters::getInt("_maxNumOfIterations");
	_maxNumOfIterationsModel = Parameters::getInt("_maxNumOfIterationsModel");
	_maxNumOfIterationsBBL = Parameters::getInt("_maxNumOfIterationsBBL");
	_maxNumOfIterationsManyStarts = Parameters::getInt("_maxNumOfIterationsManyStarts");
	_numberOfRandPointsInOptimization = Parameters::getInt("_numberOfRandPointsInOptimization");
	_numberOfRandStartPoints = Parameters::getInt("_numberOfRandStartPoints");
	_epsilonOptimizationModel = Parameters::getFloat("_epsilonOptimizationModel"); 
	_epsilonOptimizationBBL = Parameters::getFloat("_epsilonOptimizationBBL"); 
	_epsilonOptimizationIterationCycleManyStarts = Parameters::getFloat("_epsilonOptimizationIterationCycleManyStarts");

	_optimizationLevel = getOptimizationLevelTypeFromStr(Parameters::getString("_optimizationLevel"));
	_epsilonFactor_Model = Parameters::getFloat("_epsilonFactor_Model"); 
	_epsilonFactor_BBL = Parameters::getFloat("_epsilonFactor_BBL"); 
	_numIterationsFactor_Model = Parameters::getFloat("_numIterationsFactor_Model"); 
	_numIterationsFactor_BBL = Parameters::getFloat("_numIterationsFactor_BBL");

	_epsilonOptimizationIterationCycle = Parameters::getFloat("_epsilonOptimizationIterationCycle");
	_epsilonOptForPostExpSimFactor = Parameters::getFloat("_epsilonOptForPostExpSimFactor");	
	_numOfIterationsOptForPostExpSimFactor = Parameters::getFloat("_numOfIterationsOptForPostExpSimFactor");	
	_loss2gainRatioToSim = Parameters::getFloat("_loss2gainRatioToSim");	
	
	_userAlphaGain = Parameters::getFloat("_userAlphaGain");
	_userBetaGain = Parameters::getFloat("_userBetaGain");
	_userProbInvariantGain = Parameters::getFloat("_userProbInvariantGain");
	_userAlphaLoss = Parameters::getFloat("_userAlphaLoss");
	_userBetaLoss = Parameters::getFloat("_userBetaLoss");
	_userProbInvariantLoss = Parameters::getFloat("_userProbInvariantLoss");
	_userProbInvariantRate = Parameters::getFloat("_userProbInvariantRate");
	_userRateInvariantVal = Parameters::getFloat("_userRateInvariantVal");
	_userAlphaRate = Parameters::getFloat("_userAlphaRate");	
	_userBetaRate = Parameters::getFloat("_userBetaRate");	

	_userAlphaGainMax = Parameters::getFloat("_userAlphaGainMax");
	_userBetaGainMax = Parameters::getFloat("_userBetaGainMax");
	_userProbInvariantGainMax = Parameters::getFloat("_userProbInvariantGainMax");
	_userAlphaLossMax = Parameters::getFloat("_userAlphaLossMax");
	_userBetaLossMax = Parameters::getFloat("_userBetaLossMax");
	_userProbInvariantLossMax = Parameters::getFloat("_userProbInvariantLossMax");
	_userProbInvariantRateMax = Parameters::getFloat("_userProbInvariantRateMax");
	_userAlphaRateMax = Parameters::getFloat("_userAlphaRateMax");	
	_userBetaRateMax = Parameters::getFloat("_userBetaRateMax");
	_userGainMax = Parameters::getFloat("_userGainMax");
	_userLossMax = Parameters::getFloat("_userLossMax");
	_userThetaMax = Parameters::getFloat("_userThetaMax");

	_userAlphaGain = Parameters::getFloat("_userAlphaGain");
	_userBetaGain = Parameters::getFloat("_userBetaGain");
	_userProbInvariantGain = Parameters::getFloat("_userProbInvariantGain");
	_userAlphaLoss = Parameters::getFloat("_userAlphaLoss");
	_userBetaLoss = Parameters::getFloat("_userBetaLoss");
	_userProbInvariantLoss = Parameters::getFloat("_userProbInvariantLoss");
	_userProbInvariantRate = Parameters::getFloat("_userProbInvariantRate");
	_userAlphaRate = Parameters::getFloat("_userAlphaRate");	
	_userBetaRate = Parameters::getFloat("_userBetaRate");
	_userGain = Parameters::getFloat("_userGain");
	_userLoss = Parameters::getFloat("_userLoss");
	_userTheta = Parameters::getFloat("_userTheta");

	_probCutOffPrintEvent = Parameters::getFloat("_probCutOffPrintEvent");
	_probCutOffCounts = Parameters::getFloat("_probCutOffCounts");
	_isFewCutOffCounts = (Parameters::getInt("_isFewCutOffCounts") == 1) ? true : false;
	
	_calculateRate4site = (Parameters::getInt("_calculateRate4site") == 1) ? true : false;
	_calculeGainLoss4site = (Parameters::getInt("_calculeGainLoss4site") == 1) ? true : false;
	_printLikelihoodLandscape = (Parameters::getInt("_printLikelihoodLandscape") == 1) ? true : false;
	_likelihoodLandscapeIncrement = Parameters::getFloat("_likelihoodLandscapeIncrement");
	_printLikelihoodLandscapeAlphaRate = (Parameters::getInt("_printLikelihoodLandscapeAlphaRate") == 1) ? true : false;
	_printLikelihoodLandscapeGainLoss = (Parameters::getInt("_printLikelihoodLandscapeGainLoss") == 1) ? true : false;
	_printLikelihoodLandscapeTheta = (Parameters::getInt("_printLikelihoodLandscapeTheta") == 1) ? true : false;
	_optAlphaInIteration = (Parameters::getInt("_optAlphaInIteration") == 1) ? true : false;
	_optBBL_LS_InIteration = (Parameters::getInt("_optBBL_LS_InIteration") == 1) ? true : false;
	_optBBL_EM_InIteration = (Parameters::getInt("_optBBL_EM_InIteration") == 1) ? true : false;
	_printP11forgain = (Parameters::getInt("_printP11forgain") == 1) ? true : false;

	_printTree = (Parameters::getInt("_printTree") == 1) ? true : false;
	_printSeq = (Parameters::getInt("_printSeq") == 1) ? true : false;
	_printPij_t = (Parameters::getInt("_printPij_t") == 1) ? true : false;
	_printLofPos = (Parameters::getInt("_printLofPos") == 1) ? true : false;
	_printLofPosBothModels = (Parameters::getInt("_printLofPosBothModels") == 1) ? true : false;
	_performOptimizations = (Parameters::getInt("_performOptimizations") == 1) ? true : false;
	_correctOptimizationEpsilon = (Parameters::getInt("_correctOptimizationEpsilon") == 1) ? true : false;

	_isInitGainLossByEmpiricalFreq = (Parameters::getInt("_isInitGainLossByEmpiricalFreq") == 1) ? true : false;
	_isBBLEMwithSimpleSpBeforeFullOptimization = (Parameters::getInt("_isBBLEMwithSimpleSpBeforeFullOptimization") == 1) ? true : false;
	_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately = (Parameters::getInt("_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately") == 1) ? true : false;
	_isOptimizeInvariantCategoryProb = (Parameters::getInt("_isOptimizeInvariantCategoryProb") == 1) ? true : false;
	_isUpdateOnlyGainBetaForRatio = (Parameters::getInt("_isUpdateOnlyGainBetaForRatio") == 1) ? true : false;
	_isComputeLikelihoodDuringInit = (Parameters::getInt("_isComputeLikelihoodDuringInit") == 1) ? true : false;

	_performOptimizationsROOT = (Parameters::getInt("_performOptimizationsROOT") == 1) ? true : false;
	_initParamsAtRandPointsInOptimization = (Parameters::getInt("_initParamsAtRandPointsInOptimization") == 1) ? true : false;
	_gainLossDistPlusInvariant = (Parameters::getInt("_gainLossDistPlusInvariant") == 1) ? true : false;
	_isHGT_normal_Pij = (Parameters::getInt("_isHGT_normal_Pij") == 1) ? true : false;
	_isHGT_with_Q = (Parameters::getInt("_isHGT_with_Q") == 1) ? true : false;
	_initParamsAtRandPoints = (Parameters::getInt("_initParamsAtRandPoints") == 1) ? true : false;
	_calculePosteriorExpectationOfChange = (Parameters::getInt("_calculePosteriorExpectationOfChange") == 1) ? true : false;
	_modelOptimizationSimPostExp = (Parameters::getInt("_modelOptimizationSimPostExp") == 1) ? true : false;
	_BBLOptimizationSimPostExp = (Parameters::getInt("_BBLOptimizationSimPostExp") == 1) ? true : false;
	_initParamsAtRandPointsInSimPostExp = (Parameters::getInt("_initParamsAtRandPointsInSimPostExp") == 1) ? true : false;
	_initRootFreqAtRandPointsInSimPostExpEachPos = (Parameters::getInt("_initRootFreqAtRandPointsInSimPostExpEachPos") == 1) ? true : false;
	_isFlatTreeBeforOpt = (Parameters::getInt("_isFlatTreeBeforOpt") == 1) ? true : false;
	_isbBLEMwithSimpleSpSimulatePostExp = (Parameters::getInt("_isbBLEMwithSimpleSpSimulatePostExp") == 1) ? true : false;
	_noiseLevelInGammaSimulation = Parameters::getFloat("_noiseLevelInGammaSimulation");
	_isTheataFromObservedFreq = (Parameters::getInt("_isTheataFromObservedFreq") == 1) ? true : false;
	_isRootFreqEQstationaryInSimulations = (Parameters::getInt("_isRootFreqEQstationaryInSimulations") == 1) ? true : false;
	_isMatrixGainLossFromRatioInSimulations = (Parameters::getInt("_isMatrixGainLossFromRatioInSimulations") == 1) ? true : false;
	_isFlatSpBeforeOpt = (Parameters::getInt("_isFlatSpBeforeOpt") == 1) ? true : false;
	_printTreesWithProbabilityValuesAsBP = (Parameters::getInt("_printTreesWithProbabilityValuesAsBP") == 1) ? true : false;
	_printTreesWithExpectationValuesAsBP = (Parameters::getInt("_printTreesWithExpectationValuesAsBP") == 1) ? true : false;
	_printTreesWithAncestralReconstructAsBP = (Parameters::getInt("_printTreesWithAncestralReconstructAsBP") == 1) ? true : false;
	_printAncestralReconstructFullData = (Parameters::getInt("_printAncestralReconstructFullData") == 1) ? true : false;
	_printDEBUGinfo = (Parameters::getInt("_printDEBUGinfo") == 1) ? true : false;
	_printPropExpOfChangeFullData = (Parameters::getInt("_printPropExpOfChangeFullData") == 1) ? true : false;
	_printExpPerPosPerBranchMatrix = (Parameters::getInt("_printExpPerPosPerBranchMatrix") == 1) ? true : false;
	_printComputedCorrelations = (Parameters::getInt("_printComputedCorrelations") == 1) ? true : false;
	_performParametricBootstapCorrelation = (Parameters::getInt("_performParametricBootstapCorrelation") == 1) ? true : false;
	_usePosSpecificSimulations = (Parameters::getInt("_usePosSpecificSimulations") == 1) ? true : false;
	_isConsiderNegativeCorrelations = (Parameters::getInt("_isConsiderNegativeCorrelations") == 1) ? true : false;
	_isDivideBinsByRange = (Parameters::getInt("_isDivideBinsByRange") == 1) ? true : false;
	_isSortVectorOfCorrelationsBinsByLowerRateBound = (Parameters::getInt("_isSortVectorOfCorrelationsBinsByLowerRateBound") == 1) ? true : false;
	_isSortVectorOfCorrelationsBinsByMidRateBound = (Parameters::getInt("_isSortVectorOfCorrelationsBinsByMidRateBound") == 1) ? true : false;
	_relativeSizeOfOverLappedBins = Parameters::getFloat("_relativeSizeOfOverLappedBins");

	_isPrintpairWiseCorrelationsAndNmin = (Parameters::getInt("_isPrintpairWiseCorrelationsAndNmin") == 1) ? true : false;
	_isPrintCorrelationsOfAllPairs_Corr = (Parameters::getInt("_isPrintCorrelationsOfAllPairs_Corr") == 1) ? true : false;
	_isPrintCorrelationsOfAllPairs_pVal = (Parameters::getInt("_isPrintCorrelationsOfAllPairs_pVal") == 1) ? true : false;

	_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH = (Parameters::getInt("_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH") == 1) ? true : false;
	_isAllCorrTypeReqruiedToBeSignificant = (Parameters::getInt("_isAllCorrTypeReqruiedToBeSignificant") == 1) ? true : false;
	_isNminBasedOnCountBranchesOverCutOff = (Parameters::getInt("_isNminBasedOnCountBranchesOverCutOff") == 1) ? true : false;

	_numOfBinsInParametricBootstrapSimulations = Parameters::getInt("_numOfBinsInParametricBootstrapSimulations");
	_isAddSimulationsWithLowRate = (Parameters::getInt("_isAddSimulationsWithLowRate") == 1) ? true : false;
	_isFDRcorrectionForPValInCorrelation = (Parameters::getInt("_isFDRcorrectionForPValInCorrelation") == 1) ? true : false;
	_isComputeQVals = (Parameters::getInt("_isComputeQVals") == 1) ? true : false;
	_pValueCutOffForBootStrap = Parameters::getFloat("_pValueCutOffForBootStrap");
	_minExpThresholdForPValComputationForCorrelatingPair = Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
	_isUpdateMinExpThresholdGivenSimulaitonsQuantile = (Parameters::getInt("_isUpdateMinExpThresholdGivenSimulaitonsQuantile") == 1) ? true : false;
	_isUpdateMinExpThresholdGivenRealDataQuantile = (Parameters::getInt("_isUpdateMinExpThresholdGivenRealDataQuantile") == 1) ? true : false;
	_updateMinExpThresholdGivenRealDataQuantileVal = Parameters::getFloat("_updateMinExpThresholdGivenRealDataQuantileVal");
	_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel = (Parameters::getInt("_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel") == 1) ? true : false;
	_isCompExtremeValDistribution = (Parameters::getInt("_isCompExtremeValDistribution") == 1) ? true : false;

	_minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair = Parameters::getFloat("_minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair");

	_isCorrelateWithPearson = (Parameters::getInt("_isCorrelateWithPearson") == 1) ? true : false;
	_isCorrelateWithSpearman = (Parameters::getInt("_isCorrelateWithSpearman") == 1) ? true : false;
	_isCorrelationsBasedOnMaxParsimonyMapping = (Parameters::getInt("_isCorrelationsBasedOnMaxParsimonyMapping") == 1) ? true : false;
	_isAlsoCorrelateWithLoss = (Parameters::getInt("_isAlsoCorrelateWithLoss") == 1) ? true : false;
	_isAlsoCorrelateWithBoth = (Parameters::getInt("_isAlsoCorrelateWithBoth") == 1) ? true : false;
	_isOnlyCorrelateWithBoth = (Parameters::getInt("_isOnlyCorrelateWithBoth") == 1) ? true : false;
	_isUseRateForSiteAsNminForCorrelations = (Parameters::getInt("_isUseRateForSiteAsNminForCorrelations") == 1) ? true : false;
	_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur = (Parameters::getInt("_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur") == 1) ? true : false;
	_isRemoveSimulatedPositionsBasedOnMP = (Parameters::getInt("_isRemoveSimulatedPositionsBasedOnMP") == 1) ? true : false;
	_minNumOfMPEvent2RemoveSimulatedPositions = Parameters::getFloat("_minNumOfMPEvent2RemoveSimulatedPositions");
	_isUpdateminNumOfMPEvent2RemoveSimulatedPositions = (Parameters::getInt("_isUpdateminNumOfMPEvent2RemoveSimulatedPositions") == 1) ? true : false;

	_printComputedCorrelationsAllSites = (Parameters::getInt("_printComputedCorrelationsAllSites") == 1) ? true : false;
	_isIgnoreCorrelationAmongSelectedSites = (Parameters::getInt("_isIgnoreCorrelationAmongSelectedSites") == 1) ? true : false;
	_isNormalizeForBranchExpInCorrCompute = (Parameters::getInt("_isNormalizeForBranchExpInCorrCompute") == 1) ? true : false;
	_isNormalizeByExpectationPerBranch = (Parameters::getInt("_isNormalizeByExpectationPerBranch") == 1) ? true : false;

	_selectedSitesForCorrelation = Parameters::getString("_selectedSitesForCorrelation");
	_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation = (Parameters::getInt("_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation") == 1) ? true : false;
	_checkCoEvolWithUnionPAP_against_pos = Parameters::getInt("_checkCoEvolWithUnionPAP_against_pos");


	_calculateAncestralReconstruct = (Parameters::getInt("_calculateAncestralReconstruct") == 1) ? true : false;
	_calculeBranchLegthDiffFactor = (Parameters::getInt("_calculeBranchLegthDiffFactor") == 1) ? true : false;
	_initRandomGammaMixuteParam = (Parameters::getInt("_initRandomGammaMixuteParam") == 1) ? true : false;
	_incrementFactorForGain = (Parameters::getInt("_incrementFactorForGain") == 1) ? true : false;
	_slopeFactorForGain = Parameters::getFloat("_slopeFactorForGain");
	_isStartWithTheta = (Parameters::getInt("_isStartWithTheta") == 1) ? true : false;
	_isSkipGainOptimization = (Parameters::getInt("_isSkipGainOptimization") == 1) ? true : false;
	_epsilonOptimizationThetaFactor = Parameters::getFloat("_epsilonOptimizationThetaFactor");
	
	_isAlphaLimit = (Parameters::getInt("_isAlphaLimit") == 1) ? true : false;
	_isGainLimit = (Parameters::getInt("_isGainLimit") == 1) ? true : false;
	//_probCutOffSum = Parameters::getFloat("_probCutOffSum");
	_maxRateForML = Parameters::getFloat("_maxRateForML");
	_minBranchLength = Parameters::getFloat("_minBranchLength");
	_maxBranchLength = Parameters::getFloat("_maxBranchLength");
	_epsilonForReRootFactor = Parameters::getFloat("_epsilonForReRootFactor");
	_percentOfImprovManySarts = Parameters::getFloat("_percentOfImprovManySarts");
	_percentOfImprov = Parameters::getFloat("_percentOfImprov");
	_accountForMissingData = (Parameters::getInt("_accountForMissingData") == 1) ? true : false;

	_findCoEvolvingSitesOldNotWorking = (Parameters::getInt("_findCoEvolvingSitesOldNotWorking") == 1) ? true : false;
	_saveProbChanges_PosNodeXY = (Parameters::getInt("_saveProbChanges_PosNodeXY") == 1) ? true : false;
	_isComputeDistanceFromRootForRecent = (Parameters::getInt("_isComputeDistanceFromRootForRecent") == 1) ? true : false;
	_printAncestralReconstructPosterior = (Parameters::getInt("_printAncestralReconstructPosterior") == 1) ? true : false;
	_numberOfSequences2simulateForCoEvol = (Parameters::getInt("_numberOfSequences2simulateForCoEvol"));

	_simulationType = getSimulationTypeFromStr(Parameters::getString("_simulationType"));
	_isMPratio = (Parameters::getInt("_isMPratio") == 1) ? true : false;
	_isInitGainLossByEmpiricalFreqSimulatePostExp = (Parameters::getInt("_isInitGainLossByEmpiricalFreqSimulatePostExp") == 1) ? true : false;
	_is3states = (Parameters::getInt("_is3states") == 1) ? true : false;
	_3statesGain = Parameters::getFloat("_3statesGain");
	_3statesMore = Parameters::getFloat("_3statesMore");
	_3statesLess = Parameters::getFloat("_3statesLess");
	_3statesLoss = Parameters::getFloat("_3statesLoss");
	_3states0 = Parameters::getFloat("_3states0");
	_3states1 = Parameters::getFloat("_3states1");

	_simulateSequences = (Parameters::getInt("_simulateSequences") == 1) ? true : false;
	_useTheSameSpForSim = (Parameters::getInt("_useTheSameSpForSim") == 1) ? true : false;
	_isReversibleSim = (Parameters::getInt("_isReversibleSim") == 1) ? true : false;
	_numberOfSequences2simulate = Parameters::getInt("_numberOfSequences2simulate");
	_numberOfPositions2simulate = Parameters::getInt("_numberOfPositions2simulate");
	_numberOfIterations2simulate = Parameters::getInt("_numberOfIterations2simulate");
	_numberOfIterationsForPrintResults = Parameters::getInt("_numberOfIterationsForPrintResults");
	_percentileOfNminWithCorr1RequiredForLastIteration = Parameters::getFloat("_percentileOfNminWithCorr1RequiredForLastIteration");

	_rateDistributionTypeSim = getDistributionType(Parameters::getString("_rateDistributionTypeSim"));
	_gainEQlossSim = (Parameters::getInt("_gainEQlossSim") == 1) ? true : false;
	_calculateRate4siteSim = (Parameters::getInt("_calculateRate4siteSim") == 1) ? true : false;

	_isOnlyParsimony = (Parameters::getInt("_isOnlyParsimony") == 1) ? true : false;
	_calculeMaxParsimonyChange = (Parameters::getInt("_calculeMaxParsimonyChange") == 1) ? true : false;
	_calculeMaxParsimonyChangeSeveralGainLossRatios = (Parameters::getInt("_calculeMaxParsimonyChangeSeveralGainLossRatios") == 1) ? true : false;
	_costMatrixType = getCostMatrixTypeFromStr(Parameters::getString("_costMatrixType"));
	_costMatrixfile = Parameters::getString("_costMatrixfile");
	_costMatrixGainLossRatio = Parameters::getFloat("_costMatrixGainLossRatio");	
	if(_calculateRate4siteSim || _findCoEvolvingSitesOldNotWorking){
		_writeSeqSim = true;
		Parameters::updateParameter("_writeSeqSim","1");			
	}
	_writeSeqSim = (Parameters::getInt("_writeSeqSim") == 1) ? true : false;	

	if(_rateDistributionType == GAMMA_MIXTURE){	// TEMP - not DEBBUGED
		if(_performOptimizationsManyStarts){
			cout<<"For GAMMA_MIXTURE - OptimizationsManyStarts is not fully functional.";
	//		_performOptimizationsManyStarts =0;
	//		Parameters::updateParameter("_performOptimizationsManyStarts","0");
		}	
	}
}



/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::getOutDirFromFile(const string& paramFileName)
{
	_outDir = "RESULTS";
	Parameters::addParameter("_outDir", _outDir);

	ifstream params(paramFileName.c_str());
	if(params.good())
		Parameters::readParameters(params);
	params.close();
	_outDir = Parameters::getString("_outDir");
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::verifyConsistParams()
{

	if((_isReversible || gainLossOptions::_isRootFreqEQstationary) && // fixedRoot
		(gainLossOptions::_isBblEMbeforeLSWithMissSpecifiedModel || !gainLossOptions::_isBblLS) && // BBL-EM
		(gainLossOptions::_gainDistributionType==gainLossOptions::UNIFORM) // UNIFORM
		)
		errorMsg::reportError("BBL-EM fixedRoot is not working with UNIFORM");

	if(gainLossOptions::_isAlsoCorrelateWithLoss)
		LOGnOUT(3,<<"WARN: compute correlatins for co-Loss, printComputedCorrelationsData() is problematic (not all pair will have both co-gain and co-loss defined) \n");

	if(_gainLossDist == true && !(_rateDistributionType == UNIFORM)){
		cout<<"WARNING:!!! In params: _gainLossDist == 1 but _rateDistributionType != UNIFORM (update to UNIFORM)\n";
		_rateDistributionType = UNIFORM;
		Parameters::updateParameter("_rateDistributionType","UNIFORM");		
	}
	//if(gainLossOptions::_isReversible && gainLossOptions::_calculePosteriorExpectationOfChange)
	//	errorMsg::reportError("calculePosteriorExpectationOfChange is not implemented for Reversible process");
	//if((gainLossOptions::_rateDistributionType == UNIFORM) && gainLossOptions::_calculePosteriorExpectationOfChange)
	//	errorMsg::reportError("calculePosteriorExpectationOfChange is not implemented for UNIFORM rate");
	//if(gainLossOptions::_gainLossDist && gainLossOptions::_printLikelihoodLandscape)
	//	errorMsg::reportError("LikelihoodLandscape is not implemented for spVVec(gainLossDist)");
	//if(gainLossOptions::_gainLossDist && gainLossOptions::_performOptimizationsBBL)
	//	errorMsg::reportError("BBL is not implemented for spVVec(gainLossDist)");
	//if(gainLossOptions::_accountForMissingData && (gainLossOptions::_rateDistributionType == GAMMA_MIXTURE))
	//	errorMsg::reportError("accountForMissingData is not implemented with GAMMA_MIXTURE");
}

/********************************************************************************************
Updates... Verify consistencies
*********************************************************************************************/
void gainLossOptions::updateDependencies(){
	if(_simulatedAnnealing){
		cout<<"In params: _simulatedAnnealing -> double the normal epsilons\n";
		updateOptimizationLevel(low);
	}
	updateGainLossDist();
	updateAccountForMissingData();
	updateRemoveComputationNotSuiteForModels();
	updateSimulatePosteriorExpectationOfChange();
	updateInitParamsAtRandPointsInSimPostExp();
	updateGainEQloss();	
	updateGainLossAsFreq();
	updateUserGainLossRatio(_userGainLossRatio);
	updateKeepUserGainLossRatio();
	updateOnlyComputeLikelihood();
	updateFlatUserParameters();
	updateNoBBL();
	updateNoOptimization();
	updateNoBranchLengthDiffComputation();
	updateOptimizationLevel(_optimizationLevel); // should be after updateNoBBL
	updatNoSeq();
	updateParsimonyRun();
	if(_performParametricBootstapCorrelation)
		updatParametericBootstrapComputationOfCorrelation();
}


/********************************************************************************************
Updates... Verify consistencies
*********************************************************************************************/
void gainLossOptions::updateOptimizationLevel(optimizationLevel level)
{
	MDOUBLE epsilonFactor = 1;

	if(level == mid)
		return; // no change
	switch (level)	//	enum optimizationLevel {VVVlow,VVlow, Vlow, low, mid, high, Vhigh};
	{
		case VVVlow:
			epsilonFactor = 10;
			_maxNumOfIterations = 1;				
			_maxNumOfIterationsModel = 1;
			_maxNumOfIterationsBBL = 1;
			_numberOfRandPointsInOptimization = 1;
			_numberOfRandStartPoints = 10;
			_percentOfImprov = 0.0002;
			_correctOptimizationEpsilon = 1;
			_isOptimizeInvariantCategoryProb = false;
			break;
		case VVlow:
			epsilonFactor = 8;
			_maxNumOfIterations = 1;				
			_maxNumOfIterationsModel = 1;
			_maxNumOfIterationsBBL = 1;
			_numberOfRandPointsInOptimization = 2;
			_numberOfRandStartPoints = 20;
			_percentOfImprov = 0.0001;
			_correctOptimizationEpsilon = 1;
			_isOptimizeInvariantCategoryProb = false;
			break;
		case Vlow:
			epsilonFactor = 5;
			_maxNumOfIterations = 1;				
			_maxNumOfIterationsModel = 1;			  
			_maxNumOfIterationsBBL = 1;
			_numberOfRandPointsInOptimization = 3;
			_numberOfRandStartPoints = 30;
			_percentOfImprov = 0.00002;
			_correctOptimizationEpsilon = 1;
			_isOptimizeInvariantCategoryProb = false;
			break;
		case low: // same as Vlow
			epsilonFactor = 5;
			_maxNumOfIterations = 1;				
			_maxNumOfIterationsModel = 1;			  
			_maxNumOfIterationsBBL = 1;
			_numberOfRandPointsInOptimization = 3;
			_numberOfRandStartPoints = 30;
			_percentOfImprov = 0.00002;
			_correctOptimizationEpsilon = 1;
			_isOptimizeInvariantCategoryProb = false;
			break;
		case mid:			
			break;
		case high:
			epsilonFactor = 0.5;
			break;
		case Vhigh:
			epsilonFactor = 0.1;
			//_isBblLS = true;
			//Parameters::updateParameter("_isBblLS","0");
			_isbblLSWhenbblEMdontImprove = true;
			Parameters::updateParameter("_isbblLSWhenbblEMdontImprove","1");

			break;
	}
	cout<<"In params: updateOptimizationLevel -> multiply the normal epsilons by "<<epsilonFactor<<"\n";

	_epsilonOptimizationIterationCycle *=epsilonFactor;		
	Parameters::updateParameter("_epsilonOptimizationIterationCycle",double2string(_epsilonOptimizationIterationCycle).c_str());
	_epsilonOptimizationModel *=epsilonFactor;		
	Parameters::updateParameter("_epsilonOptimizationModel",double2string(_epsilonOptimizationModel).c_str());
	_epsilonOptimizationBBL *=epsilonFactor;			
	Parameters::updateParameter("_epsilonOptimizationBBL",double2string(_epsilonOptimizationBBL).c_str());
	
	Parameters::updateParameter("_maxNumOfIterations",double2string(_maxNumOfIterations).c_str());
	Parameters::updateParameter("_maxNumOfIterationsModel",double2string(_maxNumOfIterationsModel).c_str());
	Parameters::updateParameter("_maxNumOfIterationsBBL",double2string(_maxNumOfIterationsBBL).c_str());
	
	Parameters::updateParameter("_numberOfRandPointsInOptimization",double2string(_numberOfRandPointsInOptimization).c_str());
	Parameters::updateParameter("_numberOfRandStartPoints",double2string(_numberOfRandStartPoints).c_str());

	// lowering the epsilon seems problematic - alternative?
	Parameters::updateParameter("_percentOfImprov",double2string(_percentOfImprov).c_str());
	Parameters::updateParameter("_correctOptimizationEpsilon",int2string(_correctOptimizationEpsilon).c_str());
	Parameters::updateParameter("_isOptimizeInvariantCategoryProb",int2string(_isOptimizeInvariantCategoryProb).c_str());

	if(level < 4){ // 4 is mid, 
		_isBblLS = false;
		Parameters::updateParameter("_isBblLS","0");
		_isbblLSWhenbblEMdontImprove = false;
		Parameters::updateParameter("_isbblLSWhenbblEMdontImprove","0");
		if(level < 2){ // VVVlow and VVlow - no BBL
			_performOptimizationsBBL = false;
			Parameters::updateParameter("_performOptimizationsBBL","0");
			_performOptimizationsManyStarts =false;
			Parameters::updateParameter("_performOptimizationsManyStarts","0");
			//if(level > 0){	// other than VVVlow - prior simple BBLEM performed
			//	_isBBLEMwithSimpleSpBeforeFullOptimization = true;
			//	Parameters::updateParameter("_isBBLEMwithSimpleSpBeforeFullOptimization","1");
			//}
		}		
	}

}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateUserGainLossRatio(MDOUBLE _userGainLossRatio)
{	
	if(!(_userGainLossRatio<VERYBIG)) // then it was not given
		return;	
	cout<<"In params: _userGainLossRatio -> Change gain, loss, Beta, Theta to adapt by" <<_userGainLossRatio<<"\n";
	MDOUBLE basicRate = 1; // there is no need for this parameter...
	_userGain = basicRate*sqrt(_userGainLossRatio);
	
	Parameters::updateParameter("_userGain",double2string(_userGain).c_str());
	if(_userGainLossRatio == 0)
		_userLoss =1;
	else
        _userLoss = basicRate*sqrt(1/_userGainLossRatio);	
	Parameters::updateParameter("_userLoss",double2string(_userLoss).c_str());

	//MDOUBLE computedTheta = 0.5/(_userGainLossRatio/0.1);
	MDOUBLE computedTheta = _userGain/(_userGain+_userLoss);

	if(computedTheta<1 && computedTheta>0)
		_userTheta = computedTheta;	// in case _userGainLossRatio is smaller then 0.05
	else
		_userTheta = _userThetaMax;
	//_userTheta = _userGainLossRatio/(1+_userGainLossRatio);	// ???	
	Parameters::updateParameter("_userTheta",double2string(_userTheta).c_str());
	
	//_isStartWithTheta = true; // why is it required?
	//Parameters::updateParameter("_isStartWithTheta","1");

	if(_gainLossDist == 1 && (_userGainLossRatio<pow(10.0,-10.0) || _userGainLossRatio>pow(10.0,10.0)))
		LOGnOUT(3,<<"WARN: with Mixture model, no extreme gain/loss ratios are possible\n");
	MDOUBLE gainLossRatioToCompleteByBeta = _userGainLossRatio*(_userAlphaLoss/_userAlphaGain);
	if(_userGainLossRatio == 0)
		_userBetaGain = VERYBIG;
	else
		if(_isUpdateOnlyGainBetaForRatio)
			_userBetaGain =_userBetaLoss/gainLossRatioToCompleteByBeta;			// AlphaGain = 0.35
		else
			_userBetaGain =sqrt(1/gainLossRatioToCompleteByBeta);			// AlphaGain = 0.35		
	Parameters::updateParameter("_userBetaGain",double2string(_userBetaGain).c_str());

	if(!_isUpdateOnlyGainBetaForRatio){
		if(_userGainLossRatio == 0)
			_userBetaGain = VERYSMALL;
		else
			_userBetaLoss =sqrt(gainLossRatioToCompleteByBeta);				// AlphaLoss = 0.9		
		Parameters::updateParameter("_userBetaLoss",double2string(_userBetaLoss).c_str());
	}
	_isInitGainLossByEmpiricalFreq = false;
	Parameters::updateParameter("_isInitGainLossByEmpiricalFreq","0");	
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateGainLossAsFreq()
{	
	if(!_gainLossRateAreFreq)
		return;
	cout<<"In params: _gainLossRateAreFreq -> adapt g+l=1, max val = 1\n";
	_userGain= 0.4;
	Parameters::updateParameter("_userGain","0.4");
	_userLoss = 0.6;
	Parameters::updateParameter("_userLoss","0.6");
	_userGainMax = 0.9999;
	Parameters::updateParameter("_userGainMax","0.9999");
	_userLossMax = 0.9999;
	Parameters::updateParameter("_userLossMax","0.9999");
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateParamsInRangeOverrideParamFile()
{	
	_userTheta = max(_userTheta,1e-06);
	_userTheta = min(_userTheta,1-1e-06);
	Parameters::updateParameter("_userTheta",double2string(_userTheta).c_str());
	//_userGain = max(_userGain,1e-06);
	//Parameters::updateParameter("_userGain",double2string(_userGain).c_str());


}



/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updatNoSeq()
{	
	if(_seqFile!="")
		return;
	cout<<"In params: no Seq file -> \n";
	_isTheataFromObservedFreq= false;
	Parameters::updateParameter("_isTheataFromObservedFreq","0");
	_characterFreqEval = FiftyFifty;
	Parameters::updateParameter("_characterFreqEval","FiftyFifty");
	if(_simulationType == MPestEmp || _simulationType == SMestEmp)
		errorMsg::reportError("The simulation scenario based on real data, in _simulationType=MPestEmp or SMestEmp requires input Seq.\n");

}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updatParametericBootstrapComputationOfCorrelation()
{	
	cout<<"In params: ParametericBootstrapComputationOfCorrelation -> \n";

	_calculePosteriorExpectationOfChange= true;
	Parameters::updateParameter("_calculePosteriorExpectationOfChange","1");
	_printComputedCorrelations = true;
	Parameters::updateParameter("_printComputedCorrelations","1");
	if(gainLossOptions::_selectedSitesForCorrelation==""){
		_printComputedCorrelationsAllSites = true;
		Parameters::updateParameter("_printComputedCorrelationsAllSites","1");
	}
	_calculateRate4site = false;
	Parameters::updateParameter("_calculateRate4site","0");
	_calculeGainLoss4site = false;
	Parameters::updateParameter("_calculeGainLoss4site","0");
	_calculeMaxParsimonyChange = false;
	Parameters::updateParameter("_calculeMaxParsimonyChange","0");
	_calculateAncestralReconstruct = false;
	Parameters::updateParameter("_calculateAncestralReconstruct","0");
	_printLofPos = false;
	Parameters::updateParameter("_printLofPos","0");
	_isNormalizeQandTreeafterOpt = true;								// with NoOpt - false is the default
	Parameters::updateParameter("_isNormalizeQandTreeafterOpt","1");
	//_performOptimizationsBBL = false;
	//Parameters::updateParameter("_performOptimizationsBBL","0");
	_calculeBranchLegthDiffFactor = false;
	Parameters::updateParameter("_calculeBranchLegthDiffFactor","0");

	if(_usePosSpecificSimulations){
		_isOnlySimulateSeq = true;
		Parameters::updateParameter("_isOnlySimulateSeq","1");
		_simulationType = Gamma;
		Parameters::updateParameter("_simulationType", "Gamma");
		_numberOfSequences2simulate = 1;
		Parameters::updateParameter("_numberOfSequences2simulate", "1");
	}
	if(_isSortVectorOfCorrelationsBinsByLowerRateBound){
		_isSortVectorOfCorrelationsBinsByMidRateBound = false;
		Parameters::updateParameter("_isSortVectorOfCorrelationsBinsByMidRateBound","0");
		_numOfBinsInParametricBootstrapSimulations = 20;
		Parameters::updateParameter("_numOfBinsInParametricBootstrapSimulations","10");
	}
	if(_isSortVectorOfCorrelationsBinsByMidRateBound){
		_isSortVectorOfCorrelationsBinsByLowerRateBound = false;
		Parameters::updateParameter("_isSortVectorOfCorrelationsBinsByLowerRateBound","0");
		_numOfBinsInParametricBootstrapSimulations = 10;
		Parameters::updateParameter("_numOfBinsInParametricBootstrapSimulations","10");
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateNoBranchLengthDiffComputation()
{
	if(_performOptimizationsBBL == 0 || _performOptimizations == 0){
		cout<<"In params: _performOptimizationsBBL =false -> _calculeBranchLegthDiffFactor =false\n";
		_calculeBranchLegthDiffFactor = false;
		Parameters::updateParameter("_calculeBranchLegthDiffFactor","0");	
	}
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateNoBBL()
{	
	if(_performOptimizationsBBL)
		return;
	cout<<"In params: _performOptimizationsBBL =false -> _isBBLEMwithSimpleSpBeforeFullOptimization =false\n";
	_isBBLEMwithSimpleSpBeforeFullOptimization = false;
	Parameters::updateParameter("_isBBLEMwithSimpleSpBeforeFullOptimization","0");
}


/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateGainEQloss()
{	
	if(!_gainEQloss)
		return;
	cout<<"In params: _gainEQloss -> FiftyFifty, and Reversible\n";
	_characterFreqEval = FiftyFifty;
	Parameters::updateParameter("_characterFreqEval","FiftyFifty");	
	_isReversible = true;
	Parameters::updateParameter("_isReversible","1");	
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateKeepUserGainLossRatio()
{	
	if(!_keepUserGainLossRatio)
		return;
	cout<<"In params: _keepUserGainLossRatio -> No _isInitGainLossByEmpiricalFreq\n";
	_isInitGainLossByEmpiricalFreq = false;
	Parameters::updateParameter("_isInitGainLossByEmpiricalFreq","0");	
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateRemoveComputationNotSuiteForModels()
{	
	if(_isReversible){
		cout<<"In params: _isReversible -> _calculePosteriorExpectationOfChange = false\n";
		_calculePosteriorExpectationOfChange = false;
		Parameters::updateParameter("_calculePosteriorExpectationOfChange","0");
	}
	if(_rateDistributionType == UNIFORM && !_gainLossDist){	// TEMP - not DEBBUGED
		//cout<<"In params: rateDistributionType == UNIFORM -> _calculePosteriorExpectationOfChange and _calculateAncestralReconstruct = false\n";
		//_calculePosteriorExpectationOfChange = false;
		//Parameters::updateParameter("_calculePosteriorExpectationOfChange","0");
		//_calculateAncestralReconstruct = false;
		//Parameters::updateParameter("_calculateAncestralReconstruct","0");		
	}

}
/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateGainLossDist()
{	
	if(_gainLossDist){
		cout<<"In params: _gainLossDist == 1 -> _rateDistributionType = UNIFORM (prevent to option for inner complex stochastic process)\n";
		_rateDistributionType = UNIFORM;
		Parameters::updateParameter("_rateDistributionType","UNIFORM");		
		_calculateRate4site = false;
		Parameters::updateParameter("_calculateRate4site","0");
		//_isBblLS = true;
		//Parameters::updateParameter("_isBblLS","1");
		
		if((_gainDistributionType == GENERAL_GAMMA_PLUS_INV) || (_lossDistributionType == GAMMA_PLUS_INV)){
			_gainLossDistPlusInvariant = true;
			Parameters::updateParameter("_gainLossDistPlusInvariant","1"); 
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateAccountForMissingData()
{	
	if(!_accountForMissingData){
		_minNumOfOnes =0;
		_minNumOfZeros =0;
		Parameters::updateParameter("_minNumOfOnes","0");
		Parameters::updateParameter("_minNumOfZeros","0");
	}
	if(_accountForMissingData && _minNumOfOnes ==0 && _minNumOfZeros ==0){
		_accountForMissingData =false;
		Parameters::updateParameter("_accountForMissingData","0");
	}
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateInitParamsAtRandPointsInSimPostExp()
{	
	//if(_initParamsFromMPEstimation || _initParamsFromMPratio || _initParamsFromTrueEstimation || _initParamsFromGammaWithNoise){
	//	_initParamsAtRandPointsInSimPostExp = false;
	//	Parameters::updateParameter("_initParamsAtRandPointsInSimPostExp","0");
	//}
	
	if(_simulationType == Gamma){
		_isFlatSpBeforeOpt = true;
		Parameters::updateParameter("_isFlatSpBeforeOpt","1");
	}
	if(_simulationType == GammaNoise){
		_modelOptimizationSimPostExp = false;
		Parameters::updateParameter("_modelOptimizationSimPostExp","0");
	}
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateSimulatePosteriorExpectationOfChange()
{	
	if(!_simulatePosteriorExpectationOfChange)
		return;
	cout<<"In params: _simulatePosteriorExpectationOfChange -> no Opt, no Calculations ...\n";
	_performOptimizations = false;
	Parameters::updateParameter("_performOptimizations","0");
	_calculateAncestralReconstruct = false;
	Parameters::updateParameter("_calculateAncestralReconstruct","0");
	_calculateRate4site = false;
	Parameters::updateParameter("_calculateRate4site","0");
	_calculeGainLoss4site = false;
	Parameters::updateParameter("_calculeGainLoss4site","0");
	//_calculePosteriorExpectationOfChange = false;
	//Parameters::updateParameter("_calculePosteriorExpectationOfChange","0");	// required for SMestEmp
	_printLofPos = false;
	Parameters::updateParameter("_printLofPos","0");
	//_printSeq = false;
	//Parameters::updateParameter("_printSeq","0");
	_printTree = false;
	Parameters::updateParameter("_printTree","0");
	_lossBiggerGainLimit = true;
	Parameters::updateParameter("_lossBiggerGainLimit","1");
	_printPropExpOfChangeFullData = 1;
	Parameters::updateParameter("_printPropExpOfChangeFullData","1");
	_probCutOffPrintEvent = 0;
	Parameters::updateParameter("_probCutOffPrintEvent","0");
	_calculeMaxParsimonyChangeSeveralGainLossRatios =1;
	Parameters::updateParameter("_calculeMaxParsimonyChangeSeveralGainLossRatios","1");

	//_isRootFreqEQstationary =1;
	//Parameters::updateParameter("_isRootFreqEQstationary","1");
	
	_isbblLSWhenbblEMdontImprove =0;
	Parameters::updateParameter("_isbblLSWhenbblEMdontImprove","0");

	if(_seqFile==""){
		_isInitGainLossByEmpiricalFreq = 0;
		Parameters::updateParameter("_isInitGainLossByEmpiricalFreq","0");
	}	
	// Note: if tree is not Flatned (branches) there is no need for skip
	//_isSkipFirstParamsOptimization =1;
	//Parameters::updateParameter("_isSkipFirstParamsOptimization","1");
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateOnlyComputeLikelihood()
{	
	if(!_isOnlyComputeLikelihood)
		return;
	cout<<"In params: _isOnlyComputeLikelihood -> only Opt, no Calculations ...\n";
	_calculateRate4site = false;
	Parameters::updateParameter("_calculateRate4site","0");
	_calculeGainLoss4site = false;
	Parameters::updateParameter("_calculeGainLoss4site","0");
	_calculePosteriorExpectationOfChange = false;
	Parameters::updateParameter("_calculePosteriorExpectationOfChange","0");
	_calculeMaxParsimonyChange = false;
	Parameters::updateParameter("_calculeMaxParsimonyChange","0");
	_calculateAncestralReconstruct = false;
	Parameters::updateParameter("_calculateAncestralReconstruct","0");
	_calculeBranchLegthDiffFactor =false;
	Parameters::updateParameter("_calculeBranchLegthDiffFactor","0");
	_printSeq = false;
	Parameters::updateParameter("_printSeq","0");

	_printLofPos = true;
	Parameters::updateParameter("_printLofPos","1");
}
/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateNoOptimization()
{	
	if(_performOptimizations)
		return;
	cout<<"In params: _performOptimizations = F -> no Opt\n";
	_isMultipleAllBranchesByFactorAtStart = false;
	Parameters::updateParameter("_isMultipleAllBranchesByFactorAtStart","0");
	_isBBLEMwithSimpleSpBeforeFullOptimization = false;
	Parameters::updateParameter("_isBBLEMwithSimpleSpBeforeFullOptimization","0");
	_isNormalizeAtStart = false;
	Parameters::updateParameter("_isNormalizeAtStart","0");
	_isAlphaEqBetaManipulation = false;
	Parameters::updateParameter("_isAlphaEqBetaManipulation","0");
	_isNormalizeQandTreeafterOpt = false;
	Parameters::updateParameter("_isNormalizeQandTreeafterOpt","0");
	_isInitGainLossByEmpiricalFreq = false;
	Parameters::updateParameter("_isInitGainLossByEmpiricalFreq","0");
}

/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateParsimonyRun()
{	
	if(!_isCorrelationsBasedOnMaxParsimonyMapping && !_isOnlyParsimony)
		return;
	
	cout<<"In params: _performOptimizations = F -> no Opt\n";

	_performOptimizations = false;
	Parameters::updateParameter("_performOptimizations","0");
	_isMultipleAllBranchesByFactorAtStart = false;
	Parameters::updateParameter("_isMultipleAllBranchesByFactorAtStart","0");
	_isBBLEMwithSimpleSpBeforeFullOptimization = false;
	Parameters::updateParameter("_isBBLEMwithSimpleSpBeforeFullOptimization","0");
	_isNormalizeAtStart = false;
	Parameters::updateParameter("_isNormalizeAtStart","0");
	_isAlphaEqBetaManipulation = false;
	Parameters::updateParameter("_isAlphaEqBetaManipulation","0");
	_isNormalizeQandTreeafterOpt = false;
	Parameters::updateParameter("_isNormalizeQandTreeafterOpt","0");
	_isInitGainLossByEmpiricalFreq = false;
	Parameters::updateParameter("_isInitGainLossByEmpiricalFreq","0");
	_isComputeLikelihoodDuringInit = false;
	Parameters::updateParameter("_isComputeLikelihoodDuringInit","0");

}



/********************************************************************************************
*********************************************************************************************/
void gainLossOptions::updateFlatUserParameters()
{	
	if(!_isFlatUserParameters)
		return;
	cout<<"In params: _isFlatUserParameters -> all user paramas are 1.\n";
	_userGain = 1.0;
	Parameters::updateParameter("_userGain","1.0");
	_userLoss = 1.0;
	Parameters::updateParameter("_userLoss","1.0");
	_userTheta =0.5;
	Parameters::updateParameter("_userTheta","0.5");
	_userAlphaGain =1.0;
	Parameters::updateParameter("_userAlphaGain","1.0");		
	_userBetaGain =1.0;
	Parameters::updateParameter("_userBetaGain","1.0");		
	_userAlphaLoss =1.0;
	Parameters::updateParameter("_userAlphaLoss","1.0");
	_userBetaLoss =1.0;
	Parameters::updateParameter("_userBetaLoss","1.0");
	_userAlphaRate =1.0;
	Parameters::updateParameter("_userAlphaRate","1.0");
	_userBetaRate =1.0;
	Parameters::updateParameter("_userBetaRate","1.0");
}






/********************************************************************************************
Types
enum optimizationLevel {Vlow, low, mid, high, Vhigh};
*********************************************************************************************/
string gainLossOptions::getOptimizationLevelType(optimizationLevel type)
{
	string res = "";
	switch (type) //{VVlow, Vlow, low, mid, high, Vhigh}
	{
	case VVVlow:
		res = "VVVlow";
		break;	
	case VVlow:
		res = "VVlow";
		break;	
	case Vlow:
		res = "Vlow";
		break;
	case low:
		res = "low";
		break;
	case mid:
		res = "mid";
		break;
	case high:
		res = "high";
		break;
	case Vhigh:
		res = "Vhigh";
		break;
	default:
		errorMsg::reportError("unknown type in optimizationLevel - {VVVlow,VVlow,Vlow, low, mid, high, Vhigh}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::optimizationLevel gainLossOptions::getOptimizationLevelTypeFromStr(const string& str)
{
	optimizationLevel returnType;
	if (str == "VVVlow")
		returnType = VVVlow;
	else if (str == "VVlow")
		returnType = VVlow;
	else if (str == "Vlow")
		returnType = Vlow;
	else if (str=="low")
		returnType = low;
	else if (str=="mid")
		returnType = mid;
	else if (str=="high")
		returnType = high;
	else if (str=="Vhigh")
		returnType = Vhigh;
	else
		errorMsg::reportError("unknown type in gainLossOptions::optimizationLevel- {VVVlow,VVlow,Vlow, low, mid, high, Vhigh}");
	return returnType;
}

/********************************************************************************************
enum costMatrixType {file,fitch,diff,diffSquare,gainLossCost};
*********************************************************************************************/
string gainLossOptions::getCostMatrixType(costMatrixType type)
{
	string res = "";
	switch (type)
	{
	case file:
		res = "file";
		break;
	case fitch:
		res = "fitch";
		break;
	case diff:
		res = "diff";
		break;
	case diffSquare:
		res = "diffSquare";
		break;
	case gainLossCost:
		res = "gainLossCost";
		break;
	default:
		errorMsg::reportError("unknown type in gainLossOptions::getCostMatrixType - {file,fitch,diff,diffSquare,gainLossCost}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::costMatrixType gainLossOptions::getCostMatrixTypeFromStr(const string& str)
{
	costMatrixType returnType;
	if (str == "file")
		returnType = file;
	else if (str=="fitch")
		returnType = fitch;
	else if (str=="diff")
		returnType = diff;
	else if (str=="diffSquare")
		returnType = diffSquare;
	else if (str=="gainLossCost")
		returnType = gainLossCost;
	else
		errorMsg::reportError("unknown type in MPoptions::getCostMatrixTypeFromStr- {file,fitch,diff,diffSquare,gainLossCost}");
	return returnType;
}


/********************************************************************************************
enum distributionType {GAMMA, GENERAL_GAMMA, UNIFORM,GAMMA_PLUS_INV, GENERAL_GAMMA_PLUS_INV, GAMMA_FIXED_CATEGORIES,GENERAL_GAMMA_FIXED_CATEGORIES, GAMMA_MIXTURE};
*********************************************************************************************/
string gainLossOptions::getDistributionType(distributionType type) 
{
	string res = "";
	switch (type)
	{
	case GAMMA_MIXTURE:
		res = "GAMMA_MIXTURE";
		break;
	case GAMMA_PLUS_INV:
		res = "GAMMA_PLUS_INV";
		break;
	case GENERAL_GAMMA_PLUS_INV:
		res = "GENERAL_GAMMA_PLUS_INV";
		break;
	case GAMMA_FIXED_CATEGORIES:
		res = "GAMMA_FIXED_CATEGORIES";
		break;
	case GENERAL_GAMMA_FIXED_CATEGORIES:
		res = "GENERAL_GAMMA_FIXED_CATEGORIES";
		break;
	case GENERAL_GAMMA:
		res = "GENERAL_GAMMA";
		break;
	case GAMMA:
		res = "GAMMA";
		break;
	case UNIFORM:
		res = "UNIFORM";
		break;

	default:
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType - {GAMMA, GENERAL_GAMMA, UNIFORM,GAMMA_PLUS_INV, GENERAL_GAMMA_PLUS_INV, GAMMA_FIXED_CATEGORIES,GENERAL_GAMMA_FIXED_CATEGORIES, GAMMA_MIXTURE}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::distributionType gainLossOptions::getDistributionType(const string& str) 
{
	if (str == "GAMMA_MIXTURE")
		return GAMMA_MIXTURE;
	if (str == "GAMMA_FIXED_CATEGORIES")
		return GAMMA_FIXED_CATEGORIES;
	if (str == "GENERAL_GAMMA_FIXED_CATEGORIES")
		return GENERAL_GAMMA_FIXED_CATEGORIES;
	if (str == "GENERAL_GAMMA_PLUS_INV")
		return GENERAL_GAMMA_PLUS_INV;
	if (str == "GAMMA_PLUS_INV")
		return GAMMA_PLUS_INV;
	if (str == "GENERAL_GAMMA")
		return GENERAL_GAMMA;
	else if (str == "GAMMA")
		return GAMMA;
	else if (str == "UNIFORM")
		return UNIFORM;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType - {GAMMA, GENERAL_GAMMA, UNIFORM,GAMMA_PLUS_INV, GENERAL_GAMMA_PLUS_INV, GAMMA_FIXED_CATEGORIES,GENERAL_GAMMA_FIXED_CATEGORIES, GAMMA_MIXTURE}");
	return GENERAL_GAMMA;
}
/********************************************************************************************
enum discretizationType {FIXED, QUANTILE, LAGUERRE};
*********************************************************************************************/
string gainLossOptions::getDiscretizationType(discretizationType type) 
{
	string res = "";
	switch (type)
	{
	case FIXED:
		res = "FIXED";
		break;
	case QUANTILE:
		res = "QUANTILE";
		break;
	case LAGUERRE:
		res = "LAGUERRE";
		break;
	default:
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType - {FIXED, QUANTILE, LAGUERRE}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::discretizationType gainLossOptions::getDiscretizationType(const string& str) 
{
	if (str == "FIXED")
		return FIXED;
	else if (str == "QUANTILE")
		return QUANTILE;
	else if (str == "LAGUERRE")
		return LAGUERRE;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType - {FIXED, QUANTILE, LAGUERRE}");
	return QUANTILE;
}
/********************************************************************************************
enum gammmaMixtureOptimizerAlgType {EM, ONE_DIM};
*********************************************************************************************/
string gainLossOptions::getGammmaMixtureOptimizerAlgType(gammmaMixtureOptimizerAlgType type) 
{
	string res = "";
	switch (type)
	{
	case ONE_DIM:
		res = "ONE_DIM";
		break;
	case EM:
		res = "EM";
		break;

	default:
		errorMsg::reportError("unknown type in gainLossOptions::getGammmaMixtureOptimizerAlgType - {EM, ONE_DIM}");
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
gainLossOptions::gammmaMixtureOptimizerAlgType gainLossOptions::getGammmaMixtureOptimizerAlgType(const string& str) 
{
	if (str == "ONE_DIM")
		return ONE_DIM;
	else if (str == "EM")
		return EM;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getGammmaMixtureOptimizerAlgType - {EM, ONE_DIM}");
	return EM;
}
/********************************************************************************************
enum treeSearchAlgType {njJC,njML,njJCOLD};
*********************************************************************************************/
string gainLossOptions::getTreeSearchAlgType(treeSearchAlgType type) 
{
	string res = "";
	switch (type)
	{
	case njJC:
		res = "njJC";
		break;
	case njML:
		res = "njML";
		break;
	case njJCOLD:
		res = "njJCOLD";
		break;

	default:
		errorMsg::reportError("unknown type in gainLossOptions::getTreeSearchAlgType - {njJC,njML,njJCOLD}");
	}
	return res;
}

//////////////////////////////////////////////////////////////////////////
gainLossOptions::treeSearchAlgType gainLossOptions::getTreeSearchAlgType(const string& str) 
{
	if (str == "njJC")
		return njJC;
	else if (str == "njML")
		return njML;
	else if (str == "njJCOLD")
		return njJCOLD;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getTreeSearchAlgAlgType - {njJC,njML,njJCOLD}");
	return njML;
}
/********************************************************************************************
enum characterFreqEvalType {FiftyFifty, LeavesAve, optimizeOverTree};
*********************************************************************************************/
string gainLossOptions::getCharacterFreqEvalType(characterFreqEvalType type) 
{
	string res = "";
	switch (type)
	{
	case optimizeOverTree:
		res = "optimizeOverTree";
		break;
	case LeavesAve:
		res = "LeavesAve";
		break;
	case FiftyFifty:
		res = "FiftyFifty";
		break;

	default:
		errorMsg::reportError("unknown type in gainLossOptions::getCharacterFreqEvalType - {FiftyFifty, LeavesAve, optimizeOverTree}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::characterFreqEvalType gainLossOptions::getCharacterFreqEvalType(const string& str) 
{
	if (str == "optimizeOverTree")
		return optimizeOverTree;
	else if (str == "LeavesAve")
		return LeavesAve;
	else if (str == "FiftyFifty")
		return FiftyFifty;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionTypeStr - {FiftyFifty, LeavesAve, optimizeOverTree}");
	return optimizeOverTree;
}

/********************************************************************************************
enum rateEstimationMethodType {ebExp, mlRate};
*********************************************************************************************/
string gainLossOptions::getRateEstimationMethodType(rateEstimationMethodType type) 
{
	string res = "";
	switch (type)
	{
	case mlRate:
		res = "mlRate";
		break;
	case ebExp:
		res = "ebExp";
		break;

	default:
		errorMsg::reportError("unknown type in gainLossOptions::getRateEstimationMethodType - {ebExp, mlRate}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::rateEstimationMethodType gainLossOptions::getRateEstimationMethodType(const string& str) 
{
	if (str == "ebExp")
		return ebExp;
	else if (str == "mlRate")
		return mlRate;
	else
		errorMsg::reportError("unknown type in gainLossOptions::getRateEstimationMethodType - {ebExp, mlRate}");
	return ebExp;
}

/********************************************************************************************
Types
enum simulationType {Uniform, Normal, Gamma, MPestEmp GammaNoise, MPratio}
*********************************************************************************************/
string gainLossOptions::getSimulationType(simulationType type)
{
	string res = "";
	switch (type) 
	{
	case Uniform:
		res = "Uniform";
		break;
	case Normal:
		res = "Normal";
		break;
	case Gamma:
		res = "Gamma";
		break;
	case MPestEmp:
		res = "MPestEmp";
		break;
	case SMestEmp:
		res = "SMestEmp";
		break;
	case GammaNoise:
		res = "GammaNoise";
		break;
	case EQ_gEql:
		res = "EQ_gEql";
		break;
	case EQ_gVrl:
		res = "EQ_gVrl";
		break;
	case Gam_gEql:
		res = "Gam_gEql";
		break;
	case Gam_gVrl:
		res = "Gam_gVrl";
		break;
	default:
		errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
gainLossOptions::simulationType gainLossOptions::getSimulationTypeFromStr(const string& str)
{
	simulationType returnType;
	if (str == "Uniform")
		returnType = Uniform;
	else if (str=="Normal")
		returnType = Normal;
	else if (str=="Gamma")
		returnType = Gamma;
	else if (str=="MPestEmp")
		returnType = MPestEmp;
	else if (str=="SMestEmp")
		returnType = SMestEmp;
	else if (str=="GammaNoise")
		returnType = GammaNoise;
	else if (str=="EQ_gEql")
		returnType = EQ_gEql;
	else if (str=="EQ_gVrl")
		returnType = EQ_gVrl;
	else if (str=="Gam_gEql")
		returnType = Gam_gEql;
	else if (str=="Gam_gVrl")
		returnType = Gam_gVrl;
	else
		errorMsg::reportError("unknown type in gainLossOptions::optimizationLevel- {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise}");
	return returnType;
}



