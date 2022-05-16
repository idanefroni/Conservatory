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


#ifndef __gainLossOptionsParams_OPTION
#define __gainLossOptionsParams_OPTION

#include "definitions.h"
#include <string>
#include <fstream>
using namespace std;

/*
--- utilize CLASS:Parameters ---
USAGE: SETTING DEFAULT PARAMETERS
Note that the type of the parameter is set according to the addParameter arguments. 
e.g., If a parameter is set using addParameter with an integer argument then subsequent updates (using updateParameter) 
to the same parameter will all be stored as integers. 
Therefore the following code should output a 0:
EXAMPLE
Parameters::addParameter("Dummy", 3);
Parameters::updateParameter("Dummy", "This should set it to zero");
cout << Parameters::getString("Dummy");
END
Note also that when setting default values of float parameters always use
a decimal point or else these parameters will be added as integers. 

USAGE: READING PARAMETERS FROM FILE
The readParameters method receives an input stream from which parameters are to be read. 
Files are structured so that each line specifies the value of a parameter. 
Each line gives the parameter name, a white space and then the parameter value. 
Lines whose first non white-space charachter is # are ignored. 
A basic schema for using the Parameters class is to set the default
values using addParameter calls and then calling readParameters to read in
parameters with other values or new parameters. 
EXAMPLE
Parameters::addParameter("CubeSize", 1.0);
Parameters::addParameter("MinVote", 8);
ifstream params("params");
Parameters::readParameters(params);
params.close();
Parameters::dump(cout);
END
With the following parameters file:
EXAMPLE
CubeSize 0.5
File  pdb4hhb.ent
END
The following output should result:
EXAMPLE
CubeSize (Float) 0.5
File     (Str)   pdb4hhb.ent
MinVote  (Int)   8
END

USAGE: SUBCLASSING AND PERFORMANCE
The Parameters engine keeps the parameters in a sorted list.
The correct usage would have been to inherit: e.g., class ProgParams : protected Parameters

*/


class gainLossOptions{
public:
	enum discretizationType {FIXED, QUANTILE, LAGUERRE};
	enum distributionType {GAMMA, GENERAL_GAMMA, UNIFORM,GAMMA_PLUS_INV, GENERAL_GAMMA_PLUS_INV, GAMMA_FIXED_CATEGORIES,GENERAL_GAMMA_FIXED_CATEGORIES, GAMMA_MIXTURE};
	enum treeSearchAlgType {njJC,njML,njJCOLD};
	enum rateEstimationMethodType {ebExp, mlRate};
	enum characterFreqEvalType {FiftyFifty, LeavesAve, optimizeOverTree};
	enum gammmaMixtureOptimizerAlgType {EM, ONE_DIM};
	enum costMatrixType {file,fitch,diff,diffSquare,gainLossCost};
	enum optimizationLevel {VVVlow,VVlow, Vlow, low, mid, high, Vhigh};
	enum simulationType {Uniform, Normal, Gamma, MPestEmp,SMestEmp,GammaNoise  ,EQ_gEql,EQ_gVrl,Gam_gEql,Gam_gVrl};

	//enum optimizeBranchLengthsType {noBBL, mlBBLUniform, mlAndAlphaBBL};
	
public:
	virtual ~gainLossOptions();
	
	static void initOptions(const string& paramFileName);
	static void initDefault();
	static void readParameters(const string& paramFileName);
	static void getParamsFromFile(const string& paramFileName);
	static void getOutDirFromFile(const string& paramFileName);
	static void verifyConsistParams();
	ostream& out() const {return *_outPtr;};
	
// conversions from enum to (from) string 
	static string getDistributionType(distributionType type);
	static distributionType getDistributionType(const string& str);
	static characterFreqEvalType getCharacterFreqEvalType(const string& str); 
	static string getCharacterFreqEvalType(characterFreqEvalType type);
	static string getRateEstimationMethodType(rateEstimationMethodType type);
	static rateEstimationMethodType getRateEstimationMethodType(const string& str);
	static string getGammmaMixtureOptimizerAlgType(gammmaMixtureOptimizerAlgType type);
	static gammmaMixtureOptimizerAlgType getGammmaMixtureOptimizerAlgType(const string& str);
	static string getTreeSearchAlgType(treeSearchAlgType type);
	static treeSearchAlgType getTreeSearchAlgType(const string& str);
	static string getDiscretizationType(discretizationType type); 
	static discretizationType getDiscretizationType(const string& str);
	static string getCostMatrixType(costMatrixType type);
	static costMatrixType getCostMatrixTypeFromStr(const string& str);	
	static string getOptimizationLevelType(optimizationLevel type);
	static optimizationLevel getOptimizationLevelTypeFromStr(const string& str);

	static string getSimulationType(simulationType type);
	static simulationType getSimulationTypeFromStr(const string& str);

	static void readFromParameters2gainLossOptions();
// update parameters by dependencies
	static void updateDependencies();
	static void updateOptimizationLevel(optimizationLevel level);
	static void updateUserGainLossRatio(MDOUBLE gainLossRatio);
	static void updateGainLossAsFreq();
	static void updateGainEQloss();
	static void updateKeepUserGainLossRatio();
	static void updateRemoveComputationNotSuiteForModels();
	static void updateGainLossDist();
	static void updateAccountForMissingData();
	static void updateInitParamsAtRandPointsInSimPostExp();
	static void updateSimulatePosteriorExpectationOfChange();
	static void updateOnlyComputeLikelihood();
	static void updateFlatUserParameters();
	static void updateNoBBL();
	static void updateNoBranchLengthDiffComputation();
	static void updateNoOptimization();
	static void updatNoSeq();
	static void updateParamsInRangeOverrideParamFile();
	static void updatParametericBootstrapComputationOfCorrelation();
	static void updateParsimonyRun();



public:
//################### Basic parameters:
// input (general)
	static string _seqFile;		// essential - fasta file with presence(1)/absence(0) for each species over all gene families (positions)
	static string _treeFile;	// basic	 - if not given - calculated based on distanceTable
	static string _treeFileOrig;	// // used for branchDiff calc. functionality

	static string _rootAt;		// name of node to be root (the tree must contain names of internal nodes)
	static string _referenceSeq; // the results are printed with this seq in each positions. (default - first)
	//static string _mainType;

// output
	static string _outDir;		// _outDir = "RESULTS", concatenated after current dir location 'pwd'
	static string _logFile;		// print-outs of the running progress including the estimated parameters optimization
	static int _logValue;		// verbosity level - ~4 - normal, >7 - load of info
	static string _treeOutFile;				// "TheTree.ph" - tree after BBL and other changes - 
// all of these files are still part of the output, but names are fixed
	//static string _outFile;		// Rate4Site results (normalized - Ave=0, Sd=1)
	//static string _outFileNotNormalize;		// Rate4Site results (original)
	//static string _outFileGain4Site;		// gain4Site results
	//static string _outFileLoss4Site;		// loss4Site results
	//static string _outFileLikeofPos;		// compare to model with gainRate=0
	//static string _outFilePosteriorExpectationOfChange;		// exp01, exp10 per gene



//################################################## Model params
	static int _alphabet_size;	// 2 - presence(1)/absence(0)
	static bool _gainLossDist;				// GLM (mixture)
	static bool _accountForMissingData;		// for phyletic patterns - must be true
	static int _minNumOfOnes;				// for COG and EggNOG only patterns with 3 or more are observable
	static int _minNumOfZeros;				// for indels, there is no position with only 1s => minNumOfZeros=1

	static bool _gainEQloss;				// M1 (the basic model)
	static bool _isReversible;				// if _isReversible = False -> the root is fixed
	static bool _isRootFreqEQstationary;	// same "-"
	static bool _gainLossDistPlusInvariant;	// Automatically True if GENERAL_GAMMA_PLUS_INV or GAMMA_PLUS_INV
	static bool _gainLossRateAreFreq;		// test parameter where gain+loss = 1, and the "r_Q" is external

//Each of the rates governing the stochastic process are assumed to be sampled from a prior distribution.
	static distributionType _rateDistributionType;
	static distributionType _gainDistributionType; //(only for the mixture models - _gainLossDist 1) 
	static distributionType _lossDistributionType; //(only for the mixture models - _gainLossDist 1)
	static int _numberOfGainCategories;		//	gain 3-5 - the overall number of stochasticProcess 9-25
	static int _numberOfLossCategories;		//	loss 3-5
	static int _numberOfRateCategories;		//	discretization usually 4-16
	static int _numberOfRateComponents;		//	gammaMix
	static discretizationType _rateDiscretizationType;		// QUANTILE, LAGUERRE - only in use for gammaMix


//################################################## computations
	static bool _calculateRate4site;
	static rateEstimationMethodType _rateEstimationMethod;	// mlRate (only option for UNIFORM) or posteriorBayesianExpectation
	static bool _calculeGainLoss4site;
	static bool _calculePosteriorExpectationOfChange;
	static bool _calculateAncestralReconstruct;
	static bool _simulatePosteriorExpectationOfChange;		// simulate PostExp (To test to accuracy of the stochastic mapping)
	static bool _isOnlySimulateSeq;							// no mapping or parsimony is done

	static bool _simulateSequences;							// Test the rate4site computation
	static bool _calculateRate4siteSim;						// Test the rate4site computation
	static bool _calculeBranchLegthDiffFactor;				// if BBL is used for each branch - compare length before/after
	static bool _findCoEvolvingSitesOldNotWorking;						// for the co evolving project
	static bool _printAncestralReconstructPosterior;
	static bool _saveProbChanges_PosNodeXY;					// used for AnsetralReconstruc - posterior
	static bool _isComputeDistanceFromRootForRecent;		// used to classify branches

//################################################## Prints
	static bool _printLikelihoodLandscapeAlphaRate;
	static bool _printLikelihoodLandscapeGainLoss;
	static bool _printLikelihoodLandscapeTheta;
	static bool _optAlphaInIteration;
	static bool _optBBL_LS_InIteration;
	static bool _optBBL_EM_InIteration;
	static bool _printTree;
	static bool _printSeq;
	static bool _printPij_t;
	static bool _printLofPos;
	static bool _printLofPosBothModels;
	static bool _printTreesWithProbabilityValuesAsBP;	// tree for each position
	static bool _printTreesWithExpectationValuesAsBP;	// tree for each position
	static bool _printTreesWithAncestralReconstructAsBP;// tree for each position
	static bool _printPropExpOfChangeFullData;			// huge file...
	static bool _printExpPerPosPerBranchMatrix;			// Used as input for COMAP
	static bool _printComputedCorrelations;				// Correlation
	static bool _performParametricBootstapCorrelation;	// Correlation with simulation as correction
	static bool _usePosSpecificSimulations;				// pos-specific simulation using startSimultePosteriorExpectationOfChange
	static bool _isAddSimulationsWithLowRate;			// Correlation with simulation as correction
	static bool _isFDRcorrectionForPValInCorrelation;	//
	static bool _isComputeQVals;						// qVals are printed
	static MDOUBLE _pValueCutOffForBootStrap;	//0.05, 0.01
	static bool _isConsiderNegativeCorrelations;
	static int _numOfBinsInParametricBootstrapSimulations;
	static bool _isDivideBinsByRange;					// if true, each bin will get different number of samples, but the rate(Nmin) is eq-partitioned
	static bool _isSortVectorOfCorrelationsBinsByLowerRateBound;	// it true, each pair pVal is computed according to all simulation with Nmin >= that of pair ()
	static bool _isSortVectorOfCorrelationsBinsByMidRateBound;		// if ture, the bins are overlapping
	static MDOUBLE _relativeSizeOfOverLappedBins; // if 0.5, 50% of samples per bin

	static bool _isPrintpairWiseCorrelationsAndNmin;	// util, for statistics
	static bool _isPrintCorrelationsOfAllPairs_Corr;	// Huge file
	static bool _isPrintCorrelationsOfAllPairs_pVal;	// Huge file

	static bool _isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH; // only pairs with PVal significant after BH will be printed
	static bool _isAllCorrTypeReqruiedToBeSignificant;		// if true, sufficiet that one corType results with pVal>BH[corType] not to print
	static bool _isNminBasedOnCountBranchesOverCutOff;		// if true, Nmin is based on numOfEvent>cutoff, not total expectation
	static MDOUBLE _minExpThresholdForPValComputationForCorrelatingPair;	// 0, 2,3,..
	static bool	 _isUpdateMinExpThresholdGivenSimulaitonsQuantile;		// After simulation, minR is defined by 0.25 quantile in simulation (updated only if higher)
	static bool	 _isUpdateMinExpThresholdGivenRealDataQuantile;			// Given real data, minR is defined by the 0.1 percentile (updated only is higher)
	static MDOUBLE _updateMinExpThresholdGivenRealDataQuantileVal;			// if 0.2, Nmin is for sites above the 0.2 percentile rate

	static bool	_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel;	// After correlation of simulated data is computed minR is elevated to P(corr=1)<
	static bool _isCompExtremeValDistribution;							// pValue is also estimated assuming EVD distribution

	static MDOUBLE _minExpThresholdAsPercentFromNumOfSpeciesForPValComputationForCorrelatingPair;	// e.g., if =2, with 500 species, minT = 10	

	static bool _isCorrelateWithPearson;			// Pearson or Spearman's correlation computed for CoEvolution
	static bool _isCorrelateWithSpearman;			
	static bool _isCorrelationsBasedOnMaxParsimonyMapping;				

	static bool _isAlsoCorrelateWithLoss;				// additionally to gain, compute with loss vectors
	static bool _isAlsoCorrelateWithBoth;				// additionally to gain and loss, compute with a gain . loss concatenated vectors
	static bool _isOnlyCorrelateWithBoth;				// compute with a gain . loss concatenated vectors, only
	static bool _isUseRateForSiteAsNminForCorrelations;
	static bool _isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur; // Remove simulated position with too low/high occur to save later computation time (quick and (very)dirty)
	static bool _isRemoveSimulatedPositionsBasedOnMP;						// Remove simulated positions with less than 2 events based on max parsimony (quick and dirty)
	static MDOUBLE _minNumOfMPEvent2RemoveSimulatedPositions;					// If 1 then gain+loss events must be >=1
	static bool _isUpdateminNumOfMPEvent2RemoveSimulatedPositions;				// If true, add 0.2 events for every sqrt(num Of species)


	static bool _printComputedCorrelationsAllSites;		// all-against-all, in STRING format
	static string _selectedSitesForCorrelation;		// in this file, for each position, the correlation with all other positions if computed.
	static bool _isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation; // the last is a trait (with possible unknown)
	static int  _checkCoEvolWithUnionPAP_against_pos;			// PAP will be modified to union (1 in either) with selected position

	static bool _isIgnoreCorrelationAmongSelectedSites;
	static bool _isNormalizeForBranchExpInCorrCompute;
	static bool _isNormalizeByExpectationPerBranch;	// else, by branch length


	static bool _printAncestralReconstructFullData;		// huge file...
	static bool _printDEBUGinfo;						// huge file...
	static bool _printLikelihoodLandscape;				// test purpose (Ad-hoc)
	static MDOUBLE _likelihoodLandscapeIncrement;
	static bool _printP11forgain;						// test purpose (Ad-hoc)
	
//################################################## optimizations	
	static bool _isInitGainLossByEmpiricalFreq;				// the sp is initialized with the empirical 0 and 1 freq
	static bool _isBBLEMwithSimpleSpBeforeFullOptimization;	// before optimization - BBL-EM is performed with simplified sp
	static bool _isSkipFirstParamsOptimization;
	static bool _isOptimizeParamsWithLogMinMax;		// when the parameter is a positive and values are e.g., [0.01,100] brent works better for [-2,2]


	static bool _performOptimizations;
	static bool _performOptimizationsBBL;
	static bool _performOptimizationsBBLOnlyOnce;
	static bool _isLongAndAccurateOptimization;

	static bool _isBblLS;
	static bool _isbblLSWhenbblEMdontImprove;
	static bool _isSkipBblEMWhenbblEMdontImprove;

	static bool _isBblEMbeforeLSWithMissSpecifiedModel;

	static bool _isBblForceFactorCorrection;
	static MDOUBLE _BblFactorCorrection;

	static bool _isOptimizeGainLossRatioInsteadOfGainAndLossSeperately;
	static bool _isOptimizeInvariantCategoryProb;	
	static bool _isUpdateOnlyGainBetaForRatio;		// currently, not in use
	static bool _isComputeLikelihoodDuringInit;		// true, unless fast/parsimony run is performed	

	static bool _isMultipleAllBranchesByFactorAtStart;
	static bool _isNormalizeAtStart;

	static bool _performOptimizationsROOT;
	static bool _performOptimizationsManyStarts;
	static bool _performOptimizationsBBLManyStarts;
	static bool _correctOptimizationEpsilon;		// according to dataset size (initial likelihood)
	static bool _simulatedAnnealing;				// epsilon is lowered with iterations
	static MDOUBLE _simulatedAnnealingMinEpsilonFactor;		// to lower to normal epsilons (Model, BBL, Both)
	static MDOUBLE _simulatedAnnealingCoolingFactor;		// to lower epsilons each iteration

	static gammmaMixtureOptimizerAlgType _gammmaMixtureOptimizerAlg; // ONE_DIM or EM (not fully functional)
	static characterFreqEvalType _characterFreqEval;		// "-F option" the estimation of freq at root: FiftyFifty, LeavesAve, optimizeOverTree

	static bool _isStartWithTheta;			// the optimization loop of the parameter will start with Theta
	static bool _isSkipGainOptimization;	// 
	
	static MDOUBLE _epsilonOptimizationThetaFactor;			// the optimization loop of the parameter will start with Theta
	static bool _isAlphaLimit;				// 0.3 - for Alpha <<0.3, the following computations are erroneous [BUG?]
	static bool _isGainLimit;				// 0.1 - for Gain <<0.1, the following computations are erroneous [BUG?]
	static bool _isHGT_normal_Pij;			// test parameter - 
	static bool _isHGT_with_Q;				// test parameter - 
	static bool _incrementFactorForGain;	// test parameter - 
	static bool _lossBiggerGainLimit;		// test parameter - 
	static MDOUBLE _slopeFactorForGain;		// test parameter - limit growth in gain estimation

	static optimizationLevel _optimizationLevel;	// change all epsilons and related parameters
	static MDOUBLE _epsilonOptimizationIterationCycle;		//if the log-likelihood after optimization is lower than this threshold - then optimize again.
	static MDOUBLE _epsilonOptimizationModel; 
	static MDOUBLE _epsilonOptimizationBBL;
	static MDOUBLE _epsilonOptimizationIterationCycleManyStarts;

	static MDOUBLE _epsilonFactor_Model;
	static MDOUBLE _epsilonFactor_BBL;
	static MDOUBLE _numIterationsFactor_Model;
	static MDOUBLE _numIterationsFactor_BBL;

	static int _maxNumOfIterations;			//	over Model,Root, and BBL
	static int _maxNumOfIterationsModel;	 
	static int _maxNumOfIterationsBBL;
	static int _maxNumOfIterationsManyStarts;	// the basic number of manyStarts option (Model and BBL factors are used)

	static MDOUBLE _epsilonForReRootFactor;			// only for substantial improvement the tree will be re-rooted
	static MDOUBLE _percentOfImprovManySarts;		// epsilonOptimization = abs(logL)*_percentOfImprovManySarts
	static MDOUBLE _percentOfImprov;				// epsilonOptimization = abs(logL)*_percentOfImprov

	static bool _initParamsAtRandPoints;
	static bool _initParamsAtRandPointsInOptimization;
	static bool _initRandomGammaMixuteParam;
	static int	_numberOfRandPointsInOptimization;
	static int	_numberOfRandStartPoints;

// all the model parameters can be given by the user
	static MDOUBLE _userGainLossRatio;
	static bool _keepUserGainLossRatio;
	static MDOUBLE _userGain;
	static MDOUBLE _userLoss;
	static MDOUBLE _userTheta;						// default 0.5 - otherwise, counting is done prior to optimization
	static MDOUBLE _userAlphaGain;
	static MDOUBLE _userBetaGain;
	static MDOUBLE _userProbInvariantGain;
	static MDOUBLE _userAlphaLoss;
	static MDOUBLE _userBetaLoss;
	static MDOUBLE _userProbInvariantLoss;
	static MDOUBLE _userAlphaRate;
	static MDOUBLE _userBetaRate;
	static MDOUBLE _userProbInvariantRate;
	static MDOUBLE _userRateInvariantVal;			// The low (~10-8) value that corresponds to rate=0

// for initRand - Rand(x){min<x<max}
	static MDOUBLE _userGainMax;
	static MDOUBLE _userLossMax;	
	static MDOUBLE _userThetaMax;
	static MDOUBLE _userAlphaGainMax;	
	static MDOUBLE _userBetaGainMax;
	static MDOUBLE _userProbInvariantGainMax;
	static MDOUBLE _userAlphaLossMax;
	static MDOUBLE _userBetaLossMax;
	static MDOUBLE _userProbInvariantLossMax;
	static MDOUBLE _userAlphaRateMax;
	static MDOUBLE _userBetaRateMax;
	static MDOUBLE _userProbInvariantRateMax;
	
	static MDOUBLE _userGainMin;
	static MDOUBLE _userLossMin;
	static MDOUBLE _userThetaMin;
	static MDOUBLE _userAlphaGainMin;
	static MDOUBLE _userBetaGainMin;
	static MDOUBLE _userProbInvariantGainMin;
	static MDOUBLE _userAlphaLossMin;
	static MDOUBLE _userBetaLossMin;
	static MDOUBLE _userProbInvariantLossMin;
	static MDOUBLE _userAlphaRateMin;
	static MDOUBLE _userBetaRateMin;
	static MDOUBLE _userProbInvariantRateMin;

//################################################## PostExp (Stochastic mapping based Counting)
	static int _numOfSimulationsForPotExp;	// the counting (expectation) is based on simulations - val: >1000 - accurate enough
	//static MDOUBLE _probCutOffSum;				// the cutOff to sum count (0.5) "ProbabilityPerPos.txt", "ProbabilityPerPosPerBranch.txt"
	static bool _isFewCutOffCounts;				// Few Cut offs, not just one
	static MDOUBLE _probCutOffCounts;			// the cutOff to estimate HGT count (0.6) "gainLossProbExpCountPerPos.txt"
	static MDOUBLE _probCutOffPrintEvent;		// the cutOff for perPosperBranch (so that file is not too big) (0.05)


//################################################## simulate PostExp (To test to accuracy of the stochastic mapping)
	static simulationType _simulationType;	// {Uniform, Normal, Gamma, MPestEmp, SMestEmp}
	static bool _isMPratio;
	static int _numberOfPositions2simulate;
	static int _numberOfIterations2simulate;
	static int _numberOfIterationsForPrintResults;	// if =3, each 3 simulation iterations, results are updated (thus, temp results are available)
	static MDOUBLE _percentileOfNminWithCorr1RequiredForLastIteration;

	static bool _modelOptimizationSimPostExp;
	static bool _BBLOptimizationSimPostExp;
	static MDOUBLE  _epsilonOptForPostExpSimFactor;			// reduce optimization run-time in simulations
	static MDOUBLE  _numOfIterationsOptForPostExpSimFactor;	// reduce optimization run-time in simulations
	static MDOUBLE  _loss2gainRatioToSim;
	static bool _isInitGainLossByEmpiricalFreqSimulatePostExp;		// the sp is initialized with the empirical 0 and 1 freq
	static bool _is3states;
	static MDOUBLE _3statesGain;
	static MDOUBLE _3statesMore;
	static MDOUBLE _3statesLess;
	static MDOUBLE _3statesLoss;
	static MDOUBLE _3states0;
	static MDOUBLE _3states1;

	//Used as.... enum simulationType {GAMMA, UNI, MP}
	static bool _isFlatTreeBeforOpt;						// Flat the tree before model-based estimation
	static bool _isbBLEMwithSimpleSpSimulatePostExp;
	static MDOUBLE _noiseLevelInGammaSimulation;
	static bool _initParamsAtRandPointsInSimPostExp;		// gain, loss rates are sampled uniform distribution
	static bool _isMatrixGainLossFromRatioInSimulations;		// 

	static bool _initRootFreqAtRandPointsInSimPostExpEachPos;	// not required
	static bool _isTheataFromObservedFreq;					// The theta is taken from observed freq +random perturbation
	static bool _isRootFreqEQstationaryInSimulations;					
	static bool _isFlatSpBeforeOpt;							// need to change to T when performing initParamsFromTrueEstimation

//################################################## CoEvolvingSites
	static int _numberOfSequences2simulate;
	static int _numberOfSequences2simulateForCoEvol; // number of simulations used in the co-evoving computations val: >1000 - accurate enough
	static bool _useTheSameSpForSim;
	static bool _isReversibleSim;
	static distributionType _rateDistributionTypeSim;
	static bool _gainEQlossSim;
	static bool _writeSeqSim;

//################################################## Misc.
	static MDOUBLE _maxRateForML;
	static MDOUBLE _minBranchLength;
	static MDOUBLE _maxBranchLength;

	static treeSearchAlgType _treeSearchAlg;	// To construct tree from distanceTable (JC or others)
	static Vdouble* _weights;	// positions are weighted (not in use)
	static bool _isSequenceUniqPattern;
	static bool _isRemovePositionsWithHighPercentOfMissingData;
	static MDOUBLE  _fractionOfMissingDataToRemove;

	static bool _isOnlyComputeLikelihood;
	static bool _isAnaliticComputeJumps;
	static bool _isNormalizeQ;
	static bool _isNormalizeQinSpVVec;
	static bool _isNormalizeQandTreeafterOpt;
	static bool _isFlatUserParameters;
	static bool _isAlphaEqBetaManipulation;		// Turn GeneralGamma into Gamma -> Alpha=Beta
	static bool _calculeBranchLegthDiffFactorFromInputTrees;	// input 2 trees - compute logL diff per branch length
	static bool _intersectTreeAndSeq;	// input tree and seq (not the same taxa) - intersect, write seq and tree

	static bool _isOnlyParsimony;
	static bool _calculeMaxParsimonyChange;
	static bool _calculeMaxParsimonyChangeSeveralGainLossRatios;
	static string _costMatrixfile;
	static costMatrixType _costMatrixType;
	static MDOUBLE _costMatrixGainLossRatio;

private:
	static ostream* _outPtr;
	//static ofstream _out_f;

};
#endif
