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


#ifndef ___GAIN_LOSS_
#define ___GAIN_LOSS_

#include "aaJC.h"
#include "bblEM.h"
#include "bestAlpha.h"
#include "chebyshevAccelerator.h"
#include "checkcovFanctors.h"
#include "checkcovFanctorsWithFactors.h"
#include "definitions.h"
#include "distanceTable.h"
#include "distributionPlusInvariant.h"
#include "errorMsg.h"
#include "evaluateCharacterFreq.h"
#include "fastStartTree.h"
#include "gainLossAlphabet.h"
#include "gainLossModel.h"
#include "gainLossUtils.h"
#include "ancestralReconstructStates.h"
#include "gammaDistribution.h"
#include "generalGammaDistribution.h"
#include "generalGammaDistributionPlusInvariant.h"
#include "jcDistance.h"
#include "likeDist.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include "logFile.h"
#include "matrixUtils.h"
#include "nj.h"
#include "nucJC.h"
#include "numRec.h"
#include "optimizeGainLossModel.h"
#include "optimizeGainLossModelVV.h"
#include "readDatMatrix.h"
#include "recognizeFormat.h"
#include "seqContainerTreeMap.h"
#include "sequence.h"
#include "sequenceContainer.h"
#include "siteSpecificRate.h"
#include "someUtil.h"
#include "stochasticProcess.h"
#include "tree.h"
#include "treeIt.h"
#include "trivialAccelerator.h"
#include "uniDistribution.h"
#include "unObservableData.h"

#include <cassert>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <time.h>
#include <algorithm>

#ifdef WIN32	
#include <process.h>
#else
#include <unistd.h>
#endif

class gainLoss {

public:
	explicit gainLoss();
	virtual ~gainLoss();
	void run();

	
private:
	void initialize(bool isComputeLikelihood=true);
	void initializeBranchLengthDiff();
	void initializeUnObservableData();

	void fillOptionsParameters(int argc, char* argv[]);
	void printOptionParameters(ostream& out= cout);
	
	void startSequenceContainer();
	void checkMinNumOfOnesOrZeros(sequenceContainer&  sc, int minNumOfOnes, int  minNumOfZeros, bool isRemovePosNotWithinMinMax=false, bool isReportRemovedPos=false);
	void produceUnionPAP_against_pos(sequenceContainer&  sc, int pos_for_union, bool is_ignore_last_pos=true);

	void startSequenceContainerUniqPatterns();
	void countOccurPerPos();
	void removePositionsWithHighPercentOfMissingData(MDOUBLE PercentOfMissingDataToRemove);

	void startStochasticProcess(bool gainLossDist);
	void setRootFreq();
	void startStochasticProcess();
	stochasticProcess*  startStochasticProcessGeneric(gainLossOptions::distributionType rateDistributionType, const bool isReversible);
	void startStochasticProcessVec();

	void startEvolTreeTopology(ostream& out=cout);
	void startOptimizations();
	void startRate4Site(sequenceContainer& sc, tree& tr, stochasticProcess* sp,  string& outDir, unObservableData* unObservableData_p);
	void startGainLoss4Site(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> > spVVec,distribution* gainDist,distribution* lossDist,
		string& outDir, unObservableData* unObservableData_p);

	void computePosteriorExpectationOfChangeRunOnly();
	void startComputePosteriorExpectationOfChange();	
	void startComputePosteriorExpectationOfChange(sequenceContainer& sc, tree& tr, stochasticProcess* sp, VVdouble LpostPerCat, unObservableData* unObservableData_p, string& outDir,MDOUBLE distanceFromNearestOTUForRecent,bool isUpdateMPPerPos=true);
	void startComputePosteriorExpectationOfChange(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution* gainDist, distribution* lossDist, VVVdouble& LpostPerSpPerCat,unObservableData* unObservableData_p, string& outDir,MDOUBLE distanceFromNearestOTUForRecent,bool isUpdateMPPerPos=true);

	void startComputeAmongSitesCorrelations();
	void computeCoEvolutionScoresBasedOnSimulatedData(sequenceContainer& scSimulated);
	void startParametricBootstapCorrelation();
	int computeCoEvolutionScoresBasedOnSimulatedDataCoMap(sequenceContainer& scSimulated,tree& trSampled ,MDOUBLE qNminOfRealData, bool& isLastIteration, int& numOfpairsWithRateAboveMinRequiredExp, MDOUBLE& T_BH_prev, ofstream* simCorrelStream);


	void startMaxParsimonyChange(bool isUpdateMPPerPos=true);
	void startMaxParsimonyChange(sequenceContainer& sc, tree& tr, string& outDir,MDOUBLE costMatrixGainLossRatio, MDOUBLE distanceFromRootForRecent,bool isUpdateMPPerPos=true);

	void startSimulateSequences(int numOfSequenceSets, int seqLengthInSet); // if default=0, take length for input sequence
	
	void startSimultePosteriorExpectationOfChange(int numOfSequenceSets=5, const int numOfRepeats=1);
	MDOUBLE ComputeEmpiricalExpectedQforStationaryProcess(VVVdouble& EmpPerPos, MDOUBLE minRate=0.01);	
	
	//void simultePhyleticData(const int numOfSequenceSets, string strSeqFirst,MDOUBLE loss2gainRatioToSim, gainLossOptions::simulationType simulationType
	//	, MDOUBLE AlphaGain, MDOUBLE BetaGain, MDOUBLE AlphaLoss, MDOUBLE BetaLoss, MDOUBLE AlphaRate);

	void FlatSpBeforeOpt(stochasticProcess& sp , unObservableData* unObservableData_p);
	void FlatSpBeforeOpt(vector<vector<stochasticProcess*> >& spVVec,distribution * gainDist, distribution * lossDist, unObservableData* unObservableData_p);

	void getStartingTreeFromTreeFile();
	void getStartingTreeNJ_fromDistances(const VVdouble& disTab,const vector<string>& vNames);

	void fillReferenceSequence();
	Vdouble computeFreq();
	
	void optimizationsManyStarts(const MDOUBLE epsilonOptimization, const int numIterations);
	void optimizationsManyStartsNoVec(const MDOUBLE epsilonOptimization, const int numIterations);

	void optimizationsVVManyStarts(const MDOUBLE epsilonOptimization, const int numIterations);

	void optimizations(ostream& out =cout);
	void printModellValuesOfParams();
	void printModellValuesOfParams(stochasticProcess* sp, tree& tr);
	void printModellValuesOfParams(tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution* gainDist, distribution* lossDist);


	void optimizationsSPvv(ostream& out =cout);
	MDOUBLE  optimizeParameters(ostream& out =cout);
	MDOUBLE  optimizeParametersSPvv(ostream& out =cout);
	MDOUBLE optimizeBranchLengths();
	void normalizeQandTree(bool isComputeLikelihood=true, bool isMultipleAllBranchesByNormFactor= true);	// normalizeQ or normalizeMatrices and the corresponding tree
	void convertGainLossRatesToFreq();
	void AlphaEqBetaManipulation();

	void printPij_t(MDOUBLE dist=0.1,ostream& out= cout);
	void printQ(ostream& out= cout);
	void printTreeLikelihoodAllPosAlphTheSame(bool isLOGnOUT = true,ostream& out =cout);
	void printLofPos();
	MDOUBLE printLofPos(ostream& out);
	void printLofPosBothModels();
	MDOUBLE printLofPosBothModels(ostream& out);

	void printLikelihoodLandscape(stochasticProcess* sp);
	void printLikelihoodLandscapeStatFreqRatioAndRootFreqRatio();

	void computeAveAndStd();
	void normalizeRates();
	void printRatesML(ostream& out, const Vdouble & rate2print);
	void printRatesBayes(ostream& out, const Vdouble & rate2print);
	void printAveAndStd(ostream& out= cout);
	
	Vdouble computeRate4site(); // needed also for computePosteriorExpectationOfChangePerSite (if not run befor)
	void printRates(ostream & out, const Vdouble & rate2print); // needed also for gammaMix


	void printGainLossBayes(ostream& out, const Vdouble& rate2printV, const Vdouble& lowerBoundV, const Vdouble& upperBoundV,const VVdouble& posteriorV, const distribution* dist);

	void initParamsAtRandPoints(int numOfRandPoints, stochasticProcess* sp, unObservableData* currUnObs, ostream& out=cout);
	void initParamsAtRandPointsSPvv(int numOfRandPoints, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist, unObservableData* currUnObs,ostream& out =cout);
	//void initParamsAtIntervalPoints(int pointIndex,int numOfRandPoints, stochasticProcess* sp, unObservableData* currUnObs, ostream& out);


	void computePosteriorExpectationOfChangePerSite(Vdouble& expV01, Vdouble& expV10);
	void initMixtureParams(Vdouble& initAlphaRates, Vdouble& initBetaRates, Vdouble& initCompProbRates, int numOfGammaComp,
		MDOUBLE initAlphaRate=1, MDOUBLE initBetaRate=1, MDOUBLE initCompProbRate=1);
	void printGainLossProbabilityPerPosPerBranch(int pos, MDOUBLE probCutOff, VVVdouble& probChanges, ostream& out=cout, ostream& outCount=cout);
	void printGainLossExpectationPerBranch(VVVdouble& probChanges, ostream& out=cout);
	void computeBranchLegthDiffFactor(ostream& out=cout);
	//void initMissingDataInfo();
	vector<sequenceContainer>  simulateSequences(int numOfSequenceSets, int seqLengthInSet, bool writeSeq,
		bool useTheSame, bool isReversible, bool isGeqL, gainLossOptions::distributionType rateDistributionTypeSim);
	sequenceContainer  simulateSequencesForParametricBootstrap(int seqLengthInSet, sequenceContainer& scSimulated, tree& trSampled, bool writeSeq=true,	bool useTheSame=true);
	void ancestralReconstructor();
	void ancestralReconstructorBasedOnJoint();
	Vdouble getRatesVector(){return _rates;};

	// co evol functions
	void findCoEvolvingSites(const int numberOfSequences2simulateForCoEvol);
	MDOUBLE computeCorrelationBetweenVis(const VVVdouble & VIpos_i, const VVVdouble & VIpos_j);

	MDOUBLE computeDistanceFromRootForRecent(tree& tr);	// 
	MDOUBLE computeDistanceNearestOTUforRecent(tree& tr);	// 
	//void bBLEMwithSimpleSpBeforeFullOptimization(tree& tr);
	void bBLEMwithSimpleSpBeforeFullOptimization(tree& tr, const sequenceContainer& sc, stochasticProcess* spSimple,													   
		stochasticProcess* sp,	
		const vector<vector<stochasticProcess*> >& spVVec,const distribution * gainDist, const distribution * lossDist,													   
		unObservableData *unObservableData_p); 

	void updateSetLofMissingData();
	void multipleAllBranchesByFactorAtStart(MDOUBLE epsilonOptimization);
	void multipleAllBranchesByFactorAtStartByMaxParsimonyCost(int costOfTreeMP);
	void RemoveSeqWithUnknownForSelectedSiteForCorrelation(sequenceContainer&  sc, tree& tr);




private:
	stochasticProcess *_sp; 
	vector<vector<stochasticProcess*> > _spVVec; //save stochasticProcess for each category
	stochasticProcess *_spSimple;
	Vdouble _freq;

	VVVdouble _postProbPerSpPerCatPerPos; // the posterior probability for each stochastic process for each rate Cat for each site

	distribution* _gainDist;
	distribution* _lossDist;
    tree _tr;
	tree _trOrig;	// used for diff(Branch length comparisons)
	tree _trGain;
	tree _trLoss;

	MDOUBLE _gainExp;
	MDOUBLE _lossExp;

	MDOUBLE _meanGain;
	MDOUBLE _meanLoss;
	MDOUBLE _medianGain;
	MDOUBLE _medianLoss;

	sequenceContainer _sc;
	sequenceContainer _scUniqPatterns;	// to contain a non-redundant set of patterns with _weights
	sequenceContainer _scWithFullLength;	//
	sequenceContainer _scFilterMissingData;	//

	vector<int> _alphVecDist;	// number of each letter
	
	//sequenceContainer _scZero;	
	//MDOUBLE _logLforMissingData;
	//MDOUBLE* _plogLforMissingData;
	//Vdouble _LforMissingDataPerCat;		// used foreach rate category
	//Vdouble* _pLforMissingDataPerCat;
	unObservableData*  _unObservableData_p;
	Vdouble* _weightsUniqPatterns;

	MDOUBLE _logL;
	MDOUBLE _distanceFromRootForRecent;
	MDOUBLE _distanceFromNearestOTUForRecent;

	sequence* _refSeq; // the reference sequence	
	VVVVdouble _jointProb_PosNodeXY; // store the information from computePosteriorOfChangeGivenTerminals
	VVVdouble _MPPerPos;	// The MP estimation of gain and loss events _MPPerPos[i][0][1] - gain events in i position
	int _CostOfTreeMP;
	VVVdouble _SMPerPos;	// The Stochastic mapping estimation of gain and loss events _SMPerPos[i][0][1] - gain events in i position
	VVVVdouble _MP_PosNodeXY;		// _MP_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb

	Vint	_occurPerPos;	// # 1
	Vint	_unknownPerPos;	// # ?

	Vdouble _gainPerPos;	// The Stochastic mapping estimation of gain and loss events _SMPerPos[i] - gain events in i position
	Vdouble _lossPerPos;	// The Stochastic mapping estimation of gain and loss events _SMPerPos[i] - loss events in i position
	Vdouble _lossMPPerPos; // Maximum Parsimony
	Vdouble _gainMPPerPos;

	Vdouble _gainPerPosCorr;	// either_SMPerPos[i], or _MPPerPos[i]
	Vdouble _lossPerPosCorr;	

	Vdouble _rates;// the rates themselves
	Vdouble _Lrate;// the log likelihood of each position
	VVdouble _postProbPerCatPerPos; // the posterior probability for each category and each site
	Vdouble _normalizedRates; // the rates when their ave = 0 and std = 1.
	MDOUBLE _ave; // the average over all rates.
	MDOUBLE _std; // the std over all rates.
	Vdouble _BayesianSTD;// the std of the Bayesian rates
	Vdouble _BayesianLowerBound;// lower bound of rate in Bayesian inference
	Vdouble _BayesianUpperBound;// upper bound of rate in Bayesian inference
	MDOUBLE _alphaConf; // the alpha confidence interval of Bayesian rates (set to 0.5). interval - rates that are in the 95% area under the curve.
	
	VVVVdouble _expChanges_PosNodeXY;		// expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb
	VVVVdouble _expChanges_PosNodeXYSampledData;		// expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb

	// correlation vectors
	VVVdouble _correlationsPerSitePerPosVec;
	VVVdouble _correlationsPerSitePerPosVecSampledData;
	vector<vector<bool> > _isComputePairWithRateAboveNim;  // not dependent on correlation type
	Vint _selectedSites;	// either all or selected sited (e.g., test correlation with specific traits)
	Vint _evolvingSites;	// sub-set of all sites in the sequence (e.g., with >=2 Event By MP) e.g., from seqLen = 5 _evolvingSites=[0,1,4]
	Vint _numOfGapsTillSite; // sub-set of all sites in the sequence (e.g., with >=2 Event By MP), e.g., _numOfGapsTillSite=[0,0,2]

	sequenceContainer _scEvolvingSites;

	map<int, map<int, map<string,  map<string, MDOUBLE > > > > _correlationsData; // _correlationsData["i"]["j"]["type"]["R" / "pVal" / "qVal" / "Nmin"]
};

#endif // ___GAIN_LOSS_
