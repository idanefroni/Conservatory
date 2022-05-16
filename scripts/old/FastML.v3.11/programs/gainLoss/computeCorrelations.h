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


#ifndef ___computeCorrelations___GL
#define ___computeCorrelations___GL

#include "definitions.h"
#include "replacementModel.h"
#include "Parameters.h"
#include "gainLoss.h"
#include "extremeValDistribution.h"

/********************************************************************************************
rate4siteGL
*********************************************************************************************/
class computeCorrelations{
public:
	explicit computeCorrelations(tree& tr,  string& outDir, VVVVdouble* expChanges_PosNodeXY, VVVVdouble* expChanges_PosNodeXY_B=NULL);
	virtual ~computeCorrelations() ;

	computeCorrelations(const computeCorrelations& other) {*this = other;}	
	computeCorrelations& operator=(const computeCorrelations &other);

	void runComputeCorrelations(const Vint& selectedPositions, const Vint& numOfGapsTillSite, const bool isNormalizeForBranch = false);
	void printComputedCorrelations(const Vint& selectedPositions,const Vint& evolvingSites, const bool isNormalizeForBranch = false, const bool correlationForZscore = false, VVVdouble* correlationsVec=NULL, string* valType=NULL);
	//void computeMeanAndSdPerBranch(Vdouble& meanEventsPerBranch01,	Vdouble& meanEventsPerBranch10,	Vdouble& sdEventsPerBranch01,Vdouble& sdEventsPerBranch10);
	void fillMapValPerPosPerBranch(VVdouble& expEventsPerPosPerBranch,const string type, VVVVdouble& expChanges_PosNodeXY,const bool isNormalizeForBranch = true, MDOUBLE* cutOff_p =NULL);
	void fillCorrPerSelectedSites(Vdouble& correlationPerPos,VVdouble& expEventsPerPosPerBranch,VVdouble& expEventsPerPosPerBranch_B,const int selectedSite, const bool isPearson=true);
	void sumExpectationPerBranch(VVVVdouble& expChanges_PosNodeXY, VVVdouble& map_NodeXY);
	MDOUBLE computeNminPerPair(const int site_A, const int site_B, const int typeIndex, const VVVdouble&  exp_PosXY);


	void computedCorrelationsRankBasedOnSimulatedData(const Vint& selectedPositions, VVVdouble& correlationPerSitePerPos, VVVdouble& correlationPerSitePerPos_B, VVVdouble& correlationPerSitePerPos_Pval);
	void computedCorrelationsPValBasedOnSimulatedDataCoMap(VVVdouble& correlationPerSitePerPosReal,VVVVdouble& expChanges_PosXYReal, VVVdouble& correlationPerSitePerPos_Pval);
	int computedCorrelationsPValBasedOnSimulatedDataCoMapBins(VVVdouble& correlationPerSitePerPosReal,vector<vector<bool> >& isComputePairWithRateAboveNim,VVVVdouble& expChanges_PosXYReal, VVVdouble& correlationPerSitePerPos_Pval
		,map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData,  Vdouble& rate4siteReal, Vint& selectedSites, Vint& numOfGapsTillSite, Vint& evolvingSites, bool isLastIteration);
	void printComputedCorrelationsData(const bool isNormalizeForBranch, const bool correlationForZscore
		,map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData, Vdouble& T_BH, bool isPairsAboveBH = false);
	void printCorrelationsFrequencies(Vdouble& correlationsVecSorted, ofstream* simCorrelStream=NULL);

	
	int produceSortedVectorsOfCorrelationsBinedByRate(MDOUBLE medianNminOfRealData, ofstream* simCorrelStream);	

	void produceSortedVectorsOfAllCorrelations(Vdouble& rate4siteSim);
	
	VVVdouble pVals2qVals(VVVdouble& correlationsVec,map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData
		, vector<vector<bool> >& isComputePairWithRateAboveNim, Vdouble& T_BH, Vint& selectedSites, Vint& evolvingSites);

	void produceSymeticMatrix(VVVdouble& correlationPerSitePerPos_Pval, bool isMin=true);
	void produceSortedVectorsOfAllCorrelations(const VVVdouble& correlationPerSitePerPos, Vdouble& pairWiseCorrelations, Vdouble& NminForPairsInPairWiseCorrelations);

	VVVdouble getcorrelationPerSitePerPosVec(){return _correlationsPerSitePerPosVec;};


protected:
//members
	int _alphabetSize;
	tree _tr;
	//sequenceContainer _sc;

	sequence* _refSeq; // the reference sequence
	string _outDir;
	bool _isSilent;


	VVVVdouble _expChanges_PosNodeXY;		// Input, expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb
	VVVdouble _expChanges_NodeXY;			// Summed from _expChanges_PosNodeXY - to expChanges_NodeXY[nodeID][fatherState][sonState]
	VVVdouble _exp_PosXY;			// Summed from _expChanges_PosNodeXY - to expChanges_PosXY[Pos][fatherState][sonState]

	bool _isTwoSetsOfInputForCorrelation;	// when B is given
	VVVVdouble _expChanges_PosNodeXY_B;		// Input B (optional), expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb
	VVVdouble _expChanges_NodeXY_B;			// Summed from _expChanges_PosNodeXY - to expChanges_NodeXY[nodeID][fatherState][sonState]

	// V required for correlation analysis
	VVVdouble _expPerPosPerBranchVec;	//  expChanges_PosNodeXY[type][pos][nodeID], for specific type of event (from, to), may be adjusted for branch expectation
	VVVdouble _expPerPosPerBranchVec_B;

	// correlation vectors
	VVVdouble _correlationsPerSitePerPosVec;
	vector<bool> _isPearson; // [true, false]
	vector<string> _EventTypes; // ['gain', 'loss', 'both']
	map<string, int> _EventTypesMap;
	map<string, map<string, int> > _EventTypesFromTo;

	//vector< vector< map<string, MDOUBLE> > > _pairWiseCorrelationsAndNminSim; // pairWiseCorrelationsAndNmin[corrIndex][pairIndex][CorOrNmin][val], if CorOrNmin=0, val=correlation, if =1, val=Nmin
	VVVdouble _pairWiseCorrelationsAndNminSim; // pairWiseCorrelationsAndNmin[corrIndex][0/1][pairIndex][val], if CorOrNmin=0, val=correlation, if =1, val=Nmin
	VVdouble _corrVector;				
	VVdouble _NminSortedSim;			// _NminSortedSim[CorType][], vector of all Nmins = Rates
	vector<vector<extremeValDistribution > > _extremeValDistributions; // _NminSortedSim[CorType][], vector of all distributions, per bin

	VVVdouble _correlationsSubSets;		// to be filled by produceSortedVectorsOfCorrelationsBinedByRate
	VVdouble _correlationSubSetsNminLimitValues; // to be filled by produceSortedVectorsOfCorrelationsBinedByRate

	Vint _selectedSites;

	int _numOfSamplesInLowRateFirstBin;  // thus, the lowest p-value for correlations with low rate (below simulations) is limited
};

#endif




