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


#ifndef ___computeCountsGL___GL
#define ___computeCountsGL___GL

#include "definitions.h"
#include "replacementModel.h"
#include "gainLoss.h"

/********************************************************************************************
rate4siteGL
*********************************************************************************************/
class computeCountsGL{
public:
	explicit computeCountsGL(sequenceContainer& sc, tree& tr, stochasticProcess* sp,   string& outDir, VVdouble& LpostPerCat, MDOUBLE distanceFromNearestOTUForRecent, bool isSilent =false);
	explicit computeCountsGL(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution* gainDist, distribution* lossDist, string& outDir, VVVdouble& LpostPerSpPerCat, MDOUBLE distanceFromNearestOTUForRecent, bool isSilent=false);
	virtual ~computeCountsGL() ;

	computeCountsGL(const computeCountsGL& other) {*this = other;}	
	computeCountsGL& operator=(const computeCountsGL &other);
	void run();
	void computePosteriorOfChangeGivenTerminalsPerCat();
	void computePosteriorOfChangeGivenTerminalsPerSpPerCat();

	void printProbExp();
	void printProbabilityPerPosPerBranch();
	void printProbExpPerPosPerBranch(MDOUBLE probCutOff =0.0,MDOUBLE countsCutOff= 0.2);
	void printExpPerPosPerBranchMatrix(const int from, const int to);

	void printProbExpPerPosPerBranchFewCutOffs(MDOUBLE probCutOff);

	void produceExpectationPerBranch();
	void printExpectationPerBranch();
	void updateTreeByGainLossExpectationPerBranch(tree& tr, int from, int to);

	void printTreesWithExpectationValuesAsBP();
	void printTreesWithProbabilityValuesAsBP();

	//void computedCorrelations(const Vint& selectedPositions, const bool isNormalizeForBranch = false);
	//void printComputedCorrelations(const Vint& selectedPositions, const bool isNormalizeForBranch = false, const bool correlationForZscore = false);
	////void computeMeanAndSdPerBranch(Vdouble& meanEventsPerBranch01,	Vdouble& meanEventsPerBranch10,	Vdouble& sdEventsPerBranch01,Vdouble& sdEventsPerBranch10);
	//void fillMapValPerPosPerBranch(VVdouble& expEventsPerPosPerBranch,const int from, const int to, VVVVdouble& map_PosNodeXY
	//	,const bool isNormalizeForBranch = true, MDOUBLE* cutOff_p =NULL);
	//void fillCorrPerSelectedSites(Vdouble& correlationPerPos,VVdouble& expEventsPerPosPerBranch,const int selectedSite, const bool isPearson=true);


	Vdouble get_expV01(){return _expV01;};
	Vdouble get_expV10(){return _expV10;};
	VVVdouble get_expV(){return _expV;};

	Vdouble get_probV01(){return _probV01;};
	Vdouble get_probV10(){return _probV10;};
	VVVdouble get_probV(){return _probV;};

	VVVVdouble getExpChanges(){return _expChanges_PosNodeXY;};		// expChanges_PosNodeXY[pos][nodeID][x][y]
	VVVVdouble getProbChanges(){return _probChanges_PosNodeXY;};	// probChangesForBranch[pos][nodeID][x][y]
	VVVVdouble getJointProb(){return _jointProb_PosNodeXY;};		// _jointProb_PosNodeXY[pos][nodeID][x][y]

	
	//VVdouble getPerPosPerBranch01(){return _expPerPosPerBranch01;};		
	//VVdouble getPerPosPerBranch10(){return _expPerPosPerBranch10;};		
	//VVdouble getPerPosPerBranch(){return _expPerPosPerBranch;};	// vector of both, concatenated		

	//VVdouble getcorrelationPerSitePerPosGainGainSpearman(){return _correlationPerSitePerPosGainGainSpearman;};		
	//VVdouble getcorrelationPerSitePerPosLossLossSpearman(){return _correlationPerSitePerPosLossLossSpearman;};		
	//VVdouble getcorrelationPerSitePerPosBothSpearman(){return _correlationPerSitePerPosBothSpearman;};		

	//VVdouble getcorrelationPerSitePerPosGainGainPearson(){return _correlationPerSitePerPosGainGainPearson;};		
	//VVdouble getcorrelationPerSitePerPosLossLossPearson(){return _correlationPerSitePerPosLossLossPearson;};		
	//VVdouble getcorrelationPerSitePerPosBothPearson(){return _correlationPerSitePerPosBothPearson;};		



protected:
//func
	void printGainLossProbabilityPerPosPerBranch(int pos, MDOUBLE probCutOff, VVVdouble& probChanges, ostream& out, ostream& outCount);
	void printGainLossExpectationPerBranch(VVVdouble& expectChanges, ostream& out);

	void printGainLossProbExpPerPosPerBranch(int pos, MDOUBLE probCutOff, MDOUBLE countCutOff, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outCount);
	void printGainLossProbExpPerPosPerBranchFewCutOffs(int pos, MDOUBLE probCutOff, 
		MDOUBLE countCutOffLow,MDOUBLE countCutOffIncrem, MDOUBLE countCutOffHigh, VVVdouble& probChanges, VVVdouble& expChanges, ostream& out, ostream& outSum);


protected:
//members
	stochasticProcess *_sp;
	int _alphabetSize;
	tree _tr;
	sequenceContainer _sc;

	vector<vector<stochasticProcess*> > _spVVec; //save stochasticProcess for each category
	distribution* _gainDist;
	distribution* _lossDist;

	sequence* _refSeq; // the reference sequence
	string _outDir;
	bool _isSilent;

	VVdouble _postProbPerCatPerPos; // the posterior probability for each position for each rate category
	VVVdouble _postProbPerSpPerCatPerPos; // _LpostPerSpPerCat[sp][rateCat][pos]
	MDOUBLE _distanceFromNearestOTUForRecent;

	Vdouble _expV01;
	Vdouble _expV10;
	VVVdouble _expV;

	Vdouble _probV01;
	Vdouble _probV10;
	VVVdouble _probV;

	//VVVVdouble _posteriorsGivenTerminals;	// posteriorsGivenTerminals[pos][nodeID][x][y]
	VVVVdouble _probChanges_PosNodeXY;		// probChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations
	VVVVdouble _expChanges_PosNodeXY;		// expChanges_PosNodeXY[pos][nodeID][fatherState][sonState] - after simulations and postProb
	VVVdouble _expChanges_NodeXY;			// Summed from _expChanges_PosNodeXY - to expChanges_NodeXY[nodeID][fatherState][sonState]
	VVVVdouble _jointProb_PosNodeXY;		// probJoint_PosNodeXY[pos][nodeID][fatherState][sonState] - after computePosteriorOfChangeGivenTerminals

	//// required for correlation analysis
	//VVdouble _expPerPosPerBranch01;
	//VVdouble _expPerPosPerBranch10;
	//VVdouble _expPerPosPerBranch;
	//// correlation vectors
	//VVdouble _correlationPerSitePerPosGainGainSpearman;
	//VVdouble _correlationPerSitePerPosLossLossSpearman;
	//VVdouble _correlationPerSitePerPosBothSpearman;

	//VVdouble _correlationPerSitePerPosGainGainPearson;
	//VVdouble _correlationPerSitePerPosLossLossPearson;
	//VVdouble _correlationPerSitePerPosBothPearson;

};

#endif





















