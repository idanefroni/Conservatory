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


#ifndef ___GAINLOSS_UTILS__
#define ___GAINLOSS_UTILS__

#include "definitions.h"
#include "gainLossAlphabet.h"
#include "gammaDistribution.h"
#include "gammaDistributionFixedCategories.h"
#include "GamMixtureOptimizer.h"
#include "generalGammaDistributionPlusInvariant.h"
#include "logFile.h"
#include "matrixUtils.h"
#include "mixtureDistribution.h"
#include "someUtil.h"
#include "tree.h"
#include "treeIt.h"
#include "evaluateCharacterFreq.h"
#include "trivialAccelerator.h"
#include <math.h>

const string PROG_INFO = static_cast<string>("Version: gainLoss.VR01.266 - last updated 14.10.2013");
const MDOUBLE MINIMUM_PROB_PARAM = static_cast<MDOUBLE>(0.001);	
const MDOUBLE MAXIMUM_PROB_PARAM = static_cast<MDOUBLE>(0.999);
const MDOUBLE MINIMUM_FREQ_PARAM = static_cast<MDOUBLE>(0.001);	//0.05
const MDOUBLE MAXIMUM_FREQ_PARAM = static_cast<MDOUBLE>(0.999);	//0.95 
const MDOUBLE MINIMUM_GAIN_PARAM = static_cast<MDOUBLE>(0.0);	//0.01	
const MDOUBLE MAXIMUM_GAIN_PARAM = static_cast<MDOUBLE>(5.0);	
const MDOUBLE MINIMUM_LOSS_PARAM = static_cast<MDOUBLE>(0.01);	
const MDOUBLE MAXIMUM_LOSS_PARAM = static_cast<MDOUBLE>(10.0);

const MDOUBLE MINMUM_GAIN_LOSS_RATIO_PARAM = static_cast<MDOUBLE>(0.01);
const MDOUBLE MAXIMUM_GAIN_LOSS_RATIO_PARAM = static_cast<MDOUBLE>(100.0);

const int PRECISION = static_cast<int>(4);	// Used for print-outs
const int LOW_PRECISION = static_cast<int>(2);	// Used for print-outs, AncestralRec

void printTree (tree &tr, string treeFile);
void printTree (tree &tr,ostream &out);
void printTree (tree &tr);

void printTreeWithValuesAsBP(ostream &out, tree &tr, Vstring values, VVVdouble *probs=NULL ,bool printGains=true) ;
void printTreeWithValuesAsBP(ostream &out, const tree::nodeP &myNode, Vstring values, VVVdouble *probs=NULL ,bool printGains=true) ;

void printTreeStatesAsBPValues(ostream &out, Vint &states, tree &tr, VVVdouble *probs=NULL ,bool printGains=true) ;
void printTreeStatesAsBPValues(ostream &out, Vint &states, const tree::nodeP &myNode, VVVdouble *probs=NULL ,bool printGains=true) ;

void printTreeStatesAsBPValues(ostream &out, Vdouble &states, tree &tr, 
							   VVVdouble *probs=NULL ,bool printGains=true) ;
void printTreeStatesAsBPValues(ostream &out, Vdouble &states, const tree::nodeP &myNode, 
							   VVVdouble *probs=NULL ,bool printGains=true) ;

// --->> into somaUtils
//int fromIndex2gainIndex(const int i, const int gainCategories, const int lossCategories);
//int fromIndex2lossIndex(const int i, const int gainCategories, const int lossCategories);

MDOUBLE factorial (MDOUBLE num);

void printHelp();
void printProgramInfo();

bool isAlphaOptimization(distribution* dist);
bool isBetaOptimization(distribution* dist);
bool isMixOptimization(distribution* dist);
bool isInvariantOptimization(distribution* dist, bool onlyForPrintVal=false);
bool isThetaOptimization();

MDOUBLE getRateAlpha(distribution* dist);
MDOUBLE getRateBeta(distribution* dist);
//MDOUBLE getInvProbability(distribution* dist);

void setRateAlpha(distribution* dist, MDOUBLE paramAlpha);
void setRateBeta(distribution* dist, MDOUBLE paramBeta);

void updateGainAlpha(MDOUBLE param,
					 vector<vector<stochasticProcess*> >& spVVec,
					 distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateGainBeta(MDOUBLE param, 
					vector<vector<stochasticProcess*> >& spVVec,
					distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateGainProbInvariant(MDOUBLE param, distribution* gainDist);

void updateLossAlpha(MDOUBLE param, 
					 vector<vector<stochasticProcess*> >& spVVec,
					 distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateLossBeta(MDOUBLE param, 
					vector<vector<stochasticProcess*> >& spVVec,
					distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateLossProbInvariant(MDOUBLE param, distribution* lossDist);

void updateRateAlpha(MDOUBLE param, 
					 vector<vector<stochasticProcess*> >& spVVec,
					 distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateRateProbInvariant(MDOUBLE param, 
							 vector<vector<stochasticProcess*> >& spVVec,
							 distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);
void updateTheta(MDOUBLE param, 
				 vector<vector<stochasticProcess*> >& spVVec,
				 distribution * gainDist, distribution * lossDist, bool isNormalizeQ=true);

void cloneSpVVec(vector<vector<stochasticProcess*> >& spVVec, vector<vector<stochasticProcess*> >& neWspVVec);
void deleteSpVVec(vector<vector<stochasticProcess*> >* spVVec_p);

void clearVVVV(VVVVdouble& vetor);	
void clearVVV(VVVdouble& vetor);	

void resizeVVVV(int dim1, int dim2, int dim3, int dim4,  VVVVdouble& vetor);
void resizeVVV(int dim1, int dim2, int dim3, VVVdouble& vetor);

//MDOUBLE getDistance2ROOT(const tree::nodeP &myNode);
//MDOUBLE getMinimalDistance2OTU(const tree::nodeP &myNode); // Only for binary trees

//void startZeroSequenceContainer(const sequenceContainer &sc, sequenceContainer &scZero, gainLossAlphabet &alph);
void fillVnames(Vstring& Vnames,const tree& tr);
void P11forgain(ostream& out=cout)  ;

MDOUBLE normalizeQ(vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist);
MDOUBLE sumPijQijVec(vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist);
void normVec(const MDOUBLE scale, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist);
MDOUBLE normalizeQ(stochasticProcess*  sp);

MDOUBLE computeExpectationOfStationaryFrequency(distribution* gainDist, distribution* lossDist);
MDOUBLE computeExpectationOfGainLossRatio(distribution* gainDist, distribution* lossDist);
MDOUBLE computeExpOfGainByExpOfLossRatio(distribution* gainDist, distribution* lossDist);
MDOUBLE rateExpectation(distribution* dist);

void printMixtureParams(stochasticProcess*  sp);
stochasticProcess* startStochasticProcessSimpleGamma(MDOUBLE init_gain, MDOUBLE init_loss, Vdouble& freq, int numberOfRateCategories=4);

void readIntegersFromFileIntoVector(Vint& intVector, const int maxAllowed, const int minAllowed, string* inFile=NULL,Vint* evolvingSites=NULL);

void FlatTree(tree& trForSM , MDOUBLE defaultBranchLength=0.3);

void computeRateValPerPos(VVVVdouble& expChanges_PosNodeXY, VVVdouble& map_PosXY);

MDOUBLE computeNminRforCorrelWithGainAndLoss(MDOUBLE gainVal, MDOUBLE lossVal);


#endif

