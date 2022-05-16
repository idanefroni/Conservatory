
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

#include "computePosteriorExpectationOfChange.h"
#include "gainLoss.h"
#include "gainLossOptimizer.h"
#include "gainLossOptions.h"
#include "gainLossUtils.h"
#include "gammaDistributionFixedCategories.h"
#include "gammaDistributionPlusInvariant.h"
#include "mixtureDistribution.h"
#include "simulateTree.h"
#include "someUtil.h"
#include "phylipFormat.h"
#include "maseFormat.h"
#include "fastaFormat.h"
#include "clustalFormat.h"
#include "rate4siteGL.h"
#include "gainLoss4site.h"
#include "computeCountsGL.h"
#include "computeCorrelations.h"

#include "simulateOnePos.h"
#include "Parameters.h"
#include "sankoffReconstructGL.h"
#include "bblLS.h"
//#include "branchScaleTree.h"
#include "bblEMfixRoot.h"
#include "bblEM.h"

#include <cstring>


using namespace std;
/********************************************************************************************
gainLoss TOC (group of functions by order of appearance):
	-constructor+destructor
	-initialize
	-run
	-start(basics):			SequenceContainer, StochasticProcess(+Generic,+Vec), EvolTreeTopology,initMissingDataInfo
	-optimize:				startOptimizations, optimizationsManyStarts(+NoVec, +VV), initParamsAtRandPoints(+SPvv)
	-start(computations):	
	-prints
	-Mixture
	-simulate
	-Old function, now inside other classes
*********************************************************************************************/
gainLoss::gainLoss(): _sp(NULL),_unObservableData_p(NULL),_lossDist(NULL), _gainDist(NULL), _refSeq(NULL), _weightsUniqPatterns(NULL)
{
	_weightsUniqPatterns = gainLossOptions::_weights; // since - no weights are used over positions, it is NULL
	_logL = 1;
	//_maxNumberOfSpeciesForFullOptimization = 200;
	//_maxSequenceLengthForFullOptimization = 20000;
	//_maxSpeciesNumSequenceLengthMultipForFullOptimization = 500000;
}
/********************************************************************************************/
gainLoss::~gainLoss() {
	if(gainLossOptions::_gainLossDist){
		for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
			for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){
				stochasticProcess* sp2delete = _spVVec[gainCategor][lossCategor];
				delete sp2delete;
			}
		}
		delete _gainDist;
		delete _lossDist;
	}
	else
		if (_sp) delete _sp;
	
	if(_unObservableData_p) delete _unObservableData_p;
	if(_weightsUniqPatterns) delete _weightsUniqPatterns;
	if(_spSimple) delete _spSimple;
}
/********************************************************************************************
*********************************************************************************************/
void gainLoss::initialize(bool isComputeLikelihood) 
{	
	printProgramInfo();
	printOptionParameters();
	if(gainLossOptions::_seqFile!=""){
		startSequenceContainer();
		fillReferenceSequence();
	}
	countOccurPerPos();
	
	if(gainLossOptions::_isRemovePositionsWithHighPercentOfMissingData)
		removePositionsWithHighPercentOfMissingData(0.5);

	if(gainLossOptions::_isSequenceUniqPattern)
		startSequenceContainerUniqPatterns();

	startStochasticProcess(gainLossOptions::_gainLossDist);
	MDOUBLE epsilon2add = 0.0;
	//if(_gainExp<1e-08)
		//epsilon2add = 1e-08;
	_spSimple =  startStochasticProcessSimpleGamma(_gainExp+epsilon2add,_lossExp,_freq); // simple initialization, based on empiricalCounting of '1' and '0'
	MDOUBLE norm_factor = normalizeQ(_spSimple);
	LOGnOUT(4,<<"Stochastic process 'simple' normalized with norm_factor="<<norm_factor<<endl);

	startEvolTreeTopology();
	
	if(gainLossOptions::_intersectTreeAndSeq){ // input tree and seq (not the same taxa) - intersect, write seq and tree
		intersectNamesInTreeAndSequenceContainer(_tr,_sc);
		LOGnOUT(4,<<"NumOfSeq= "<<_sc.numberOfSeqs()<<endl);
		LOGnOUT(4,<<"NumOfTaxa= "<<_tr.getLeavesNum()<<endl);
		bool isRemovePosNotWithinMinMax=true;
		int minNumOfOnes = Parameters::getInt("_minNumOfOnes");
		int minNumOfZeros = Parameters::getInt("_minNumOfZeros");
		//sequenceContainer&  sc, int minNumOfOnes, int  minNumOfZeros, bool isRemovePosNotWithinMinMax
		bool isReportRemovedPos = true;
		checkMinNumOfOnesOrZeros(_sc,minNumOfOnes,minNumOfZeros, isRemovePosNotWithinMinMax,  isReportRemovedPos);
		_scWithFullLength = _sc;
		_scUniqPatterns = _sc;
		_trOrig = _tr;
		string strSeqNum = gainLossOptions::_outDir + "//" + "seq" + ".fa";
		ofstream seq_out(strSeqNum.c_str());
		fastaFormat::  write(seq_out,_scWithFullLength);
		ofstream treeStream(gainLossOptions::_treeOutFile.c_str());
		_tr.output(treeStream);
		treeStream.close();
		return;
	}
	if(gainLossOptions::_seqFile!="")
		checkThatNamesInTreeAreSameAsNamesInSequenceContainer(_tr,_sc);
	if(gainLossOptions::_seqFile!="" && Parameters::getInt("_accountForMissingData") && (Parameters::getInt("_minNumOfOnes")>0 || Parameters::getInt("_minNumOfZeros")>0)){
		initializeUnObservableData();
	}
	if(gainLossOptions::_seqFile!="" && isComputeLikelihood){		
		printTreeLikelihoodAllPosAlphTheSame();	// update of _logL is done as well
	}
	if(Parameters::getInt("_isNormalizeAtStart")){
		bool isNormalizeBothQandTree = false;	// Under the assumption that the input tree was normalized, only need to start with Q
		normalizeQandTree(isComputeLikelihood, isNormalizeBothQandTree);
	}
	printTree(_tr);
	printModellValuesOfParams();

	if(gainLossOptions::_printSeq && _sc.seqLen() != _scWithFullLength.seqLen() ){
		string strSeqNum = gainLossOptions::_outDir + "//" + "seq.not.full.length.fa";
		ofstream seq_out(strSeqNum.c_str());
		fastaFormat::  write(seq_out,_sc);				// not full length
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::bBLEMwithSimpleSpBeforeFullOptimization(tree& tr, const sequenceContainer& sc, stochasticProcess* spSimple,													   
													   stochasticProcess* sp,	
													   const vector<vector<stochasticProcess*> >& spVVec,const distribution * gainDist, const distribution * lossDist,													   
													   unObservableData *unObservableData_p) 
{
	LOGnOUT(4,<<" *** Starting bbBLEMwithSimpleSpBeforeFullOptimization"<<endl);	
	//MDOUBLE oldLnoUnObservableDataCorrection = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,*spSimple,_weightsUniqPatterns,NULL);
	//MDOUBLE newLnoUnObservableDataCorrection;
	MDOUBLE oldL = _logL;
	tree oldTree = tr;
	MDOUBLE tollForPairwiseDist=0.01; // the BBL default, epsilon per branch (brent's value)
	MDOUBLE epsilonOptimizationBBLIter = max(0.1,abs(_logL)/10000);
	int maxNumOfIterationsBBL =50;
	//bool isFixedRoot = !gainLossOptions::_isReversible && !gainLossOptions::_isRootFreqEQstationary;
	//if(isFixedRoot){
	//	LOGnOUT(4,<<"*** Fix Root BBL-EM.  Likelihood="<<oldL<<endl);
	//	LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
	//	bblEMfixRoot bblEM1(_tr, _sc, *_spSimple, NULL, maxNumOfIterationsBBL , epsilonOptimizationBBLIter,tollForPairwiseDist
	//		,NULL,&oldLnoUnObservableDataCorrection);
	//	newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
	//}
	//else{
	LOGnOUT(4,<<"*** BBL-EM.  Likelihood="<<oldL<<endl);
	LOGnOUT(4,<<" BBL-EM:  tollForPairwiseDist="<<tollForPairwiseDist<<" and epsilonOptimizationBBLIter="<<epsilonOptimizationBBLIter<<endl);
	bblEM bblEM1(tr, sc, *spSimple, NULL, maxNumOfIterationsBBL , epsilonOptimizationBBLIter,tollForPairwiseDist
		,NULL,NULL); // last argument optional: &oldLnoUnObservableDataCorrection

	//newLnoUnObservableDataCorrection = bblEM1.getTreeLikelihood();
	//}
	if(_unObservableData_p){
		if(!gainLossOptions::_gainLossDist)
			_unObservableData_p->setLforMissingData(tr,sp);
		else
			_unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
	}
	if(!gainLossOptions::_gainLossDist)
		_logL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,*sp,_weightsUniqPatterns,unObservableData_p);
	else
		_logL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,_weightsUniqPatterns,unObservableData_p);
	if(_logL  < oldL){
		LOGnOUT(4,<<"Likelihood went down with simplified BBL-EM. from "<<oldL<<" to "<<_logL<<" Go back to old tree."<<endl);
		tr = oldTree;
		_logL = oldL;
	}
	else{
		LOGnOUT(4,<<"Likelihood after BBL-EM="<<_logL<<endl);
	}
	printTree(_tr);
	LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::initializeBranchLengthDiff() 
{				
	printProgramInfo();
	printOptionParameters();
	startSequenceContainer();
	fillReferenceSequence();	
	startStochasticProcess(gainLossOptions::_gainLossDist);	
	//startEvolTreeTopology();
	_tr= tree(gainLossOptions::_treeFile);
	checkThatNamesInTreeAreSameAsNamesInSequenceContainer(_tr,_sc);
	if(Parameters::getInt("_accountForMissingData")){
		_unObservableData_p = new unObservableData(_sc, _sp, gainLossAlphabet() ,Parameters::getInt("_minNumOfOnes"), Parameters::getInt("_minNumOfZeros"));
		LOGnOUT(4,<<"unObservableData object initialized with number of unObservable patterns= "<<_unObservableData_p->getNumOfUnObservablePatterns() <<endl);
		if(Parameters::getInt("_minNumOfOnes") >= _sc.numberOfSeqs())
			errorMsg::reportError("Error: number of seqs smaller than minNumOfOnes\n");
		updateSetLofMissingData();
	}
	printTreeLikelihoodAllPosAlphTheSame();
	//if(!gainLossOptions::_gainLossDist)
	//	_logL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*_sp,_weightsUniqPatterns,_unObservableData_p);
	//else{
	//	_logL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
	//}
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::initializeUnObservableData(){
	if(Parameters::getInt("_minNumOfOnes") >= _sc.numberOfSeqs())
		errorMsg::reportError("Error: number of seqs smaller than minNumOfOnes\n");
	if( (_sc.numberOfSeqs()>250) && (Parameters::getInt("_minNumOfOnes") >1) )
		LOGnOUT(4,<< "WARNING: There are more than 250 sequences. Using more than 1 unObseravable pattern will run to slow\n");

	_unObservableData_p = new unObservableData(_scWithFullLength, _sp, gainLossAlphabet() ,Parameters::getInt("_minNumOfOnes"), Parameters::getInt("_minNumOfZeros"));

	LOGnOUT(4,<<"unObservableData object initialized with number of unObservable patterns= "<<_unObservableData_p->getNumOfUnObservablePatterns() <<endl);
	if(!gainLossOptions::_gainLossDist)
		_unObservableData_p->setLforMissingData(_tr,_sp);
	else
		_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::run(){
	// Special options, partial runs
	if(gainLossOptions::_calculeBranchLegthDiffFactorFromInputTrees){	// if BBL is used for each branch - compare length before/after		
		LOGnOUT(4,<<"\n\n RUN type: calculeBranchLegthDiffFactorFromInputTrees and return \n\n"<<endl);	
		initialize();
		//initializeBranchLengthDiff();
		_trOrig= tree(gainLossOptions::_treeFileOrig);
		string branchLegthDiffFactor = gainLossOptions::_outDir + "//" + "branchLegthDiffFactor.txt"; 
		ofstream branchLegthDiffFactorStream(branchLegthDiffFactor.c_str());
		branchLegthDiffFactorStream.precision(PRECISION);
		computeBranchLegthDiffFactor(branchLegthDiffFactorStream);
		return;
	}
	if(gainLossOptions::_intersectTreeAndSeq){ 
		LOGnOUT(4,<<"\n\n RUN type: intersect Tree And Seq, write them and return \n\n"<<endl);	
		initialize();
		return;
	}
	if(gainLossOptions::_simulatePosteriorExpectationOfChange){	// Test the simulation based counting
		LOGnOUT(4,<<"\n\n RUN type: simulatePosteriorExpectationOfChange \n\n"<<endl);
		initialize();
		if(Parameters::getInt("_performOptimizations") ){ //&& ! Parameters::getInt("_keepUserGainLossRatio")
			startOptimizations();
		}
		if (gainLossOptions::_treeFile=="") 
			errorMsg::reportError("SimulatePosteriorExpectationOfChange require input tree ");
		startSimultePosteriorExpectationOfChange(gainLossOptions::_numberOfPositions2simulate,gainLossOptions::_numberOfSequences2simulate);
		return;
	}
	if(gainLossOptions::_printLikelihoodLandscape){
		LOGnOUT(4,<<"\n\n RUN type: printLikelihoodLandscape \n\n"<<endl);
		//printLikelihoodLandscape(_sp);	//UtilFunction (Ad-hoc)
		initialize();
		printLikelihoodLandscapeStatFreqRatioAndRootFreqRatio();	//UtilFunction (Ad-hoc)
		return;
	}
	if(gainLossOptions::_printP11forgain){
		LOGnOUT(4,<<"\n\n RUN type: _printP11forgain \n\n"<<endl);
		initialize();
		P11forgain();	//UtilFunction (Ad-hoc)
		return;
	}
	if(gainLossOptions::_isOnlyParsimony){
		initialize(false);
		startMaxParsimonyChange();
		//printTree(_tr); // already done in initialize
		return;        
	}
	bool Normalize =false;
	if(Normalize){
		bool isComputeL = false;
		initialize(isComputeL);
		bool isNormalizeBothQandTree = false;	// Under the assumption that the input tree was normalized, only need to start with Q
		normalizeQandTree(isComputeL, isNormalizeBothQandTree);		
		//printTree(_tr); // already done in initialize
		return;
	}

	//////////////////////////////////////////////////////////////////////////	
	////////////////// Start normal run
	initialize(gainLossOptions::_isComputeLikelihoodDuringInit);	
	printTree(_tr);	
	//printQ();
	if(gainLossOptions::_calculeMaxParsimonyChange || gainLossOptions::_isRemoveSimulatedPositionsBasedOnMP || gainLossOptions::_isCorrelationsBasedOnMaxParsimonyMapping){	// simulation based counting 
		startMaxParsimonyChange();
	}

	startOptimizations();
	if(Parameters::getInt("_isNormalizeQandTreeafterOpt")){
		normalizeQandTree();
	}
	printModellValuesOfParams();

	// intersect according to missing data
	if(gainLossOptions::_isRemoveSeqWithUnknownForLastSelectedSiteForCorrelation){
		LOGnOUT(4,<<"\n Remove Seq With Unknown For Last Position. \n Compute correlation for last position"<<endl);
		RemoveSeqWithUnknownForSelectedSiteForCorrelation(_sc,_tr);
	}

	if(gainLossOptions::_printTree){
		printTree(_tr);	// Two trees are printed. ("TheTree.INodes.ph" with internal nodes)
	}
	if(gainLossOptions::_printSeq){
		string strSeqNum = gainLossOptions::_outDir + "//" + "seq.fa";
		ofstream seq_out(strSeqNum.c_str());
		fastaFormat::  write(seq_out,_scWithFullLength);
	}
	if (gainLossOptions::_isComputeDistanceFromRootForRecent){
		LOGnOUT(4,<<"\n Estimate: distanceFromRootForRecent, distanceFromNearestOTUForRecent"<<endl);
		_distanceFromRootForRecent = computeDistanceFromRootForRecent(_tr);
		_distanceFromNearestOTUForRecent = computeDistanceNearestOTUforRecent(_tr);	
	}else{
		LOGnOUT(4,<<"\n WARN: distanceFromRootForRecent=1, distanceFromNearestOTUForRecent=0.000001 are not estimated"<<endl);
		_distanceFromRootForRecent = 1;
		_distanceFromNearestOTUForRecent = 0.000001;
	}
	if(gainLossOptions::_printPij_t){
		printPij_t(0.1);
		printQ();
	}
	if(gainLossOptions::_printLofPos){
		printLofPos();	
	}
	if(gainLossOptions::_printLofPosBothModels){
		//LOGnOUT(4,<<"_printLofPosBothModels not implemented in this version"<<endl);	
		printLofPosBothModels();	// Print Likelihood of each position for m0 and m1
	}	
	if(gainLossOptions::_calculateRate4site  && !gainLossOptions::_gainLossDist && gainLossOptions::_rateDistributionType!=gainLossOptions::UNIFORM){
		startRate4Site(_scWithFullLength,_tr,_sp,gainLossOptions::_outDir,_unObservableData_p);
	}
	if(gainLossOptions::_calculeGainLoss4site && gainLossOptions::_gainLossDist){
		startGainLoss4Site(_scWithFullLength,_tr,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_unObservableData_p);
	}
	
	
	// Note: fill  VVVVdouble _probChanges_PosNodeXY - Delete it after AncestralReconstruct?
	if(gainLossOptions::_calculePosteriorExpectationOfChange && !gainLossOptions::_isCorrelationsBasedOnMaxParsimonyMapping){
		startComputePosteriorExpectationOfChange();
	}
	if(gainLossOptions::_printComputedCorrelations){
		startComputeAmongSitesCorrelations();	
	}
	if(gainLossOptions::_performParametricBootstapCorrelation){
		startParametricBootstapCorrelation();
	}
	if(gainLossOptions::_printTree && gainLossOptions::_calculePosteriorExpectationOfChange && _SMPerPos.size()>0){ // to check if it was done
		string treeGain = gainLossOptions::_outDir + "//" + "TheTree.Gain.ph"; 
		printTree(_trGain, treeGain);	
		string treeLoss = gainLossOptions::_outDir + "//" + "TheTree.Loss.ph";
		printTree(_trLoss, treeLoss);	
	}
	if(gainLossOptions::_calculateAncestralReconstruct){
		//LOGnOUT(4,<<"_calculateAncestralReconstruct not implemented in this version"<<endl);
		ancestralReconstructorBasedOnJoint();
		//ancestralReconstructor();
	}	
	if(gainLossOptions::_calculeBranchLegthDiffFactor){	// if BBL is used for each branch - compare length before/after		
		//LOGnOUT(4,<<"_calculeBranchLegthDiffFactor not implemented in this version"<<endl);	
		string branchLegthDiffFactor = gainLossOptions::_outDir + "//" + "branchLegthDiffFactor.txt"; 
		ofstream branchLegthDiffFactorStream(branchLegthDiffFactor.c_str());
		branchLegthDiffFactorStream.precision(PRECISION);
		computeBranchLegthDiffFactor(branchLegthDiffFactorStream);
	}	
	if(gainLossOptions::_simulateSequences){					// Test the rate4site computation
		//LOGnOUT(4,<<"_simulateSequences not implemented in this version"<<endl);	
		startSimulateSequences(gainLossOptions::_numberOfSequences2simulate, _sc.seqLen());
	}
	if(gainLossOptions::_findCoEvolvingSitesOldNotWorking){
		//LOGnOUT(4,<<"_findCoEvolvingSitesOldNotWorking not implemented in this version"<<endl);	
		findCoEvolvingSites(gainLossOptions::_numberOfSequences2simulateForCoEvol);
	}	
}



// start(basics)
/********************************************************************************************
startSequenceContainer
*********************************************************************************************/
void gainLoss::startSequenceContainer(){
	LOGnOUT(4,<<"\n startSequenceContainer"<<endl);
	bool isCountUnknownChars = true; // move to options

	gainLossAlphabet alph;
	ifstream in(gainLossOptions::_seqFile.c_str());
	sequenceContainer original = recognizeFormat::read(in,&alph);
	original.changeGaps2MissingData();
	_sc = original;
	_scWithFullLength = original;
	_scUniqPatterns = original;
	if(Parameters::getInt("_accountForMissingData")){
		int minNumOfOnes = Parameters::getInt("_minNumOfOnes");
		int minNumOfZeros = Parameters::getInt("_minNumOfZeros");

		bool isRemovePosNotWithinMinMax=false;
		bool isReportRemovedPos=false;
		checkMinNumOfOnesOrZeros(_sc,minNumOfOnes,minNumOfZeros, isRemovePosNotWithinMinMax, isReportRemovedPos);		
	}
	if(gainLossOptions::_checkCoEvolWithUnionPAP_against_pos){
		produceUnionPAP_against_pos(_sc, gainLossOptions::_checkCoEvolWithUnionPAP_against_pos);
	}	



	int alphSize = alph.size();
	if(isCountUnknownChars)
		alphSize++;
	_alphVecDist.resize(alphSize);
	_alphVecDist = _sc.getAlphabetDistribution(isCountUnknownChars);
	LOGnOUT(4,<<"numberOfSeqs="<<_sc.numberOfSeqs()<<endl);	
	LOGnOUT(4,<<"seqLen="<<_sc.seqLen()<<endl);
	LOGnOUT(4,<<"Num of zeros="<<_alphVecDist[0]<<"\nNum of ones="<<_alphVecDist[1]<<endl);
	if(isCountUnknownChars)
		LOGnOUT(4,<<"Num of unKnowns="<<_alphVecDist[2]<<endl);

	//bool isOverRideDataSizeForOptimization = false;
	//if(!isOverRideDataSizeForOptimization 
	//	&& (_sc.numberOfSeqs()>maxNumberOfSpeciesForFullOptimization 
	//		|| _sc.seqLen()>maxSequenceLengthForFullOptimization 
	//		|| _sc.numberOfSeqs()* sc.seqLen()>maxSpeciesNumSequenceLengthMultipForFullOptimization) ){
	//			LOGnOUT(2,<<"WARN: optimization level is reduced with too large dataset.\n To overRide re-run with _isOverRideDataSizeForOptimization 0"<<_sc.numberOfSeqs()<<endl);
	//}

	//Parameters::updateParameter("_calculeBranchLegthDiffFactor","0"); // why is it here???

}

/********************************************************************************************
The likelihood correction, requires that unObservable patterns do not exist
*********************************************************************************************/
void gainLoss::produceUnionPAP_against_pos(sequenceContainer&  sc, int pos_for_union, bool is_ignore_last_pos){
	LOGnOUT(4,<<"produceUnionPAP_against_pos with "<<pos_for_union<<endl);
	int seq_length_to_union =  sc.seqLen();
	if(is_ignore_last_pos){
		seq_length_to_union--;
		LOGnOUT(4,<<"Ignore last pos for union. Modify positions "<<seq_length_to_union<<endl);
		seq_length_to_union--;
	}
	int pos_for_union_start_count_from_0 = pos_for_union-1;
	for (int pos = 0; pos < seq_length_to_union; ++pos){
		for (int seqID = 0; seqID < sc.numberOfSeqs(); ++seqID){
			if(sc[seqID][pos] == 1 || sc[seqID][pos_for_union_start_count_from_0] == 1  ){
				sc[seqID][pos] = 1;
			}
		}	
	}
}





/********************************************************************************************
The likelihood correction, requires that unObservable patterns do not exist
*********************************************************************************************/
void gainLoss::checkMinNumOfOnesOrZeros(sequenceContainer&  sc, int minNumOfOnes, int  minNumOfZeros, bool isRemovePosNotWithinMinMax, bool isReportRemovedPos){
	vector<int> posToRemove(sc.seqLen(),false);
	vector<int> _alphVecDist = sc.getAlphabetDistribution();
	int numOfPosBelowMinNumOfOnes = 0;
	int numOfPosBelowMinNumOfZeros = 0;
	for (int pos = 0; pos < sc.seqLen(); ++pos){		
		Vint alphVecPerPos = sc.getAlphabetDistribution(pos);
		if(alphVecPerPos[1]< minNumOfOnes){
			if(isRemovePosNotWithinMinMax || gainLossOptions::_intersectTreeAndSeq){
				posToRemove[pos] = true;
				numOfPosBelowMinNumOfOnes++;
				if(isReportRemovedPos || gainLossOptions::_intersectTreeAndSeq)
					LOGnOUT(4,<<"Belove minOnes, Remove pos="<<pos+1<<endl);
			}
			else{
				LOGnOUT(4,<<"! WARN: Illegal minNumOfOnes found in pos="<<pos+1<<" Reset to minNumOfOnes="<<alphVecPerPos[1]<<endl);				
				Parameters::updateParameter("_minNumOfOnes",int2string(alphVecPerPos[1]).c_str());
				minNumOfOnes = alphVecPerPos[1];
			}
		}
		if(alphVecPerPos[0]< minNumOfZeros){
			if(isRemovePosNotWithinMinMax || gainLossOptions::_intersectTreeAndSeq){
				posToRemove[pos] = true;
				numOfPosBelowMinNumOfZeros++;
				if(isReportRemovedPos || gainLossOptions::_intersectTreeAndSeq)					
					LOGnOUT(4,<<"Belove minZeros, Remove pos="<<pos+1<<endl);
			}
			else{
				LOGnOUT(4,<<"! WARN: Illegal minNumOfZeros found in pos="<<pos+1<<" Reset to minNumOfZeros="<<alphVecPerPos[0]<<endl);				
				Parameters::updateParameter("_minNumOfZeros",int2string(alphVecPerPos[0]).c_str());
				minNumOfZeros = alphVecPerPos[0];
			}
		}
	}
	if(minNumOfOnes==0 && minNumOfZeros==0){
		LOGnOUT(4,<<"!!! WARN: both minNumOfOnes and minNumOfZeros=0. Thus, no need perform likelihood correction (accountForMissingData)"<<endl);
		Parameters::updateParameter("_accountForMissingData","0");
	}
	if(numOfPosBelowMinNumOfOnes>0)
		LOGnOUT(4,<<"WARN: removed "<<numOfPosBelowMinNumOfOnes<<" positions below minNumOfOnes="<<minNumOfOnes<<endl);
	if(numOfPosBelowMinNumOfZeros>0)
		LOGnOUT(4,<<"WARN: removed "<<numOfPosBelowMinNumOfZeros<<" positions below minNumOfZeros="<<minNumOfZeros<<endl);
	sc.removePositions(posToRemove);
}


/********************************************************************************************
The likelihood correction, requires that unObservable patterns do not exist
*********************************************************************************************/
//void gainLoss::removePositionsBelowNmin(sequenceContainer&  sc, int minNumOfOnes, int  MinVal, bool isRemovePosNotWithinMinMax, bool isReportRemovedPos){
//	vector<int> posToRemove(sc.seqLen(),false);
//	vector<int> _alphVecDist = sc.getAlphabetDistribution();	
//	for (int pos = 0; pos < sc.seqLen(); ++pos){		
//		Vint alphVecPerPos = sc.getAlphabetDistribution(pos);
//		if(alphVecPerPos[1]< minNumOfOnes){
//			if(isRemovePosNotWithinMinMax || gainLossOptions::_intersectTreeAndSeq){
//				posToRemove[pos] = true;
//				if(isReportRemovedPos)
//				 LOGnOUT(4,<<"Belove minOnes, Remove pos="<<pos+1<<endl);
//			}
//			else{
//				LOGnOUT(4,<<"! WARN: Illegal minNumOfOnes found in pos="<<pos+1<<" Reset to minNumOfOnes="<<alphVecPerPos[1]<<endl);				
//				Parameters::updateParameter("_minNumOfOnes",int2string(alphVecPerPos[1]).c_str());
//			}
//		}
//		if(alphVecPerPos[0]< minNumOfZeros){
//			if(isRemovePosNotWithinMinMax || gainLossOptions::_intersectTreeAndSeq){
//				if(isReportRemovedPos)
//					posToRemove[pos] = true;
//				LOGnOUT(4,<<"Belove minZeros, Remove pos="<<pos+1<<endl);
//			}
//			else{
//				LOGnOUT(4,<<"! WARN: Illegal minNumOfZeros found in pos="<<pos+1<<" Reset to minNumOfZeros="<<alphVecPerPos[0]<<endl);				
//				Parameters::updateParameter("_minNumOfZeros",int2string(alphVecPerPos[0]).c_str());
//			}
//		}
//	}
//	if(minNumOfOnes==0 && minNumOfZeros==0){
//		LOGnOUT(4,<<"!!! WARN: both minNumOfOnes and minNumOfZeros=0. Thus, no need perform likelihood correction (accountForMissingData)"<<endl);
//		Parameters::updateParameter("_accountForMissingData","0");
//	}
//	sc.removePositions(posToRemove);
//}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::countOccurPerPos(){
	string occur = gainLossOptions::_outDir + "//" + "occurPerPos.txt"; 
	ofstream occurStream(occur.c_str());
	occurStream.precision(PRECISION);
	occurStream<<"POS"<<"\t"<<"occur"<<"\t"<<"unknown"<<endl;

	char missigDataChar = -2;
	char occurChar = 1;
	for(int pos=0; pos<_sc.seqLen(); pos++){
		int NumOfOccurancesPerPos = _sc.getNumOfOccurancesPerPos(pos, occurChar);
		_occurPerPos.push_back(NumOfOccurancesPerPos);
		int NumOfUnknownPerPos = _sc.getNumOfOccurancesPerPos(pos, missigDataChar);
		_unknownPerPos.push_back(NumOfUnknownPerPos);
		occurStream<<pos+1<<"\t"<<NumOfOccurancesPerPos<<"\t"<<NumOfUnknownPerPos<<endl;
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::removePositionsWithHighPercentOfMissingData(MDOUBLE fractionOfMissingDataToRemove){
	char missigDataChar = -2;
	MDOUBLE numberOfSeq =  _sc.numberOfSeqs();
	//MDOUBLE numberOfMissinPerPos;

	vector<int> posToRemove(_sc.seqLen(),false);
	for(int pos=0; pos<_sc.seqLen(); ++pos){
		int NumOfOccurancesPerPos =_unknownPerPos[pos-1]; // pre-computed 
		if( (float)NumOfOccurancesPerPos/numberOfSeq >=  fractionOfMissingDataToRemove ){
			posToRemove[pos] = true;
		}
	}
	_scFilterMissingData = _sc;
	_scFilterMissingData.removePositions(posToRemove);
	LOGnOUT(4,<<"The number of positions with missing less than "<<fractionOfMissingDataToRemove<<" is "<<_scFilterMissingData.seqLen()<<endl);
	string strSeq = gainLossOptions::_outDir + "//" + "seqFilterMissingData."+ double2string(fractionOfMissingDataToRemove) + ".fa";
	ofstream seq_out(strSeq.c_str());
	fastaFormat::write(seq_out,_scFilterMissingData);

	//int startPosInLoop = 1;
	// start the first position
	//for (int firstPos = 0; firstPos<_sc.seqLen(); ++firstPos){
	//	if((float)_sc.getNumOfOccurancesPerPos(firstPos,missigDataChar)/numberOfSeq <  fractionOfMissingDataToRemove ){
	//		_scFilterMissingData = _sc.getSubSeq(firstPos,firstPos); //start with first position
	//		LOGnOUT(5,<<"The first position with more missing data less than "<<fractionOfMissingDataToRemove<<" is "<<firstPos<<endl);
	//		startPosInLoop = firstPos+1;
	//		break;
	//	}
	//}

	//if(_scFilterMissingData.seqLen()<1){
	//	LOGnOUT(4,<<"WARN: there is no position with more than "<<fractionOfMissingDataToRemove<<endl);
	//	return;
	//}
 //   
	//for(int pos=startPosInLoop; pos<_sc.seqLen(); ++pos){
	//	if(pos%1000==0)
	//		cout<<pos<<endl; // DEB
	//	if(! ((float)_sc.getNumOfOccurancesPerPos(pos,missigDataChar)/numberOfSeq >=  fractionOfMissingDataToRemove) ){
	//		_scFilterMissingData.concatenate(_sc.getSubSeq(pos,pos));
	//	}

	//}
	//cout<<_scFilterMissingData.seqLen()<<" "<<_sc.seqLen()<<endl; //DEB
}


/********************************************************************************************
startSequenceContainer
*********************************************************************************************/
void gainLoss::startSequenceContainerUniqPatterns(){
	LOGnOUT(4,<<" *** Starting compute Unique patterns"<<endl);
	
	time_t t1;
	time(&t1);
	time_t t2;

	gainLossAlphabet alph;
	Vint scUniqPatternsNumberOfOnesPerPos;
	vector<sequenceContainer> sequenceContainerVector;
	_scUniqPatterns = _sc.getSubSeq(0,0); //start with first position
	sequenceContainerVector.push_back(_scUniqPatterns.getSubSeq(0,0));
	scUniqPatternsNumberOfOnesPerPos.push_back(_scUniqPatterns.getNumOfOccurancesPerPos(0,1));

	Vint posWeights;
	posWeights.push_back(1);

	for(int pos=1; pos<_sc.seqLen(); ++pos){
		if(pos%1000==0)
			cout<<pos<<endl; // DEB	
		bool isPosUniq = true;
		sequenceContainer seqPos(_sc.getSubSeq(pos,pos));
		int numberOfOnesPerPos = seqPos.getNumOfOccurancesPerPos(0,1);

		for(int i=0; i<_scUniqPatterns.seqLen(); ++i){
			if(scUniqPatternsNumberOfOnesPerPos.size()<i){
				scUniqPatternsNumberOfOnesPerPos.push_back(_scUniqPatterns.getNumOfOccurancesPerPos(i,1));
				sequenceContainerVector.push_back(_scUniqPatterns.getSubSeq(i,i));
			}

			if(numberOfOnesPerPos == scUniqPatternsNumberOfOnesPerPos[i] && sequenceContainerVector[i] == seqPos){
				isPosUniq = false;
				posWeights[i] = posWeights[i]+1;
				break;
			}
		}
		if(isPosUniq){
			_scUniqPatterns.concatenate(seqPos);
			posWeights.push_back(1);
		}
	}
	_weightsUniqPatterns = new Vdouble;
	_weightsUniqPatterns->resize(posWeights.size());

	string posWeightsSt = gainLossOptions::_outDir + "//" + "posWeights" + ".txt";
	ofstream posWeights_out(posWeightsSt.c_str());

	if(posWeights.size() == _scUniqPatterns.seqLen()){
		for(int i=0; i<posWeights.size(); ++i){
			posWeights_out<<posWeights[i]<<endl;	
			(*_weightsUniqPatterns)[i] = posWeights[i];
		}
	}
	else
		errorMsg::reportError("posWeights and _scUniqPatterns - Not the same length");

	string strSeq = gainLossOptions::_outDir + "//" + "seqUniq" + ".fa";
	ofstream seq_out(strSeq.c_str());
	fastaFormat::write(seq_out,_scUniqPatterns);

	//_sc =  sequenceContainer(_scUniqPatterns,&alph);
	time(&t2);
	LOGnOUT(4,<<"Computed Unique pattern RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	LOGnOUT(4,<<"seqLenUnique pattern="<<_scUniqPatterns.seqLen()<<endl);
}




/********************************************************************************************
fillReferenceSequence
*********************************************************************************************/
void gainLoss::fillReferenceSequence(){
	if (strcmp(gainLossOptions::_referenceSeq.c_str(),"non")==0) {
		_refSeq = &(_sc[0]);
	}
	else {
		int id1 = _sc.getId(gainLossOptions::_referenceSeq,true);
		_refSeq = (&_sc[id1]);
	}
}
/********************************************************************************************
startStochasticProcess
*********************************************************************************************/
void gainLoss::startStochasticProcess(bool gainLossDist){
	if(!gainLossDist)
		startStochasticProcess();
	else
		startStochasticProcessVec();	//gain, loss ~ gamma
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::setRootFreq(){
	_freq.resize(gainLossOptions::_alphabet_size);
	if(gainLossOptions::_isRootFreqEQstationary){
		_freq[1]=gainLossOptions::_userGain/(gainLossOptions::_userGain+gainLossOptions::_userLoss);
		_freq[0]=1-_freq[1];
	}else{
		if(gainLossOptions::_userTheta != 0.5){		
			_freq[1]=gainLossOptions::_userTheta;
			_freq[0]=1-_freq[1];
		}
		else if (gainLossOptions::_userGainLossRatio <VERYBIG ){ // then it was given specifically
			_freq[1]= gainLossOptions::_userGainLossRatio/(1+gainLossOptions::_userGainLossRatio);
			_freq[0]=1-_freq[1];
		}
		else{
			_freq = computeFreq();	// if user didnt provide the data computeFreq.
		}
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::startStochasticProcess()
{ 
	LOGnOUT(4,<<"\n startStochasticProcess..."<<endl);	
	MDOUBLE init_gain = gainLossOptions::_userGain;
	MDOUBLE init_loss = gainLossOptions::_userLoss;
	if(Parameters::getInt("_isInitGainLossByEmpiricalFreq")){		
		_freq = evaluateCharacterFreq(_scWithFullLength);
		init_gain = _freq[1];
		init_loss = _freq[0];
		LOGnOUT(4,<<endl<<"InitGainLossByEmpiricalFreq: freq 1=init_gain= "<<_freq[1]<<endl);
	}	
	setRootFreq();
	
	_gainExp = init_gain;
	_lossExp = init_loss;
	replacementModel* glm;
	if(!gainLossOptions::_isReversible){
		glm = new gainLossModelNonReversible(init_gain,init_loss,_freq,gainLossOptions::_isRootFreqEQstationary,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}
	else{
		glm = new gainLossModel(init_gain, _freq,gainLossOptions::_isRootFreqEQstationary, true,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}	
	trivialAccelerator* pijAcc = new trivialAccelerator(glm);	
	MDOUBLE initAlphaRate = gainLossOptions::_userAlphaRate;
	MDOUBLE initBetaRate = gainLossOptions::_userBetaRate;
	MDOUBLE initProbInvariant = gainLossOptions::_userProbInvariantRate;
	MDOUBLE initGlobalRate =1;
	MDOUBLE initRateInvariantVal = gainLossOptions::_userRateInvariantVal;
	
//Mixture
	int numOfGammaComp = gainLossOptions::_numberOfRateComponents;
	int numOfRateCategories = gainLossOptions::_numberOfRateCategories;
	Vdouble initAlphaRates;
	Vdouble initBetaRates;
	Vdouble initCompProbRates;

	distribution* rateDist;
	distribution* baseDistr;
	switch (gainLossOptions::_rateDistributionType){
		case (gainLossOptions::UNIFORM):  
			rateDist = new uniDistribution();
			LOGnOUT(4,<<"rateDist UNIFORM" <<endl);
			break;
		case (gainLossOptions::GAMMA):
			rateDist = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories); //
			LOGnOUT(4,<<"rateDist GAMMA with: initAlphaRate="<<initAlphaRate<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfRateCategories<<endl);
			break;
		case (gainLossOptions::GENERAL_GAMMA):
			rateDist = new generalGammaDistribution(initAlphaRate,initBetaRate,gainLossOptions::_numberOfRateCategories); //
			LOGnOUT(4,<<"rateDist GENERAL_GAMMA with: initAlphaRate="<<initAlphaRate<<" initBetaRate="<<initBetaRate<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfRateCategories<<endl);
			break;
		case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
			if(!( (numOfRateCategories==1)||(numOfRateCategories==2)||(numOfRateCategories==4)||(numOfRateCategories==5)||(numOfRateCategories==16) )){
				string err = "in rate distr. GAMMA_FIXED_CATEGORIES only #cat={1,2,4,5,16} is supported. Not #cat=";
				err+=int2string(numOfRateCategories);
				errorMsg::reportError(err);
			}				
			rateDist = new gammaDistributionFixedCategories(initAlphaRate, gainLossOptions::_numberOfRateCategories);
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfRateCategories<<" initAlphaRate= "<<initAlphaRate<<endl);
			break;
		case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
			if(!( (numOfRateCategories==1)||(numOfRateCategories==2)||(numOfRateCategories==4)||(numOfRateCategories==5)||(numOfRateCategories==16) )){
				string err = "in rate distr. GENERAL_GAMMA_FIXED_CATEGORIES only #cat={1,2,4,5,16} is supported. Not #cat=";
				err+=int2string(numOfRateCategories);
				errorMsg::reportError(err);
			}				
			rateDist = new generalGammaDistributionFixedCategories(initAlphaRate, initBetaRate,gainLossOptions::_numberOfRateCategories);
			LOGnOUT(4,<<"The rateDist was initialized as GENERAL_GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfRateCategories<<" initAlphaRate= "<<initAlphaRate<<" initBetaRate="<<initBetaRate<<endl);
			break;
		case (gainLossOptions::GAMMA_MIXTURE):
			if(gainLossOptions::_initRandomGammaMixuteParam){
				if(gainLossOptions::_rateDiscretizationType == gainLossOptions::QUANTILE)
					rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,QUANTILE);
				if(gainLossOptions::_rateDiscretizationType == gainLossOptions::LAGUERRE)
					rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,LAGUERRE);
			}
			else{
				initMixtureParams(initAlphaRates,initBetaRates,initCompProbRates,numOfGammaComp,initAlphaRate, initBetaRate);
				rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,initAlphaRates,initBetaRates,initCompProbRates);
			}
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_MIXTURE"<<endl);
			break;
		case (gainLossOptions::GAMMA_PLUS_INV):
			baseDistr = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories);
			rateDist = new gammaDistributionPlusInvariant(baseDistr,initProbInvariant,initGlobalRate,initRateInvariantVal);
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_PLUS_INV with: initAlphaRate="<<initAlphaRate<<" and initProbInvariant="<<initProbInvariant<<endl);
			if(baseDistr) delete baseDistr;
			break;
		case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):{
			baseDistr = new generalGammaDistribution(initAlphaRate,initBetaRate,gainLossOptions::_numberOfRateCategories);
			rateDist = new generalGammaDistributionPlusInvariant(baseDistr,initProbInvariant,initGlobalRate,initRateInvariantVal);
			LOGnOUT(4,<<"The rateDist was initialized as GENERAL_GAMMA_PLUS_INV with: initBetaRate="<<initBetaRate<<" and initProbInvariant="<<initProbInvariant<<endl);
			if(baseDistr) delete baseDistr;
		}		
	}
	_sp = new stochasticProcess(rateDist,pijAcc,gainLossOptions::_isReversible);
	MDOUBLE norm_factor = normalizeQ(_sp);
	LOGnOUT(4,<<"Stochastic process normalized with norm_factor="<<norm_factor<<endl);

	if (rateDist) delete rateDist;	//at r4s after the sp object is created all other objects dynamically constructed are deleted
	if (pijAcc) delete pijAcc;
	if (glm) delete glm;
}



/********************************************************************************************
startStochasticProcessGeneric
*********************************************************************************************/
stochasticProcess*  gainLoss::startStochasticProcessGeneric(gainLossOptions::distributionType rateDistributionType, const bool isReversible)
{ 
	LOGnOUT(4,<<"\n startStochasticProcessGeneric..."<<endl);
	
	MDOUBLE init_gain = gainLossOptions::_userGain;
	MDOUBLE init_loss = gainLossOptions::_userLoss;
	if(Parameters::getInt("_isInitGainLossByEmpiricalFreq")){		
		_freq = evaluateCharacterFreq(_scWithFullLength);
		init_gain = _freq[1];
		init_loss = _freq[0];
		LOGnOUT(4,<<endl<<"InitGainLossByEmpiricalFreq: freq 1=init_gain= "<<_freq[1]<<endl);
	}	
	setRootFreq();

	replacementModel* glm;
	if(!isReversible){
		glm = new gainLossModelNonReversible(init_gain,init_loss,_freq,gainLossOptions::_isRootFreqEQstationary,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}
	else{
		glm = new gainLossModel(init_gain, _freq,gainLossOptions::_isRootFreqEQstationary, true,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}	
	trivialAccelerator* pijAcc = new trivialAccelerator(glm);

	MDOUBLE initAlphaRate = gainLossOptions::_userAlphaRate;
	MDOUBLE initBetaRate = gainLossOptions::_userBetaRate;
	MDOUBLE initProbInvariant = gainLossOptions::_userProbInvariantRate;
	MDOUBLE initGlobalRate =1;
	MDOUBLE initRateInvariantVal = gainLossOptions::_userRateInvariantVal;
	//Mixture
	int numOfGammaComp = gainLossOptions::_numberOfRateComponents;
	int numOfRateCategories = gainLossOptions::_numberOfRateCategories;
	Vdouble initAlphaRates;
	Vdouble initBetaRates;
	Vdouble initCompProbRates;

	distribution* rateDist =NULL;
	distribution* baseDistr=NULL;
	switch (rateDistributionType){
		case (gainLossOptions::UNIFORM):  
			rateDist = new uniDistribution();
			LOGnOUT(4,<<"rateDist UNIFORM" <<endl);
			break;
		case (gainLossOptions::GAMMA):
			rateDist = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories); //
			LOGnOUT(4,<<"rateDist GAMMA with: initAlphaRate="<<initAlphaRate<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfRateCategories<<endl);
			break;
		case (gainLossOptions::GENERAL_GAMMA):
			rateDist = new generalGammaDistribution(initAlphaRate,initBetaRate,gainLossOptions::_numberOfRateCategories); //
			LOGnOUT(4,<<"rateDist GENERAL_GAMMA with: initAlphaRate="<<initAlphaRate<<" initBetaRate="<<initBetaRate<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfRateCategories<<endl);
			break;
		case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
			if(!( (numOfRateCategories==1)||(numOfRateCategories==2)||(numOfRateCategories==3)||(numOfRateCategories==5)||(numOfRateCategories==8)||(numOfRateCategories==12||(numOfRateCategories==16)||(numOfRateCategories==24)||(numOfRateCategories==32)||(numOfRateCategories==36)) )){
				string err = "in rate distr. GAMMA_FIXED_CATEGORIES only #cat={1,2,3,5,8,12,16,24,32,36} is supported. Not #cat=";
				err+=int2string(numOfRateCategories);
				errorMsg::reportError(err);
			}				
			rateDist = new gammaDistributionFixedCategories(initAlphaRate, gainLossOptions::_numberOfRateCategories);
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfRateCategories<<" initAlphaRate= "<<initAlphaRate<<endl);
			break;
		case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
			if(!( (numOfRateCategories==1)||(numOfRateCategories==2)||(numOfRateCategories==3)||(numOfRateCategories==5)||(numOfRateCategories==8)||(numOfRateCategories==12)||(numOfRateCategories==16)||(numOfRateCategories==24)||(numOfRateCategories==32)||(numOfRateCategories==36) )){
				string err = "in rate distr. GENERAL_GAMMA_FIXED_CATEGORIES only #cat={1,2,3,5,8,12,16,24,32,36} is supported. Not #cat=";
				err+=int2string(numOfRateCategories);
				errorMsg::reportError(err);
			}				
			rateDist = new generalGammaDistributionFixedCategories(initAlphaRate, initBetaRate,gainLossOptions::_numberOfRateCategories);
			LOGnOUT(4,<<"The rateDist was initialized as GENERAL_GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfRateCategories<<" initAlphaRate= "<<initAlphaRate<<" initBetaRate="<<initBetaRate<<endl);
			break;
		case (gainLossOptions::GAMMA_MIXTURE):
			if(gainLossOptions::_initRandomGammaMixuteParam){
				if(gainLossOptions::_rateDiscretizationType == gainLossOptions::QUANTILE)
					rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,QUANTILE);
				if(gainLossOptions::_rateDiscretizationType == gainLossOptions::LAGUERRE)
					rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,LAGUERRE);
			}
			else{
				initMixtureParams(initAlphaRates,initBetaRates,initCompProbRates,numOfGammaComp,initAlphaRate, initBetaRate);	// standard points
				rateDist = new mixtureDistribution(gainLossOptions::_numberOfRateComponents, gainLossOptions::_numberOfRateCategories,initAlphaRates,initBetaRates,initCompProbRates);
			}
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_MIXTURE"<<endl);
			break;
		case (gainLossOptions::GAMMA_PLUS_INV):
			baseDistr = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories);
			rateDist = new gammaDistributionPlusInvariant(baseDistr,initProbInvariant,initGlobalRate,initRateInvariantVal);
			LOGnOUT(4,<<"The rateDist was initialized as GAMMA_PLUS_INV with: initAlphaRate="<<initAlphaRate<<" and initProbInvariant="<<initProbInvariant<<endl);
			if(baseDistr) delete baseDistr;
			break;
		case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):{
			baseDistr = new generalGammaDistribution(initAlphaRate,initBetaRate,gainLossOptions::_numberOfRateCategories);
			rateDist = new generalGammaDistributionPlusInvariant(baseDistr,initProbInvariant,initGlobalRate,initRateInvariantVal);
			LOGnOUT(4,<<"The rateDist was initialized as GENERAL_GAMMA_PLUS_INV with: initBetaRate="<<initBetaRate<<" and initProbInvariant="<<initProbInvariant<<endl);
			if(baseDistr) delete baseDistr;
													   }		
	}	 
	stochasticProcess* sp = new stochasticProcess(rateDist,pijAcc,gainLossOptions::_isReversible);

	MDOUBLE norm_factor = normalizeQ(sp);
	LOGnOUT(4,<<"Stochastic process normalized with norm_factor="<<norm_factor<<endl);
	if (rateDist) delete rateDist;	//at r4s after the sp object is created all other objects dynamically constructed are deleted
	if (pijAcc) delete pijAcc;
	if (glm) delete glm;
	return sp;
}
/********************************************************************************************
startStochasticProcessVec
*********************************************************************************************/
void gainLoss::startStochasticProcessVec(){ 
	LOGnOUT(4,<<"\n startStochasticProcessVec with: GainCategories="<<gainLossOptions::_numberOfGainCategories<<" and GainCategories="<<gainLossOptions::_numberOfLossCategories<<endl);
	bool isReversible =gainLossOptions::_isReversible;
	//Vdouble freq;
	//Vdouble freqEmpirical;
	MDOUBLE init_gain;
	MDOUBLE init_loss;
	MDOUBLE initAlphaRate = gainLossOptions::_userAlphaRate;	
	MDOUBLE initAlphaGain = gainLossOptions::_userAlphaGain;
	MDOUBLE initBetaGain = gainLossOptions::_userBetaGain;		
	MDOUBLE initAlphaLoss = gainLossOptions::_userAlphaLoss;
	MDOUBLE initBetaLoss = gainLossOptions::_userBetaLoss;
	MDOUBLE initProbInvariantGain = gainLossOptions::_userProbInvariantGain;
	MDOUBLE initProbInvariantLoss = gainLossOptions::_userProbInvariantLoss;
	MDOUBLE initGlobalRate =1;
	MDOUBLE initRateInvariantVal = gainLossOptions::_userRateInvariantVal;

	setRootFreq();

	if(Parameters::getInt("_isInitGainLossByEmpiricalFreq")){
		_freq = evaluateCharacterFreq(_scWithFullLength);
		init_gain = _freq[1];
		init_loss = _freq[0];
		MDOUBLE gainLossRatioToCompleteByBeta = (init_gain/init_loss)*(initAlphaLoss/initAlphaGain);
		LOGnOUT(4,<<endl<<"InitGainLossByEmpiricalFreq: freq 1= "<<_freq[1]<<" Thus, gainLossRatioToCompleteByBeta= "<<gainLossRatioToCompleteByBeta<<endl);
		if(gainLossOptions::_isUpdateOnlyGainBetaForRatio)
			initBetaGain =(initBetaLoss/gainLossRatioToCompleteByBeta);			// AlphaGain = 0.35
		else{
			initBetaGain =sqrt(1/gainLossRatioToCompleteByBeta);			// AlphaGain = 0.35
			initBetaLoss =sqrt(gainLossRatioToCompleteByBeta);				// AlphaLoss = 0.9
		}
	}

// gain		
	switch (gainLossOptions::_gainDistributionType){
			case (gainLossOptions::UNIFORM):  
				_gainDist = new uniDistribution();
				LOGnOUT(4,<<"rateDist UNIFORM" <<endl);
				break;
			case (gainLossOptions::GAMMA):
				_gainDist = new gammaDistribution(initAlphaGain,gainLossOptions::_numberOfGainCategories); //
				LOGnOUT(4,<<"gainDist GAMMA with: initAlpha="<<initAlphaGain<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfGainCategories<<endl);
				break;
			case (gainLossOptions::GENERAL_GAMMA):
				_gainDist = new generalGammaDistribution(initAlphaGain,initBetaGain,gainLossOptions::_numberOfGainCategories); //
				LOGnOUT(4,<<"gainDist GENERAL_GAMMA with: initAlphaGain="<<initAlphaGain<<" initBetaGain="<<initBetaGain<<" and _numberOfGainCategories= "<<gainLossOptions::_numberOfGainCategories<<endl);
				break;
			case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
				//if(!( (numOfRateCategories==1)||(numOfRateCategories==2)||(numOfRateCategories==3)||(numOfRateCategories==5)||(numOfRateCategories==8)||(numOfRateCategories==12||(numOfRateCategories==16)||(numOfRateCategories==24)||(numOfRateCategories==32)||(numOfRateCategories==36)) )){
				//	string err = "in gain distr. GAMMA_FIXED_CATEGORIES only #cat={1,2,3,5,8,12,16,24,32,36} is supported. Not #cat=";
				//	err+=int2string(numOfRateCategories);
				//	errorMsg::reportError(err);
				//}				
				_gainDist = new gammaDistributionFixedCategories(initAlphaGain, gainLossOptions::_numberOfGainCategories);
				LOGnOUT(4,<<"The _gainDist was initialized as GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfGainCategories<<" initAlphaGain= "<<initAlphaGain<<endl);
				break;
			case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
				errorMsg::reportError("in gainDist. GENERAL_GAMMA_FIXED_CATEGORIES is not realized");
				break;
			case (gainLossOptions::GAMMA_MIXTURE):
				errorMsg::reportError("in gainDist. GAMMA_MIXTURE is not realized");
				break;
			case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):{
				distribution* baseDistr = new generalGammaDistribution(initAlphaGain,initBetaGain,gainLossOptions::_numberOfGainCategories);
				_gainDist = new generalGammaDistributionPlusInvariant(baseDistr,initProbInvariantGain,initGlobalRate,initRateInvariantVal);
				LOGnOUT(4,<<"The gainDist was initialized as GENERAL_GAMMA_PLUS_INV with: initAlphaGain="<<initAlphaGain<<" initBetaGain="<<initBetaGain<<" initProbInvariantGain="<<initProbInvariantGain<<" without optimization"<<endl);
				delete baseDistr;
				break;}				
			default:
				errorMsg::reportError("error in startStochasticProcessVec the distribution chosen is not implemented");
		}
// loss		
	switch (gainLossOptions::_lossDistributionType){
			case (gainLossOptions::UNIFORM):  
				_lossDist = new uniDistribution();
				LOGnOUT(4,<<"rateDist UNIFORM" <<endl);
				break;
			case (gainLossOptions::GAMMA):
				_lossDist = new gammaDistribution(initAlphaLoss,gainLossOptions::_numberOfLossCategories); //
				LOGnOUT(4,<<"lossDist GAMMA with: initAlpha="<<initAlphaLoss<<" and _numberOfRateCategories= "<<gainLossOptions::_numberOfLossCategories<<endl);
				break;
			case (gainLossOptions::GENERAL_GAMMA):
				_lossDist = new generalGammaDistribution(initAlphaLoss,initBetaLoss,gainLossOptions::_numberOfLossCategories); //
				LOGnOUT(4,<<"lossDist GENERAL_GAMMA with: initAlphaLoss="<<initAlphaLoss<<" initBetaLoss="<<initBetaLoss<<" and _numberOfLossCategories= "<<gainLossOptions::_numberOfLossCategories<<endl);
				break;
			case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
				_lossDist = new gammaDistributionFixedCategories(initAlphaLoss, gainLossOptions::_numberOfLossCategories);
				LOGnOUT(4,<<"The _gainDist was initialized as GAMMA_FIXED_CATEGORIES with: num Of Categories "<<gainLossOptions::_numberOfLossCategories<<" initAlphaLoss= "<<initAlphaLoss<<endl);
				errorMsg::reportError("in lossDist. GAMMA_FIXED_CATEGORIES is not realized");
				break;
			case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
				errorMsg::reportError("in lossDist. GENERAL_GAMMA_FIXED_CATEGORIES is not realized");
				break;
			case (gainLossOptions::GAMMA_MIXTURE):
				errorMsg::reportError("in lossDist. GAMMA_MIXTURE is not realized");
				break;
			case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):{
				distribution* baseDistr = new generalGammaDistribution(initAlphaLoss,initBetaLoss,gainLossOptions::_numberOfLossCategories);
				_lossDist = new generalGammaDistributionPlusInvariant(baseDistr,initProbInvariantLoss,initGlobalRate,initRateInvariantVal);
				LOGnOUT(4,<<"The lossDist was initialized as GENERAL_GAMMA_PLUS_INV with: initAlphaLoss="<<initAlphaLoss<<" initBetaLoss="<<initBetaLoss<<" initProbInvariantLoss="<<initProbInvariantLoss<<" without optimization"<<endl);
				delete baseDistr;  
				break;}
			default:errorMsg::reportError("error in startStochasticProcessVec the distribution chosen is not implemented");
	}

// Loop over gain and loss distributions
	distribution* baseDistr;
	_spVVec.resize(_gainDist->categories());
	for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
		_spVVec[gainCategor].resize(_lossDist->categories());
		for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){			
			replacementModel* glm;
			if(!isReversible){
				glm = new gainLossModelNonReversible(_gainDist->rates(gainCategor),_lossDist->rates(lossCategor),_freq,gainLossOptions::_isRootFreqEQstationary,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
			}
			else{
				glm = new gainLossModel(_gainDist->rates(gainCategor),_freq,gainLossOptions::_isRootFreqEQstationary, true,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
			}			
			pijAccelerator* pijAcc = new trivialAccelerator(glm);

			distribution* rateDist;
			switch (gainLossOptions::_rateDistributionType){
				case (gainLossOptions::UNIFORM):  
					rateDist = new uniDistribution();
					break;
				case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
					rateDist = new gammaDistributionFixedCategories(initAlphaRate, gainLossOptions::_numberOfRateCategories);
					break;
				case (gainLossOptions::GAMMA):
					rateDist = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories); //
					break;
				case (gainLossOptions::GAMMA_PLUS_INV):
					baseDistr = new gammaDistribution(initAlphaRate,gainLossOptions::_numberOfRateCategories);
					rateDist = new gammaDistributionPlusInvariant(baseDistr,gainLossOptions::_userProbInvariantRate,initGlobalRate,initRateInvariantVal);
					if(baseDistr) delete baseDistr;
					break;			
				default:
					errorMsg::reportError("unknown type in distributionType");
			}
			stochasticProcess* sp = new stochasticProcess(rateDist,pijAcc,gainLossOptions::_isReversible);			
			_spVVec[gainCategor][lossCategor] = sp->clone();
			
			if (rateDist) delete rateDist;	//at r4s after the sp object is created all other objects dynamically constructed are deleted
			if (pijAcc) delete pijAcc;
			if (glm) delete glm;
			if (sp) delete sp;
		}
	}
	_gainExp = rateExpectation(_gainDist);
	_lossExp = rateExpectation(_lossDist);

	MDOUBLE norm_factor = normalizeQ(_spVVec, _gainDist, _lossDist);
	LOGnOUT(4,<<"Stochastic process vector normalized with norm_factor="<<norm_factor<<endl);
	_sp = _spVVec[0][0]; // initialize the "normal _sp" data member
}
/********************************************************************************************
computeFreq
*********************************************************************************************/
Vdouble gainLoss::computeFreq(){
	Vdouble freq;
	switch (gainLossOptions::_characterFreqEval){
		case (gainLossOptions::FiftyFifty):  
			freq.push_back(0.5);	// initializing the frequency vector to 0.5
			freq.push_back(0.5);
			LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<freq[0]<<" "<<freq[1]<<endl);
			break;
		case (gainLossOptions::LeavesAve):
			freq = evaluateCharacterFreq(_scWithFullLength); //
			LOGnOUT(4,<<"frequencies are based on LeavesAve "<<freq[0]<<" "<<freq[1]<<endl);
			break;
		case (gainLossOptions::optimizeOverTree):
			freq = evaluateCharacterFreq(_scWithFullLength); // the rest will be be preformed during the optimization stage
			LOGnOUT(4,<<"frequencies are "<<freq[0]<<" "<<freq[1]<<endl);
			break;
	}
	return freq;
}
/********************************************************************************************
startingEvolTreeTopology
*********************************************************************************************/
void gainLoss::startEvolTreeTopology(ostream& out){
	//time_t ltime1;
	//time( &ltime1 );
	LOGnOUT(4,<<"\n startingEvolTreeTopology..."<<endl);
	VVdouble disTab;
	vector<string> vNames;
	if (gainLossOptions::_treeFile=="") {
		LOGnOUT(4,<<"No treeFile was given. The tree will be estimated from distance matrix"<<endl);
		distanceMethod* pDm;
		switch (gainLossOptions::_treeSearchAlg){
			case (gainLossOptions::njJC):
				pDm = new jcDistance();
				giveDistanceTable(pDm, _sc, disTab, vNames);
				break;
			case (gainLossOptions::njJCOLD):
				pDm = new jcDistanceOLD(gainLossOptions::_alphabet_size);
				giveDistanceTable(pDm, _sc,disTab, vNames);
				break;
			case (gainLossOptions::njML): {
				uniDistribution lUni;
				const pijAccelerator* lpijAcc = _spSimple->getPijAccelerator();// note this is just a copy of the pointer.
				stochasticProcess lsp(&lUni,lpijAcc);
				pDm = new likeDist(lsp,0.01);
				//pDm = new likeDist(*_spSimple);	// in this sp the gain and loss are taken from empirical freq and gamma dist is used
				giveDistanceTable(pDm,_sc,disTab,vNames);
				}
				break;
			default:
				errorMsg::reportError("this tree search mode is not yet available");
		}
		delete pDm;

		//calc distance table statistics
		MDOUBLE low_bound = VERYBIG;
		MDOUBLE upper_bound = VERYSMALL;
		MDOUBLE sum = 0.0;
		int count = 0;
		for (int i = 0; i < disTab.size(); ++i){
			for (int j = i+1; j < disTab[i].size(); ++j){
				sum += disTab[i][j];
				++count;
				if (disTab[i][j] < low_bound)
					low_bound = disTab[i][j];
				if (disTab[i][j] > upper_bound)
					upper_bound = disTab[i][j];
			}
		}
		MDOUBLE avg = sum / static_cast<MDOUBLE>(count);
		LOG(5,<<"#MSA diversity matrix"<<endl);
		LOG(5,<<"#Average pairwise distance= "<<avg<<endl);
		LOG(5,<<"#lower bound = "<<low_bound<<endl);
		LOG(5,<<"#upper bound = "<<upper_bound<<endl);
		LOG(5,<<"#end of MSA diversity matrix"<<endl);
		getStartingTreeNJ_fromDistances(disTab, vNames);
	}
	
	else 
		getStartingTreeFromTreeFile();

	if (!(gainLossOptions::_rootAt =="")){
		tree::nodeP myroot = _tr.findNodeByName(gainLossOptions::_rootAt); //returns NULL if not found
		if (myroot){
			_tr.rootAt(myroot);
			LOGnOUT(4,<<"tree rooted at "<<myroot->name()<<"\n sons of root are:"<<endl);
			for(int i = 0; i<_tr.getRoot()->getNumberOfSons(); ++i ){
				LOGnOUT(4,<<_tr.getRoot()->getSon(i)->name()<<"  ");
			}
			LOGnOUT(4,<<"\n");
			return;
		}
	}
	LOGnOUT(4,<<"default rooting used, root name is "<<_tr.getRoot()->name()<<endl);
	LOGnOUT(4,<<"sons of root are:"<<endl);
	for(int i = 0; i<_tr.getRoot()->getNumberOfSons(); ++i ){
		LOGnOUT(4,<<_tr.getRoot()->getSon(i)->name()<<"  ");
	}
	LOGnOUT(4,<<"\n");

	//return;

	if(gainLossOptions::_seqFile!="" && !_tr.getLeavesNum()==_sc.numberOfSeqs()){
		errorMsg::reportError("The number of sequence is not equal to the number of taxas in the tree");
	}
	_tr.makeSureAllBranchesAreLargerThanEpsilon(gainLossOptions::_minBranchLength);
	_trOrig = _tr;

	
	//time_t ltime2;
	//time( &ltime2 );
	//int t = static_cast<long>(ltime2 - ltime1);
	//timingsF<<"time for tree topology = "<<t<<endl;
}
/********************************************************************************************
getStartingTreeFromTreeFile
*********************************************************************************************/
void gainLoss::getStartingTreeFromTreeFile(){
	_tr= tree(gainLossOptions::_treeFile);
	//if (!_tr.withBranchLength()) _tr.createFlatLengthMatrix(0.1);	// not required, checked before
}

/********************************************************************************************
getStartingTreeNJ_fromDistances
*********************************************************************************************/
void gainLoss::getStartingTreeNJ_fromDistances(const VVdouble& disTab,const vector<string>& vNames) {
	NJalg nj1;
	_tr= nj1.computeTree(disTab,vNames);
	ofstream f;
	string fileName1=gainLossOptions::_treeOutFile;
	f.open(fileName1.c_str());
	_tr.output(f);
	f.close();
}
/********************************************************************************************
*********************************************************************************************/
//void gainLoss::initMissingDataInfo(){
//	//if(gainLossOptions::_accountForMissingData){
//	//	//gainLossAlphabet alph;		
//	//	//_scZero.startZeroSequenceContainerGL(_sc,gainLossAlphabet());
//	//	//_LforMissingDataPerCat.resize(_sp->categories());
//	//	//_pLforMissingDataPerCat = &_LforMissingDataPerCat;
//
//
//	//	//_plogLforMissingData = &_logLforMissingData;
//	//	//*_plogLforMissingData = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scZero,*_sp);
//
//	//	computePijGam pi;
//	//	pi.fillPij(_tr,*_sp);
//	//	*_pLforMissingDataPerCat = likelihoodComputationGL::getLofPosPerCat(0,_tr,_scZero,pi,*_sp);
//	//	//*_plogLforMissingData = log(likelihoodComputationGL::getLofPos(0,_tr,_scZero,pi,*_sp,*_pLforMissingDataPerCat)); // cause error in tree destructor
//	//}
//	//else{
//	//	//_plogLforMissingData = NULL;
//	//	_pLforMissingDataPerCat = NULL;	
//	//}
//}



// Optimizations
/********************************************************************************************
*********************************************************************************************/
void gainLoss::startOptimizations(){
	LOGnOUT(4,<<"\n\n *** Start Optimizations"<<endl);
	time_t t1,t2;
	time(&t1);

	//bool isScaleTree = false;
	bool isBBL = true; // in optimizer also check gainLossOptions::_isBBL. This if to differ from manyStarts
	MDOUBLE epsilonOptimizationCorrected = gainLossOptions::_epsilonOptimizationIterationCycle;
	MDOUBLE epsilonOptimizationModelCorrected  = gainLossOptions::_epsilonOptimizationModel;
	MDOUBLE epsilonOptimizationBBLCorrected  = gainLossOptions::_epsilonOptimizationBBL;;
	if(gainLossOptions::_correctOptimizationEpsilon && abs(_logL)>100 ){	// if Likelihood is computed with very bad seed - misleading
		LOGnOUT(4,<<" Modify epsilonOptimizations according to logL: logL ("<<abs(_logL)<<") * originalEpsilon * percent improve ("<<gainLossOptions::_percentOfImprov<<") "<<endl);
		epsilonOptimizationCorrected =		min(abs(_logL) * gainLossOptions::_epsilonOptimizationIterationCycle *  gainLossOptions::_percentOfImprov,  gainLossOptions::_epsilonOptimizationIterationCycle*10);
		epsilonOptimizationModelCorrected = min(abs(_logL) * gainLossOptions::_epsilonOptimizationModel *  gainLossOptions::_percentOfImprov,			gainLossOptions::_epsilonOptimizationModel*10);
		epsilonOptimizationBBLCorrected =	min(abs(_logL) * gainLossOptions::_epsilonOptimizationBBL * gainLossOptions::_percentOfImprov,				gainLossOptions::_epsilonOptimizationBBL*10);
		LOGnOUT(4,<<"eOptCorrected(cycle)=\t"<<epsilonOptimizationCorrected<<"\neOptModelCorrected=\t"<<epsilonOptimizationModelCorrected<<"\neOptBBLCorrected=\t"<<epsilonOptimizationBBLCorrected<<endl);		
	}
	if(_sc.seqLen()<50){
		LOGnOUT(4,<<"\n WARN: no branch length estimation is performed with too few positions ="<<_sc.seqLen()<<endl);
		Parameters::updateParameter("_isBBLEMwithSimpleSpBeforeFullOptimization","0");
		Parameters::updateParameter("_performOptimizationsBBL","0");
		isBBL =false;
	}
	
	// NOTE:! This is within block of _performOptimizations to allow, if explicit param request - i.e., only Tree optimization
	if(Parameters::getInt("_isMultipleAllBranchesByFactorAtStart") ){ // else unstable
		if(_sc.seqLen()<50){
			LOGnOUT(4,<<"\n WARN: Skip MultipleAllBranchesByFactorAtStart with too few number of positions "<<_sc.seqLen()<<endl);
			multipleAllBranchesByFactorAtStartByMaxParsimonyCost(_CostOfTreeMP);

		}else{
			multipleAllBranchesByFactorAtStart(epsilonOptimizationBBLCorrected);
		}		
	}
	
	if(Parameters::getInt("_isBBLEMwithSimpleSpBeforeFullOptimization") ){
		bBLEMwithSimpleSpBeforeFullOptimization(_tr,_sc,_spSimple,_sp,_spVVec,_gainDist,_lossDist,_unObservableData_p);		
	}

	if(Parameters::getInt("_performOptimizations") ){
// correctOptimizationEpsilon


// optimize one Stochastic process
		if(!gainLossOptions::_gainLossDist){
			if(gainLossOptions::_initParamsAtRandPoints) initParamsAtRandPoints(gainLossOptions::_numberOfRandStartPoints,_sp,_unObservableData_p);	
			if(gainLossOptions::_performOptimizationsManyStarts)
				optimizationsManyStarts(gainLossOptions::_epsilonOptimizationIterationCycleManyStarts,gainLossOptions::_maxNumOfIterationsManyStarts);			
			gainLossOptimizer glOpt(_tr,_sp,_scUniqPatterns,
				epsilonOptimizationCorrected,gainLossOptions::_maxNumOfIterations,
				epsilonOptimizationModelCorrected,gainLossOptions::_maxNumOfIterationsModel,
				epsilonOptimizationBBLCorrected,gainLossOptions::_maxNumOfIterationsBBL,_weightsUniqPatterns,_unObservableData_p,
				isBBL, gainLossOptions::_isbblLSWhenbblEMdontImprove);
			_tr = glOpt.getOptTree();
			_logL = glOpt.getBestL();
			LOGnOUT(4,<<"# Best likelihood after optimization="<<_logL<<endl);
		}

// optimize Mixture of Stochastic processes
		else{
			if(gainLossOptions::_initParamsAtRandPoints) initParamsAtRandPointsSPvv(gainLossOptions::_numberOfRandStartPoints,_spVVec,_gainDist,_lossDist,_unObservableData_p);	
			if(gainLossOptions::_performOptimizationsManyStarts)
				optimizationsVVManyStarts(gainLossOptions::_epsilonOptimizationIterationCycleManyStarts,gainLossOptions::_maxNumOfIterationsManyStarts);			
			gainLossOptimizer glOpt(_tr,_spVVec,_gainDist,_lossDist,_scUniqPatterns,
				epsilonOptimizationCorrected,gainLossOptions::_maxNumOfIterations,
				epsilonOptimizationModelCorrected,gainLossOptions::_maxNumOfIterationsModel,
				epsilonOptimizationBBLCorrected,gainLossOptions::_maxNumOfIterationsBBL,_weightsUniqPatterns,_unObservableData_p,
				isBBL, gainLossOptions::_isbblLSWhenbblEMdontImprove); // only one set: epsilon,Iter for Model,BBL...
			_tr = glOpt.getOptTree();
			_logL = glOpt.getBestL();
			LOGnOUT(4,<<"# Best likelihood after optimization="<<_logL<<endl);
		}
	}
// No Optimizations
	else{
		LOGnOUT(4,<<"NOTE: No optimization performed. Proceed with initialized parameter."<<endl);
	}
// common lines:
	bool _isTransferGainLossRateToFreq = false;	// there is no point, if normalizing Q later...
	if(_isTransferGainLossRateToFreq){
		convertGainLossRatesToFreq();
	}	
	if(gainLossOptions::_isAlphaEqBetaManipulation && _lossDist && _gainDist && isBetaOptimization(_lossDist) && isBetaOptimization(_gainDist) && !gainLossOptions::_isOptimizeGainLossRatioInsteadOfGainAndLossSeperately /*&& gainLossOptions::_performOptimizationsBBL*/){
		AlphaEqBetaManipulation();
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::printModellValuesOfParams()
{

	string modelParams = gainLossOptions::_outDir + "//" + "EstimatedParameters.txt"; 
	ofstream modelParamsStream(modelParams.c_str());	
	modelParamsStream<<"# Log-likelihood= "<<_logL<<endl;

	if(gainLossOptions::_gainLossDist){
		MDOUBLE bestGainAlpha=1;
		MDOUBLE bestGainBeta=1;
		if(isAlphaOptimization(_gainDist)){
			bestGainAlpha=getRateAlpha(_gainDist);
			LOGnOUT(4,<<"AlphaGain "<<bestGainAlpha <<endl);
			modelParamsStream<<"_userAlphaGain "<<bestGainAlpha<<endl;
		}
		//if(isBetaOptimization(_gainDist))LOGnOUT(4,<<" BetaGain "<<getRateBeta(_gainDist) <<endl);
		if(isBetaOptimization(_gainDist)){
			bestGainBeta=getRateBeta(_gainDist);
			LOGnOUT(4,<<"BetaGain "<<bestGainBeta<<endl);
			LOGnOUT(4,<<"  Gain Expectation = "<< rateExpectation(_gainDist)<<endl);
			LOGnOUT(4,<<"  Gain Expectancy = "<< bestGainAlpha/bestGainBeta<<endl);
			LOGnOUT(4,<<"  Gain Standard Deviation = "<< sqrt(bestGainAlpha/(bestGainBeta*bestGainBeta))<<endl);
			modelParamsStream<<"_userBetaGain "<<bestGainBeta<<endl;
			modelParamsStream<<"#  Gain Expectation = "<<rateExpectation(_gainDist)<<endl;
			modelParamsStream<<"#  Gain Expectancy = "<<bestGainAlpha/bestGainBeta<<endl;
			modelParamsStream<<"#  Gain Standard Deviation = "<<sqrt(bestGainAlpha/(bestGainBeta*bestGainBeta))<<endl;
		}
		if(isInvariantOptimization(_gainDist, true)){
			MDOUBLE probInvariantGain = static_cast<generalGammaDistributionPlusInvariant*>(_gainDist)->getInvProb();
			LOGnOUT(4,<<" ProbInvariantGain "<<probInvariantGain <<endl);
			modelParamsStream<<"_userProbInvariantGain "<<probInvariantGain<<endl;
		}
		MDOUBLE bestLossAlpha=1;
		MDOUBLE bestLossBeta=1;
		if(isAlphaOptimization(_lossDist)){
			bestLossAlpha = getRateAlpha(_lossDist);
			LOGnOUT(4,<<"AlphaLoss "<<bestLossAlpha <<endl);
			modelParamsStream<<"_userAlphaLoss "<<bestLossAlpha<<endl;
		}			
		if(isBetaOptimization(_lossDist)){
			bestLossBeta=getRateBeta(_lossDist);
			LOGnOUT(4,<<"BetaLoss "<<bestLossBeta<<endl);
			LOGnOUT(4,<<"  Loss Expectation = "<< rateExpectation(_lossDist)<<endl);
			LOGnOUT(4,<<"  Loss Expectancy = "<< bestLossAlpha/bestLossBeta<<endl);
			LOGnOUT(4,<<"  Loss Standard Deviation = "<< sqrt(bestLossAlpha/(bestLossBeta*bestLossBeta))<<endl);
			modelParamsStream<<"_userBetaLoss "<<bestLossBeta<<endl;
			modelParamsStream<<"#  Loss Expectation = "<<rateExpectation(_lossDist)<<endl;
			modelParamsStream<<"#  Loss Expectancy = "<<bestLossAlpha/bestLossBeta<<endl;
			modelParamsStream<<"#  Loss Standard Deviation = "<<sqrt(bestLossAlpha/(bestLossBeta*bestLossBeta))<<endl;
		}			
		if(isInvariantOptimization(_lossDist, true)){
			MDOUBLE probInvariantLoss = static_cast<generalGammaDistributionPlusInvariant*>(_lossDist)->getInvProb();
			LOGnOUT(4,<<" ProbInvariantLoss "<<probInvariantLoss <<endl);
			modelParamsStream<<"_userProbInvariantLoss "<<probInvariantLoss<<endl;
		}
		LOGnOUT(4,<<"	Expectancy(Gain)/Expectancy(Loss)  ratio by Gamma Params= "<< (bestGainAlpha/bestGainBeta)/(bestLossAlpha/bestLossBeta)<<endl);
		//LOGnOUT(4,<<"	Expectancy(Gain/Loss) ratio by computation = "<< computeExpectationOfGainLossRatio(_gainDist, _lossDist)<<endl);
		LOGnOUT(4,<<"	Expectancy(Gain)/Expectancy(Loss)  by computation = "<< computeExpOfGainByExpOfLossRatio(_gainDist, _lossDist)<<endl);
		modelParamsStream<<"# GainLossRatio Expectation "<<computeExpOfGainByExpOfLossRatio(_gainDist, _lossDist)<<endl;
		modelParamsStream<<"_userGainLossRatio "<<(bestGainAlpha/bestGainBeta)/(bestLossAlpha/bestLossBeta)<<endl;

		if (gainLossOptions::_isRootFreqEQstationary){
			MDOUBLE estimatedStationaryFreq = computeExpectationOfStationaryFrequency(_gainDist, _lossDist);
			LOGnOUT(4,<<" Stationary '1' Freq at the root - for each stochastic process g/(g+l), with expectation of "<<estimatedStationaryFreq<<endl);
			modelParamsStream<<"# Stationary '1' Freq at the root - for each stochastic process g/(g+l), with expectation of "<<estimatedStationaryFreq<<endl;
		}
		else{
			MDOUBLE thetaVal = static_cast<gainLossModel*>((*_spVVec[0][0]).getPijAccelerator()->getReplacementModel())->getTheta();
			switch (gainLossOptions::_characterFreqEval){
				case (gainLossOptions::FiftyFifty):  
					LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<thetaVal<<endl);
					modelParamsStream<<"# frequencies were set to FiftyFifty "<<thetaVal<<endl;
					break;
				case (gainLossOptions::LeavesAve):
					LOGnOUT(4,<<"frequencies are based on LeavesAve (-F option) "<<thetaVal<<endl);
					modelParamsStream<<"# frequencies are based on LeavesAve (-F option) "<<thetaVal<<endl;
					break;
				case (gainLossOptions::optimizeOverTree):
					LOGnOUT(4,<<"Theta "<<thetaVal <<endl);
					modelParamsStream<<"_userTheta "<<thetaVal<<endl;
					break;
			}
		}
	}

	else{
		if(isAlphaOptimization(_sp->distr())){
			LOGnOUT(4,<<" AlphaRate "<<getRateAlpha(_sp->distr()) <<endl);
			modelParamsStream<<"_userAlphaRate "<<getRateAlpha(_sp->distr())<<endl;
		}
		if(isBetaOptimization(_sp->distr())){
			LOGnOUT(4,<<" BetaRate "<<getRateBeta(_sp->distr()) <<endl);
			modelParamsStream<<"_userBetaRate "<<getRateBeta(_sp->distr())<<endl;
		}
		if(isInvariantOptimization(_sp->distr(), true)){
			MDOUBLE probInvariantRate = static_cast<generalGammaDistributionPlusInvariant*>(_sp->distr())->getInvProb();
			LOGnOUT(4,<<" ProbInvariantRate "<<probInvariantRate <<endl);
			modelParamsStream<<"_userProbInvariantRate "<<probInvariantRate<<endl;
		}

		MDOUBLE gain = static_cast<gainLossModel*>((*_sp).getPijAccelerator()->getReplacementModel())->getMu1();
		LOGnOUT(4,<<" Gain "<<gain <<endl);
		modelParamsStream<<"_userGain "<<gain<<endl;
		MDOUBLE loss;
		if(!gainLossOptions::_isReversible){
			loss = static_cast<gainLossModelNonReversible*>((*_sp).getPijAccelerator()->getReplacementModel())->getMu2();				
			LOGnOUT(4,<<" Loss "<<loss<<endl);
			modelParamsStream<<"_userLoss "<<loss<<endl;
		}
		else{
			loss = static_cast<gainLossModel*>((*_sp).getPijAccelerator()->getReplacementModel())->getMu2();
		}
		LOGnOUT(4,<<"	Gain/Loss  ratio= "<< gain/loss<<endl);		
		modelParamsStream<<"_userGainLossRatio "<<gain/loss<<endl;

		if (gainLossOptions::_isRootFreqEQstationary){
			LOGnOUT(4,<<" Stationary '1' Freq at the root (g/(g+l) = "<<gain/(gain+loss) <<endl);
			modelParamsStream<<"# Stationary '1' Freq at the root (g/(g+l) = "<< gain/(gain+loss)<<endl;
		}
		else{
			MDOUBLE thetaVal = static_cast<gainLossModel*>((*_sp).getPijAccelerator()->getReplacementModel())->getTheta();
			switch (gainLossOptions::_characterFreqEval){
					case (gainLossOptions::FiftyFifty):  
						LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<thetaVal<<endl);
						modelParamsStream<<"# frequencies were set to FiftyFifty "<<thetaVal<<endl;
						break;
					case (gainLossOptions::LeavesAve):
						LOGnOUT(4,<<"frequencies are based on LeavesAve (-F option) "<<thetaVal<<endl);
						modelParamsStream<<"# frequencies are based on LeavesAve (-F option) "<<thetaVal<<endl;
						break;
					case (gainLossOptions::optimizeOverTree):
						LOGnOUT(4,<<" Theta ('1' at root):"<<thetaVal <<endl);
						modelParamsStream<<"_userTheta "<<thetaVal<<endl;
						break;
			}
		}
	}
	LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
	modelParamsStream<<"# Total branch lengths:"<<_tr.getAllBranchesLengthSum()<<endl;
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::printModellValuesOfParams(tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist)
{
	MDOUBLE bestGainAlpha=1;
	MDOUBLE bestGainBeta=1;
	if(isAlphaOptimization(gainDist)){
		bestGainAlpha=getRateAlpha(gainDist);
		LOGnOUT(4,<<"AlphaGain "<<bestGainAlpha <<endl);
	}
	if(isBetaOptimization(gainDist)){
		bestGainBeta=getRateBeta(gainDist);
		LOGnOUT(4,<<"BetaGain "<<bestGainBeta<<endl);
		LOGnOUT(4,<<"  Gain Expectation = "<< rateExpectation(_gainDist)<<endl);
		LOGnOUT(4,<<"  Gain Expectancy = "<< bestGainAlpha/bestGainBeta<<endl);
		LOGnOUT(4,<<"  Gain Standard Deviation = "<< sqrt(bestGainAlpha/(bestGainBeta*bestGainBeta))<<endl);
	}
	if(isInvariantOptimization(gainDist, true)) LOGnOUT(4,<<" ProbInvariantGain "<<static_cast<generalGammaDistributionPlusInvariant*>(gainDist)->getInvProb() <<endl);
	MDOUBLE bestLossAlpha=1;
	MDOUBLE bestLossBeta=1;
	if(isAlphaOptimization(lossDist)){
		bestLossAlpha = getRateAlpha(lossDist);
		LOGnOUT(4,<<"AlphaLoss "<<bestLossAlpha <<endl);
	}			
	if(isBetaOptimization(lossDist)){
		bestLossBeta=getRateBeta(lossDist);
		LOGnOUT(4,<<"BetaLoss "<<bestLossBeta<<endl);
		LOGnOUT(4,<<"  Loss Expectation = "<< rateExpectation(_lossDist)<<endl);
		LOGnOUT(4,<<"  Loss Expectancy = "<< bestLossAlpha/bestLossBeta<<endl);
		LOGnOUT(4,<<"  Loss Standard Deviation = "<< sqrt(bestLossAlpha/(bestLossBeta*bestLossBeta))<<endl);
	}			
	if(isInvariantOptimization(lossDist, true))LOGnOUT(4,<<" ProbInvariantLoss "<<static_cast<generalGammaDistributionPlusInvariant*>(lossDist)->getInvProb() <<endl);

	LOGnOUT(4,<<"	Expectancy(Gain)/Expectancy(Loss)  ratio by Gamma Params= "<< (bestGainAlpha/bestGainBeta)/(bestLossAlpha/bestLossBeta)<<endl);
	//LOGnOUT(4,<<"	Expectancy(Gain/Loss) ratio by computation = "<< computeExpectationOfGainLossRatio(gainDist, lossDist)<<endl);
	LOGnOUT(4,<<"	Expectancy(Gain)/Expectancy(Loss)  by computation = "<< computeExpOfGainByExpOfLossRatio(gainDist, lossDist)<<endl);	
	
	if (gainLossOptions::_isRootFreqEQstationary){
		MDOUBLE estimatedStationaryFreq = computeExpectationOfStationaryFrequency(gainDist,lossDist);
		LOGnOUT(4,<<" Stationary '1' Freq at the root - for each stochastic process g/(g+l), with expectation of "<<estimatedStationaryFreq<<endl);}
	else{
		switch (gainLossOptions::_characterFreqEval){
				case (gainLossOptions::FiftyFifty):  
					LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<static_cast<gainLossModel*>((*spVVec[0][0]).getPijAccelerator()->getReplacementModel())->getTheta()<<endl);
					break;
				case (gainLossOptions::LeavesAve):
					LOGnOUT(4,<<"frequencies are based on LeavesAve (-F option) "<<static_cast<gainLossModel*>((*spVVec[0][0]).getPijAccelerator()->getReplacementModel())->getTheta()<<endl);
					break;
				case (gainLossOptions::optimizeOverTree):
					LOGnOUT(4,<<"Theta "<<static_cast<gainLossModel*>((*spVVec[0][0]).getPijAccelerator()->getReplacementModel())->getTheta() <<endl);
					break;
		}
	}
	LOGnOUT(4,<<" Total branch lengths:"<<tr.getAllBranchesLengthSum() <<endl);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::printModellValuesOfParams(stochasticProcess* sp, tree& tr)
{
	if(isAlphaOptimization(sp->distr()))LOGnOUT(4,<<" AlphaRate "<<getRateAlpha(sp->distr()) <<endl);
	if(isBetaOptimization(sp->distr()))LOGnOUT(4,<<" BetaRate "<<getRateBeta(sp->distr()) <<endl);
	MDOUBLE gain = static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->getMu1();
	LOGnOUT(4,<<" Gain "<<gain <<endl);
	MDOUBLE loss;
	if(!gainLossOptions::_isReversible){
		loss = static_cast<gainLossModelNonReversible*>((*sp).getPijAccelerator()->getReplacementModel())->getMu2();				
		LOGnOUT(4,<<" Loss "<<loss<<endl);
		LOGnOUT(4,<<"	Gain/Loss  ratio= "<< gain/loss<<endl);
	}
	else{
		loss = static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->getMu2();
	}
	if (gainLossOptions::_isRootFreqEQstationary){
		LOGnOUT(4,<<" Stationary '1' Freq at the root (g/(g+l) = "<<gain/(gain+loss) <<endl);}
	else{
		switch (gainLossOptions::_characterFreqEval){
				case (gainLossOptions::FiftyFifty):  
					LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->getTheta()<<endl);
					break;
				case (gainLossOptions::LeavesAve):
					LOGnOUT(4,<<"frequencies are based on LeavesAve (-F option) "<<static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->getTheta()<<endl);
					break;
				case (gainLossOptions::optimizeOverTree):
					LOGnOUT(4,<<" Theta ('1' at root):"<<static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->getTheta() <<endl);
					break;
		}
	}	
	LOGnOUT(4,<<" Total branch lengths:"<<tr.getAllBranchesLengthSum() <<endl);
}




/********************************************************************************************
 *********************************************************************************************/
void gainLoss::optimizationsManyStarts(const MDOUBLE epsilonOptimization, const int numIterations){
	int bestModel=0;
	MDOUBLE epsilonOptimizationCorrected = min(epsilonOptimization, abs(_logL)*gainLossOptions::_percentOfImprovManySarts);
	LOGnOUT(4,<<"\n\n --- start optimizationsManyStarts for "<<gainLossOptions::_numberOfRandPointsInOptimization<<" rand points, with epsilonIteration "<<epsilonOptimizationCorrected<<endl);

	Vdouble likeVecOpt;
	likeVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);
	vector<stochasticProcess*> spVecOpt;
	spVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);
	vector<tree> trVecOpt;
	trVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);	

	for(int i=0; i<gainLossOptions::_numberOfRandPointsInOptimization; ++i){
		LOGnOUT(4,<<"\n\n-------startOptimization "<<i+1<<endl);
		stochasticProcess* sp = _sp->clone();
		tree tr = _tr;
		unObservableData* currUnObs;
		if(_unObservableData_p)
			currUnObs = _unObservableData_p->clone();
		else
			currUnObs = NULL;

		// initialize
		initParamsAtRandPoints(gainLossOptions::_numberOfRandStartPoints,sp,currUnObs);
		// optimize
		//cout<<"before: "<<currUnObs->getlogLforMissingData()<<endl;
		bool  isbblLSWhenbblEMdontImprove = false;
		gainLossOptimizer glOpt(tr,sp,_scUniqPatterns,
			epsilonOptimizationCorrected,numIterations,
			epsilonOptimizationCorrected*gainLossOptions::_epsilonFactor_Model,
			(int)floor(numIterations*gainLossOptions::_numIterationsFactor_Model),
			epsilonOptimizationCorrected*gainLossOptions::_epsilonFactor_BBL,
			(int)floor(numIterations*gainLossOptions::_numIterationsFactor_BBL),
			_weightsUniqPatterns,
			currUnObs,(bool)Parameters::getInt("_performOptimizationsBBLManyStarts"), isbblLSWhenbblEMdontImprove);
		//cout<<"after: "<<currUnObs->getlogLforMissingData()<<endl;

		tr = glOpt.getOptTree();
		spVecOpt[i]=sp;
		trVecOpt[i]=tr;
		likeVecOpt[i]=glOpt.getBestL();
//if(currUnObs){
//	currUnObs->setLforMissingData(tr,sp);				
//}
		MDOUBLE estL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,_scUniqPatterns,*sp,_weightsUniqPatterns,currUnObs);
		if(!DEQUAL(likeVecOpt[i],estL)){
			LOGnOUT(3,<<" --- error: different likelihood after optimizeGainLossModel,diff= "<<likeVecOpt[i]-estL <<"\n");
		}
		if(likeVecOpt[i]>likeVecOpt[bestModel])
			bestModel = i;
		LOGnOUT(4,<<"-------L= "<<likeVecOpt[i]<<endl);
		if(currUnObs)	delete currUnObs;
	}
	_sp = spVecOpt[bestModel];
	_tr = trVecOpt[bestModel];
	if(Parameters::getInt("_accountForMissingData"))
		_unObservableData_p->setLforMissingData(_tr,_sp);

	LOGnOUT(4,<<" --- likelihood of All models: "<<endl);
	for(int i=0; i<gainLossOptions::_numberOfRandPointsInOptimization; ++i){
		LOGnOUT(4,<<"likelihood of model "<<i+1<<"\t"<<likeVecOpt[i]<<endl);
		if((spVecOpt[i]) && (i!=bestModel)) {delete spVecOpt[i];}
	}
	LOGnOUT(4,<<"likelihood of Best model "<<bestModel+1<<"\t"<<likeVecOpt[bestModel]<<endl);
	MDOUBLE lll = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*_sp,_weightsUniqPatterns,_unObservableData_p);
	if(!DEQUAL(likeVecOpt[bestModel],lll)){
		LOGnOUT(3,<<"ERROR: re-computed likelihood is diff by- "<<likeVecOpt[bestModel]<<lll<<endl);
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::optimizationsVVManyStarts(const MDOUBLE epsilonOptimization, const int numIterations){
	int bestModel=0;
	Vdouble likeVecOpt;
	MDOUBLE epsilonOptimizationCorrected = min(epsilonOptimization, abs(_logL)*gainLossOptions::_percentOfImprovManySarts);
	LOGnOUT(4,<<"\n\n --- start optimizationsVVManyStarts for "<<gainLossOptions::_numberOfRandPointsInOptimization<<" rand points, with epsilonIteration "<<epsilonOptimizationCorrected<<endl);

	likeVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);
	vector<vector<vector<stochasticProcess*> > > spVVVecOpt;
	spVVVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);
	vector<distribution*> gainDistVecOpt;
	gainDistVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);			
	vector<distribution*> lossDistVecOpt;
	lossDistVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);	
	vector<tree> trVecOpt;
	trVecOpt.resize(gainLossOptions::_numberOfRandPointsInOptimization);	

	for(int i=0; i<gainLossOptions::_numberOfRandPointsInOptimization; ++i){
		LOGnOUT(4,<<"\n\n-------startOptimization "<<i+1<<endl);
			tree tr = _tr;
		distribution* gainDist =_gainDist->clone();
		distribution* lossDist =_lossDist->clone();
		vector<vector<stochasticProcess*> >  spVVec;
		spVVec.resize(_gainDist->categories());
		for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
			spVVec[gainCategor].resize(_lossDist->categories());
			for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){		
				spVVec[gainCategor][lossCategor] = _spVVec[gainCategor][lossCategor]->clone();
			}
		}
		//stochasticProcess* sp = _sp->clone();
		unObservableData* currUnObs;
		if(_unObservableData_p)
			currUnObs = _unObservableData_p->clone();
		else
			currUnObs = NULL;
		//initialize random
		initParamsAtRandPointsSPvv(gainLossOptions::_numberOfRandStartPoints,spVVec,gainDist,lossDist,currUnObs);
		bool  isbblLSWhenbblEMdontImprove = false;
		gainLossOptimizer glOpt(tr,spVVec,gainDist,lossDist,_scUniqPatterns,
			epsilonOptimizationCorrected,numIterations,
			epsilonOptimizationCorrected*gainLossOptions::_epsilonFactor_Model,
			(int)floor(numIterations*gainLossOptions::_numIterationsFactor_Model),
			epsilonOptimizationCorrected*gainLossOptions::_epsilonFactor_BBL,
			(int)floor(numIterations*gainLossOptions::_numIterationsFactor_BBL),
			_weightsUniqPatterns, currUnObs,(bool)Parameters::getInt("_performOptimizationsBBLManyStarts"), isbblLSWhenbblEMdontImprove);
		tr = glOpt.getOptTree();
		spVVVecOpt[i]=spVVec;
		gainDistVecOpt[i]=gainDist;
		lossDistVecOpt[i]=lossDist;
		trVecOpt[i]=tr;
		likeVecOpt[i]=glOpt.getBestL();
		if(likeVecOpt[i]>likeVecOpt[bestModel])
			bestModel = i;
		LOGnOUT(4,<<"-------L= "<<likeVecOpt[i]<<endl);
		if(currUnObs)	delete currUnObs;
	}

	_spVVec = spVVVecOpt[bestModel];
	_gainDist = gainDistVecOpt[bestModel];
	_lossDist = lossDistVecOpt[bestModel];
	_tr = trVecOpt[bestModel];
	if(_unObservableData_p)	_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);

	for(int i=0; i<gainLossOptions::_numberOfRandPointsInOptimization; ++i){
		LOGnOUT(4,<<"likelihood of model "<<i+1<<"\t"<<likeVecOpt[i]<<endl);
		if((gainDistVecOpt[i]) && (i!=bestModel)) {delete gainDistVecOpt[i];}
		if((lossDistVecOpt[i]) && (i!=bestModel)) {delete lossDistVecOpt[i];}
		if((spVVVecOpt[i][0][0]) && (i!=bestModel)){
			for (int gainCategor=0; gainCategor<_gainDist->categories(); gainCategor++){
				for (int lossCategor=0; lossCategor<_lossDist->categories(); lossCategor++){		
					delete spVVVecOpt[i][gainCategor][lossCategor];
				}
			}
		}
	}
	LOGnOUT(4,<<"likelihood of Best model "<<bestModel+1<<"\t"<<likeVecOpt[bestModel]<<endl);

}

/********************************************************************************************
initParamsAtIntervalPoints
*********************************************************************************************/
//void gainLoss::initParamsAtIntervalPoints(int pointIndex,int numOfPoints, stochasticProcess* sp, unObservableData* currUnObs, ostream& out){
//	int numberOfParameters = 1;
//	bool optimizeAlpha = isAlphaOptimization(sp->distr());
//	bool optimizeBeta = isBetaOptimization(sp->distr());
//	bool optimizeMixture = isMixOptimization(sp->distr());
//	bool probInvariant = isInvariantOptimization(sp->distr());
//	bool evalTheta = isThetaOptimization();
//	if(optimizeAlpha)
//		++numberOfParameters;
//	if(optimizeBeta)
//		++numberOfParameters;
//	if(evalTheta)
//		++numberOfParameters;
//	if(probInvariant)
//		++numberOfParameters;
//	if(optimizeMixture)
//		++numberOfParameters;
//	if (!gainLossOptions::_isReversible)
//		++numberOfParameters;
//
//	MDOUBLE numOfPointsPerParam = (MDOUBLE)numOfPoints/numberOfParameters;
//
//
//}



/********************************************************************************************
initParamsAtRandPoints
*********************************************************************************************/
void gainLoss::initParamsAtRandPoints(int numOfRandPoints, stochasticProcess* sp, unObservableData* currUnObs, ostream& out){
	time_t t1;
	time(&t1);
	time_t t2;

	LOGnOUT(4,<<"Starting initParamsAtRandPoints with: numOfRandPoints="<<numOfRandPoints<<endl);	
	bool optimizeAlpha = isAlphaOptimization(sp->distr());
	bool optimizeBeta = isBetaOptimization(sp->distr());
	//bool optimizeMixture = isMixOptimization(sp->distr());
	bool probInvariant = isInvariantOptimization(sp->distr());
	bool evalTheta = isThetaOptimization();		
	
	MDOUBLE bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*_sp,_weightsUniqPatterns,currUnObs);
	MDOUBLE bestM1 =1; 
	MDOUBLE bestM2 =1; 
	MDOUBLE bestAlpha =1;
	MDOUBLE bestBeta =1;
	MDOUBLE bestTheta =0.5;
	MDOUBLE bestprobInvariantRate =0.05;
	bool isImprovedRandPoint = false;

	MDOUBLE L =VERYSMALL;
	MDOUBLE currM1;  
	MDOUBLE currM2; 
	MDOUBLE currAlpha;
	MDOUBLE currBeta; 
	MDOUBLE currTheta;
	MDOUBLE currprobInvariantRate;
	int i;
	for (i = 0; i < numOfRandPoints ; ++i) {
		 currM1 =talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userGainMin, gainLossOptions::_userGainMax);  
		 if (!gainLossOptions::_isReversible) currM2=talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userLossMin, gainLossOptions::_userLossMax); 
		 if(optimizeAlpha) currAlpha=talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userAlphaRateMin, gainLossOptions::_userAlphaRateMax);
		 if(optimizeBeta) currBeta=talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userBetaRateMin, gainLossOptions::_userBetaRateMax); 
		 if(evalTheta) currTheta=talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userThetaMin, gainLossOptions::_userThetaMax);
		 if(probInvariant) currprobInvariantRate =talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userProbInvariantRateMin, gainLossOptions::_userProbInvariantRateMax);
		static_cast<gainLossModel*>(sp->getPijAccelerator()->getReplacementModel())->setMu1(currM1, gainLossOptions::_isReversible);
		if (!gainLossOptions::_isReversible){
			static_cast<gainLossModelNonReversible*>(sp->getPijAccelerator()->getReplacementModel())->setMu2(currM2);	}
		if(optimizeAlpha){
			setRateAlpha(sp->distr(),currAlpha);;	}
		if(optimizeBeta){
			setRateBeta(sp->distr(),currBeta);	}
		if(evalTheta){
			static_cast<gainLossModel*>(sp->getPijAccelerator()->getReplacementModel())->setTheta(currTheta);	}
		if(probInvariant){
			static_cast<generalGammaDistributionPlusInvariant*>(sp->distr())->setInvProb(currprobInvariantRate);	}

// compute Likelihood
		MDOUBLE sumPijQij = normalizeQ(sp);
		if(currUnObs) currUnObs->setLforMissingData(_tr,sp);
		L = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*sp,_weightsUniqPatterns,currUnObs);

//print
		LOG(7,<<"--paramsSet: "<<i<<endl);
		LOG(7,<<"mu1="<<currM1<<endl);
		if (!gainLossOptions::_isReversible)	LOG(7,<<"mu2="<<currM2<<endl);
		if(optimizeAlpha)						LOG(7,<<"AlphaRate="<<currAlpha<<endl);
		if(optimizeBeta)						LOG(7,<<"AlphaBeta="<<currBeta<<endl);
		LOG(7,<<"likelihood is "<<L<<endl);

		
		if(bestL < L){
			bestM1 = currM1; 
			if (!gainLossOptions::_isReversible)bestM2 = currM2; 
			if(optimizeAlpha)bestAlpha = currAlpha;
			if(optimizeBeta)bestBeta = currBeta;
			if(evalTheta)bestTheta = currTheta;
			if(probInvariant) bestprobInvariantRate = currprobInvariantRate;
			bestL = L;
			isImprovedRandPoint = true;
			LOGnOUT(4,<<" logL improved = "<<L<<" at rand point "<<i<<endl);
		}
		if(isImprovedRandPoint && i>10)	// The loop is break after improvement and numOfRandPoints/2
			break;

	}
	// set best params after all rand points were calculated
	static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->setMu1(bestM1,gainLossOptions::_isReversible);	
	if (!gainLossOptions::_isReversible)
		static_cast<gainLossModelNonReversible*>((*sp).getPijAccelerator()->getReplacementModel())->setMu2(bestM2);
	if(optimizeAlpha)
		setRateAlpha((*sp).distr(),bestAlpha);
	if(optimizeBeta)
		setRateBeta((*sp).distr(),bestBeta);
	if(evalTheta)
		static_cast<gainLossModel*>((*sp).getPijAccelerator()->getReplacementModel())->setTheta(bestTheta);
	if(probInvariant)
		static_cast<generalGammaDistributionPlusInvariant*>(sp->distr())->setInvProb(bestprobInvariantRate);
	if(currUnObs) currUnObs->setLforMissingData(_tr,sp);
	L = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*sp,_weightsUniqPatterns,currUnObs);
	time(&t2);
	LOGnOUT(4,<<"End initParamsAtRandPoints after "<<i<<" randPoints with:\n bestGain "<<bestM1<<"\n bestLoss "<<bestM2<<
		"\n bestAlpha "<<bestAlpha<<"\n bestBeta "<<bestBeta <<"\n bestTheta "<<bestTheta <<"\n bestlogL "<<L<<endl);
	LOGnOUT(4,<<"initParamsAtRandPoints RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::initParamsAtRandPointsSPvv(int numOfRandPoints, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist,unObservableData* currUnObs, ostream& out){
	time_t t1;
	time(&t1);
	time_t t2;
	
	LOGnOUT(4,<<"Starting initParamsAtRandPointsSPvv with: numOfRandPoints="<<numOfRandPoints<<endl);
	MDOUBLE bestL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,_spVVec, _gainDist,_lossDist,_weightsUniqPatterns,currUnObs);
	stochasticProcess* sp = spVVec[0][0];

	MDOUBLE bestGainAlpha=1;	//Gain
	MDOUBLE bestGainBeta=1;
	MDOUBLE bestGainProbInvariant = 0.5;
	MDOUBLE bestLossAlpha=1; // Loss (for non-reversible model only)
	MDOUBLE bestLossBeta=1; 
	MDOUBLE bestLossProbInvariant = 0.5;
	MDOUBLE bestAlpha=1;	//Rate
	MDOUBLE bestTheta = 0.5;
	bool isImprovedRandPoint = false;

	MDOUBLE L;
	MDOUBLE currGainAlpha;	//Gain
	MDOUBLE currGainBeta;
	MDOUBLE currGainProbInvariant;
	MDOUBLE currLossAlpha; // Loss (for non-reversible model only)
	MDOUBLE currLossBeta; 
	MDOUBLE currLossProbInvariant;
	MDOUBLE currAlpha;	//Rate
	MDOUBLE currTheta;	
	
	bool optimizeAlpha = isAlphaOptimization((sp->distr()));
	bool optimizeBetaGain = isBetaOptimization(gainDist);
	bool optimizeBetaLoss = isBetaOptimization(lossDist);
	bool probInvariant = isInvariantOptimization(gainDist);	//for both
	bool evalTheta = isThetaOptimization();	

	int i;
	for (i = 0; i < numOfRandPoints ; ++i) {
//rand make
		currGainAlpha = talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userAlphaGainMin, gainLossOptions::_userAlphaGainMax); 
		if(optimizeBetaGain) currGainBeta = talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userBetaGainMin, gainLossOptions::_userBetaGainMax);
		if(probInvariant) currGainProbInvariant = talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userProbInvariantGainMin, gainLossOptions::_userProbInvariantGainMax);
		if (!gainLossOptions::_isReversible){
			currLossAlpha= talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userAlphaLossMin, gainLossOptions::_userAlphaLossMax);// Loss (for non-reversible model only)
			if(optimizeBetaLoss) currLossBeta= talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userBetaLossMin, gainLossOptions::_userBetaLossMax); 
			if(probInvariant) currLossProbInvariant= talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userProbInvariantLossMin, gainLossOptions::_userProbInvariantLossMax);
		}
		if(optimizeAlpha)
			currAlpha = talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userAlphaRateMin, gainLossOptions::_userAlphaRateMax);
		if(evalTheta)
			//currTheta = talRandom::giveRandomNumberBetweenTwoPoints(max(0.0,(currTheta-0.1)), min(1.0,(currTheta+0.1)));
			currTheta = talRandom::giveRandomNumberBetweenTwoPoints(gainLossOptions::_userThetaMin, gainLossOptions::_userThetaMax);

//set params
		updateGainAlpha(currGainAlpha,spVVec,gainDist,lossDist);
		if(optimizeBetaGain)			
			updateGainBeta(currGainBeta,spVVec,gainDist,lossDist);			
		if(probInvariant)
			updateGainProbInvariant(currGainProbInvariant,gainDist);			
		if (!gainLossOptions::_isReversible){
			updateLossAlpha(currLossAlpha,spVVec,gainDist,lossDist);
			if(optimizeBetaLoss)
				updateLossBeta(currLossBeta,spVVec,gainDist,lossDist);
			if(probInvariant)
				updateLossProbInvariant(currLossProbInvariant,lossDist);			
		}
		if(optimizeAlpha)
			updateRateAlpha(currAlpha,spVVec,gainDist,lossDist);
		if (evalTheta)
			updateTheta(currTheta,spVVec,gainDist,lossDist);		
		
// compute Likelihood
		MDOUBLE sumPijQij = normalizeQ(spVVec,gainDist,lossDist);
		if(currUnObs)	currUnObs->setLforMissingData(_tr,spVVec,gainDist,lossDist);
		L = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,spVVec, gainDist,lossDist,_weightsUniqPatterns,currUnObs);

		//print
		LOG(7,<<"--paramsSet: "<<i<<endl);
		LOG(7,<<"GainAlpha is "<<currGainAlpha<<endl);
		if(optimizeBetaGain) LOG(7,<<"GainBeta is "<<currGainBeta<<endl);
		if(probInvariant) LOG(7,<<"GainProbInvariant is "<<currGainProbInvariant<<endl);
		if (!gainLossOptions::_isReversible){
			LOG(7,<<"LossAlpha is "<<currLossAlpha<<endl);
			if(optimizeBetaLoss) LOG(7,<<"LossBeta is "<<currLossBeta<<endl);
			if(probInvariant) LOG(7,<<"LossProbInvariant is "<<currLossProbInvariant<<endl);
		}
		if(optimizeAlpha)	LOG(7,<<"Alpha is "<<currAlpha<<endl);	
		if(evalTheta)	LOG(7,<<"Theta is "<<currTheta<<endl);
		LOG(7,<<"likelihood is "<<L<<endl);

		
		if(bestL < L){
			bestGainAlpha=currGainAlpha;	
			if(optimizeBetaGain) bestGainBeta=currGainBeta;
			if(probInvariant) bestGainProbInvariant = currGainProbInvariant;
			if (!gainLossOptions::_isReversible){
				bestLossAlpha=currLossAlpha; 
				if(optimizeBetaLoss) bestLossBeta=currLossBeta; 
				if(probInvariant) bestLossProbInvariant = currLossProbInvariant;
			}
			if(optimizeAlpha)	bestAlpha=currAlpha;	
			if(evalTheta)	bestTheta = currTheta;
			bestL = L;
			isImprovedRandPoint = true;
			LOGnOUT(4,<<" logL improved = "<<L<<" at rand point "<<i<<endl);
		}
		if(isImprovedRandPoint && i>10)	// The loop is break after improvement and numOfRandPoints/2
			break;

	}
// set best params after all rand points were calculated
	updateGainAlpha(bestGainAlpha,spVVec,gainDist,lossDist);
	if(optimizeBetaGain) updateGainBeta(bestGainBeta,spVVec,gainDist,lossDist);
	if(probInvariant) updateGainProbInvariant(bestGainProbInvariant,gainDist);
	if (!gainLossOptions::_isReversible){
		updateLossAlpha(bestLossAlpha,spVVec,gainDist,lossDist);
		if(optimizeBetaLoss) updateLossBeta(bestLossBeta,spVVec,gainDist,lossDist);
		if(probInvariant) updateLossProbInvariant(bestLossProbInvariant,lossDist);
	}
	if(optimizeAlpha)
		updateRateAlpha(bestAlpha,spVVec,gainDist,lossDist);
	if(evalTheta)
		updateTheta(bestTheta,spVVec,gainDist,lossDist);
	if(currUnObs)	currUnObs->setLforMissingData(_tr,spVVec,gainDist,lossDist);
	L = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,spVVec, gainDist,lossDist,_weightsUniqPatterns,currUnObs);
	time(&t2);
	LOGnOUT(4,<<"End initParamsAtRandPointsSPvv after "<<i<<" randPoints with:\n bestGainAlpha "<<bestGainAlpha<<"\n bestGainBeta "<<bestGainBeta<<
		"\n bestLossAlpha "<<bestLossAlpha<<"\n bestLossBeta "<<bestLossBeta <<"\n bestTheta "<<bestTheta <<"\n bestlogL "<<L<<endl);
	LOGnOUT(4,<<"initParamsAtRandPointsSPvv RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
}





// computations
/********************************************************************************************
*********************************************************************************************/
void gainLoss::startRate4Site(sequenceContainer& sc, tree& tr, stochasticProcess* sp,   string& outDir, unObservableData* unObservableData_p)
{
	LOGnOUT(4,<<endl<<"Starting rate4site..."<<endl);
	time_t t1,t2;
	time(&t1);

	rate4siteGL r4s(sc,tr,sp,outDir, unObservableData_p);
	r4s.run();
	
	r4s.printRates();
	r4s.printRatesNormalized();

	_postProbPerCatPerPos = r4s.getLpostPerCat();
	_rates  = r4s.getRates();
	_normalizedRates = r4s.getNormalizedRates();
	
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}
/********************************************************************************************
*********************************************************************************************/
void gainLoss::startGainLoss4Site(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> > spVVec,distribution* gainDist,distribution* lossDist,
								  string& outDir,  unObservableData* unObservableData_p)
{
	LOGnOUT(4,<<endl<<"Starting gain4site and loss4site..."<<endl);
	time_t t1;
	time(&t1);
	time_t t2;

	gainLoss4site gl4s(sc,tr,spVVec,gainDist,lossDist,outDir,unObservableData_p);
	gl4s.computeGain4Site();
	gl4s.computeLoss4Site();
	gl4s.printGain4Site();
	gl4s.printLoss4Site();

	_postProbPerSpPerCatPerPos = gl4s.getLpostPerSpPerCat();

	time(&t2);
	LOGnOUT(4,<<"computeEB_EXP_GL4Site RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::computePosteriorExpectationOfChangeRunOnly()
{
	LOGnOUT(4,<<endl<<"Starting computePosteriorExpectationOfChangeRunOnly..."<<endl);
	time_t t1,t2;
	time(&t1);

	computeCountsGL* countsGL =NULL;

	if(!gainLossOptions::_gainLossDist){
		if(_postProbPerCatPerPos.size()==0  )
		{	
			//resizeMatrix(LpostPerCat,sp->categories(),sc.seqLen()) ;	// Not needed done with vector "=" sign later
			if(_sp->categories()>1){	// to fill LpostPerCat - run computeRate4site()
				LOGnOUT(4,<<endl<<"The required LpostPerCat is empty - run Rate4Site to compute."<<endl);
				rate4siteGL r4s(_sc,_tr,_sp,gainLossOptions::_outDir, _unObservableData_p);
				r4s.run();
				_postProbPerCatPerPos = r4s.getLpostPerCat();
			}
			else{
				_postProbPerCatPerPos.resize(1);
				_postProbPerCatPerPos[0].resize(_sc.seqLen());
				oneMatrix(_postProbPerCatPerPos);
			}
		}
		countsGL = new computeCountsGL(_sc,_tr,_sp,gainLossOptions::_outDir,_postProbPerCatPerPos, _distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
	}
	else{
		if(_postProbPerSpPerCatPerPos.size()==0  )
		{	
			LOGnOUT(4,<<endl<<"The required LpostPerSpPerCat is empty - run computeGain4Site to compute."<<endl);
			gainLoss4site gl4s(_sc,_tr,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_unObservableData_p);
			gl4s.computeGain4Site();
			//gl4s.computeLoss4Site();	// No need to run both
			_postProbPerSpPerCatPerPos = gl4s.getLpostPerSpPerCat();
		}
		countsGL = new computeCountsGL(_sc,_tr,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_postProbPerSpPerCatPerPos,_distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
	}
	countsGL->run();
	_jointProb_PosNodeXY = countsGL->getJointProb();
	if(countsGL) delete countsGL;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}	

/********************************************************************************************
Main version, if used with other data (not the gainLoss _sc,_tr,_sp,.. members) other version required
*********************************************************************************************/
void gainLoss::startComputePosteriorExpectationOfChange()
{
	LOGnOUT(4,<<endl<<"Starting calculePosteriorExpectationOfChange..."<<endl);
	time_t t1,t2;
	time(&t1);

	computeCountsGL* countsGL =NULL;
	rate4siteGL* r4s=NULL;
	gainLoss4site* gl4s=NULL;

	if(!gainLossOptions::_gainLossDist){
		if(_postProbPerCatPerPos.size()==0  )
		{	
			//resizeMatrix(LpostPerCat,sp->categories(),sc.seqLen()) ;	// Not needed done with vector "=" sign later
			if(_sp->categories()>1){	// to fill LpostPerCat - run computeRate4site()
				LOGnOUT(4,<<endl<<"The required LpostPerCat is empty - run Rate4Site to compute."<<endl);
				r4s = new rate4siteGL(_sc,_tr,_sp,gainLossOptions::_outDir, _unObservableData_p);
				r4s->run();
				_postProbPerCatPerPos = r4s->getLpostPerCat();
				if(r4s) delete r4s;
			}
			else{
				_postProbPerCatPerPos.resize(1);
				_postProbPerCatPerPos[0].resize(_sc.seqLen());
				oneMatrix(_postProbPerCatPerPos);
			}
		}
		countsGL = new computeCountsGL(_sc,_tr,_sp,gainLossOptions::_outDir,_postProbPerCatPerPos, _distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
	}
	else{
		if(_postProbPerSpPerCatPerPos.size()==0  )
		{	
			LOGnOUT(4,<<endl<<"The required LpostPerSpPerCat is empty - run computeGain4Site to compute."<<endl);
			gl4s = new gainLoss4site(_sc,_tr,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_unObservableData_p);
			gl4s->computeGain4Site();
			//gl4s.computeLoss4Site();	// No need to run both
			_postProbPerSpPerCatPerPos = gl4s->getLpostPerSpPerCat();
			if(gl4s) delete gl4s;
		}
		countsGL = new computeCountsGL(_sc,_tr,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_postProbPerSpPerCatPerPos,_distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
	}
	countsGL->run();	
	countsGL->printProbExp();						// Expectation and Probability PerPos
	countsGL->produceExpectationPerBranch();		// required before printExpectationPerBranch
	countsGL->printExpectationPerBranch();			// sum over all pos
	countsGL->updateTreeByGainLossExpectationPerBranch(_trGain,0,1);
	countsGL->updateTreeByGainLossExpectationPerBranch(_trLoss,1,0);
	countsGL->printProbabilityPerPosPerBranch();	// with probCutOff

	if(gainLossOptions::_isFewCutOffCounts)
		countsGL->printProbExpPerPosPerBranchFewCutOffs(gainLossOptions::_probCutOffPrintEvent);
	else
		countsGL->printProbExpPerPosPerBranch(gainLossOptions::_probCutOffPrintEvent,gainLossOptions::_probCutOffCounts);

	if(gainLossOptions::_printPropExpOfChangeFullData){
		MDOUBLE probCutOffPrintEvent = 0;	// if <0.05 results with a huge file
		countsGL->printProbExpPerPosPerBranch(probCutOffPrintEvent ,gainLossOptions::_probCutOffCounts);
	}
	if(gainLossOptions::_printExpPerPosPerBranchMatrix){
		countsGL->printExpPerPosPerBranchMatrix(0,1);
		countsGL->printExpPerPosPerBranchMatrix(1,0);
	}
	if(gainLossOptions::_printTreesWithExpectationValuesAsBP){
		countsGL->printTreesWithExpectationValuesAsBP();
	}
	if(gainLossOptions::_printTreesWithProbabilityValuesAsBP){
		countsGL->printTreesWithProbabilityValuesAsBP();
	}
	//if(gainLossOptions::_saveProbChanges_PosNodeXY){	// the computedProbChanges_PosNodeXY is saved to be used
		resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),_jointProb_PosNodeXY);
		_jointProb_PosNodeXY = countsGL->getJointProb();
	//}
	_SMPerPos = countsGL->get_expV();
	_expChanges_PosNodeXY = countsGL->getExpChanges();

	_gainPerPos =countsGL-> get_expV01();
	_lossPerPos = countsGL-> get_expV10();

	_meanGain = computeAverage(countsGL-> get_expV01());
	_meanLoss = computeAverage(countsGL-> get_expV10());
	_medianGain = computeMedian(countsGL-> get_expV01());
	_medianLoss = computeMedian(countsGL-> get_expV10());
	LOGnOUT(4,<<"Mean   values Gain="<<_meanGain<<"\tLoss="<<_meanLoss<<endl);
	LOGnOUT(4,<<"Median values Gain="<<_medianGain<<"\tLoss="<<_medianLoss<<endl<<endl);
	if(countsGL) delete countsGL;

	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes (calculePosteriorExpectationOfChange)"<<endl<<endl);
}

/********************************************************************************************

*********************************************************************************************/
void gainLoss::startComputeAmongSitesCorrelations()
{
	LOGnOUT(4,<<endl<<"Starting computeAmongSitesCorrelations..."<<endl);
	time_t t1,t2;
	time(&t1);

	_scEvolvingSites = _sc;
	if(gainLossOptions::_isCorrelationsBasedOnMaxParsimonyMapping){
		_expChanges_PosNodeXYSampledData = _MP_PosNodeXY;
		_gainPerPosCorr = _gainMPPerPos;
		_lossPerPosCorr = _lossMPPerPos;
	}
	else{
		_expChanges_PosNodeXYSampledData = _expChanges_PosNodeXY;
		_gainPerPosCorr = _gainPerPos;
		_lossPerPosCorr = _lossPerPos;
	}

	if(Parameters::getInt("_isUpdateminNumOfMPEvent2RemoveSimulatedPositions")){
		MDOUBLE minNumOfMPEvent2RemoveSimulatedPositions = Parameters::getFloat("_minNumOfMPEvent2RemoveSimulatedPositions");
		MDOUBLE addedMinNumOfMPEvent = (int)sqrt((double)_sc.numberOfSeqs())/5;
		if(addedMinNumOfMPEvent>0){
			Parameters::updateParameter("_minNumOfMPEvent2RemoveSimulatedPositions",double2string(addedMinNumOfMPEvent+minNumOfMPEvent2RemoveSimulatedPositions).c_str());
			LOGnOUT(4,<<"Update _minNumOfMPEvent2RemoveSimulatedPositions to  "<<addedMinNumOfMPEvent+minNumOfMPEvent2RemoveSimulatedPositions<<" from "<<minNumOfMPEvent2RemoveSimulatedPositions<<" with respect to "<<_sc.numberOfSeqs()<< " species"<<endl);
		}
	}
	if(gainLossOptions::_isRemoveSimulatedPositionsBasedOnMP){
		LOGnOUT(4,<<endl<<"Based on Maximum parsimony Remove positions from Real data (remove positions with no events)..."<<endl);		
		vector<int> posToRemove(_scEvolvingSites.seqLen(),false);
		MDOUBLE minExpT_MP = Parameters::getFloat("_minNumOfMPEvent2RemoveSimulatedPositions")/2;//  
		MDOUBLE Nmin = 0;			
		LOGnOUT(4,<<"min Number Of Max Parsimony Event to consider a Position is "<<minExpT_MP*2<<endl);
		int numOfRemovedPos=0;
		for (int pos = 0; pos<_sc.seqLen(); ++pos){
			MDOUBLE Gain = _MPPerPos[pos][0][1];
			MDOUBLE Loss = _MPPerPos[pos][1][0];		
			if(gainLossOptions::_isOnlyCorrelateWithBoth)
				Nmin = computeNminRforCorrelWithGainAndLoss(Gain,Loss);
			else
				Nmin = max(Gain,Loss); // thus, position are removed only if both their gain and loss values are below minT
			if(Nmin < minExpT_MP){
				posToRemove[pos] = true;
				_expChanges_PosNodeXYSampledData.erase(_expChanges_PosNodeXYSampledData.begin() + pos-numOfRemovedPos );
				numOfRemovedPos++;
			}else{
				_evolvingSites.push_back(pos);
				_numOfGapsTillSite.push_back(numOfRemovedPos);
			}
		}
		_scEvolvingSites.removePositions(posToRemove);
		LOGnOUT(4,<<"removed="<<numOfRemovedPos<<endl);
		int numOfSimulatedPositionsAboveMinRate = _scEvolvingSites.seqLen();
		LOGnOUT(4,<<"After remove numOfPositions="<<numOfSimulatedPositionsAboveMinRate<<endl);
	}else{
		for (int pos = 0; pos<_sc.seqLen(); ++pos){
			_evolvingSites.push_back(pos);
			_numOfGapsTillSite.push_back(0);
		}
	}	
	computeCorrelations* computeCorrel =NULL;	// required for correlation computation	
	if(_expChanges_PosNodeXYSampledData.size()==0)
		errorMsg::reportError("ERROR: Correlation request with empty mapping vector");
	else
		computeCorrel = new computeCorrelations(_tr,  gainLossOptions::_outDir, &_expChanges_PosNodeXYSampledData);
		

	if(gainLossOptions::_printComputedCorrelationsAllSites || gainLossOptions::_selectedSitesForCorrelation==""){
		LOGnOUT(4,<<"Correlate all sites (all-against-all, STRING style print)"<<endl);
		//readIntegersFromFileIntoVector(_selectedSites,_scEvolvingSites.seqLen(), 0, NULL);	// all sites in range
		_selectedSites = _evolvingSites;		
	}else{
		LOGnOUT(4,<<"printComputedCorrelations... read positions from "<<gainLossOptions::_selectedSitesForCorrelation<<endl);
		readIntegersFromFileIntoVector(_selectedSites,_sc.seqLen()-1, 0, &gainLossOptions::_selectedSitesForCorrelation, &_evolvingSites); // fill _selectedSites
		Vint numOfGapsTillSiteSelected; // fill numOfGapsTillSiteSelected => _numOfGapsTillSite
		for(int i=0;i<_selectedSites.size();++i){
			for(int j=0;j<_evolvingSites.size();++j){
				if(_selectedSites[i]==_evolvingSites[j])
					numOfGapsTillSiteSelected.push_back(_numOfGapsTillSite[j]);
			}
		}
		_numOfGapsTillSite = numOfGapsTillSiteSelected;
	}

	//bool correlationForZscore = false;
	//LOGnOUT(4,<<"Warning: isNormalizeForBranch is by branch length. correlationForZscore false by Default. Both with and withour branch"<<endl);
	computeCorrel->runComputeCorrelations(_selectedSites, _numOfGapsTillSite, gainLossOptions::_isNormalizeForBranchExpInCorrCompute);
	
	// required before print. Can't be done before - out of vec: _expChanges_PosNodeXYSampledData index
	if(gainLossOptions::_isPrintCorrelationsOfAllPairs_Corr)
		computeCorrel->printComputedCorrelations(_selectedSites,_evolvingSites, gainLossOptions::_isNormalizeForBranchExpInCorrCompute);

	//if(gainLossOptions::_performParametricBootstapCorrelation){	// later use these values to print rank according to simulations
		_correlationsPerSitePerPosVec = computeCorrel->getcorrelationPerSitePerPosVec();
		_correlationsPerSitePerPosVecSampledData = _correlationsPerSitePerPosVec;
	//}
	//else{ // else we'll print it later, while taking into account simulations
	//	computeCorrel->printComputedCorrelations(selectedSites, true/*, correlationForZscore*/);
	//}	
	if(computeCorrel) delete computeCorrel;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes (computeAmongSitesCorrelations)"<<endl<<endl);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::startComputePosteriorExpectationOfChange(sequenceContainer& sc, tree& tr, stochasticProcess* sp
														, VVdouble LpostPerCat, unObservableData* unObservableData_p, string& outDir,MDOUBLE distanceFromNearestOTUForRecent ,bool isUpdateMPPerPos)
{
	LOGnOUT(4,<<endl<<"Starting calculePosteriorExpectationOfChange..."<<endl);
	time_t t1,t2;
	time(&t1);
	
	computeCountsGL* countsGL;	
	if(LpostPerCat.size()==0  )
	{	
		LOGnOUT(4,<<endl<<"The required LpostPerCat is empty - run Rate4Site to compute."<<endl);
		//resizeMatrix(LpostPerCat,sp->categories(),sc.seqLen()) ;	// Not needed done with vector "=" sign later
		if(sp->categories()>1){	// to fill LpostPerCat - run computeRate4site()
			rate4siteGL r4s(sc,tr,sp,outDir, unObservableData_p);
			r4s.run();
			LpostPerCat = r4s.getLpostPerCat();
		}
		else{
			oneMatrix(LpostPerCat);
		}
	}
	countsGL = new computeCountsGL(sc,tr,sp,outDir,LpostPerCat,distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent

	countsGL->run();	
	countsGL->printProbExp();						// Expectation and Probability PerPos
	countsGL->produceExpectationPerBranch();		// required before printExpectationPerBranch
	countsGL->printExpectationPerBranch();			// sum over all pos
	countsGL->updateTreeByGainLossExpectationPerBranch(_trGain,0,1);
	countsGL->updateTreeByGainLossExpectationPerBranch(_trLoss,1,0);
	countsGL->printProbabilityPerPosPerBranch();	// with probCutOff
	if(gainLossOptions::_isFewCutOffCounts)
		countsGL->printProbExpPerPosPerBranchFewCutOffs(gainLossOptions::_probCutOffPrintEvent);
	else
		countsGL->printProbExpPerPosPerBranch(gainLossOptions::_probCutOffPrintEvent,gainLossOptions::_probCutOffCounts);

	if(gainLossOptions::_printPropExpOfChangeFullData){
		MDOUBLE probCutOffPrintEvent = 0.0;	// if <0.05 results with a huge file
		countsGL->printProbExpPerPosPerBranch(probCutOffPrintEvent ,gainLossOptions::_probCutOffCounts);
	}
	if(gainLossOptions::_printExpPerPosPerBranchMatrix){
		countsGL->printExpPerPosPerBranchMatrix(0,1);
		countsGL->printExpPerPosPerBranchMatrix(1,0);
	}
	if(gainLossOptions::_printTreesWithExpectationValuesAsBP){
		countsGL->printTreesWithExpectationValuesAsBP();
	}
	if(gainLossOptions::_printTreesWithProbabilityValuesAsBP){
		countsGL->printTreesWithProbabilityValuesAsBP();
	}
	if(isUpdateMPPerPos)
		_SMPerPos = countsGL->get_expV();
	if(countsGL) delete countsGL;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}
/********************************************************************************************
*********************************************************************************************/
void gainLoss::startComputePosteriorExpectationOfChange(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> >& spVVec, distribution* gainDist, distribution* lossDist
														, VVVdouble& LpostPerSpPerCat,unObservableData* unObservableData_p, string& outDir,MDOUBLE distanceFromNearestOTUForRecent,bool isUpdateMPPerPos)
{
	LOGnOUT(4,<<endl<<"Starting calculePosteriorExpectationOfChange..."<<endl);
	time_t t1,t2;
	time(&t1);
	computeCountsGL* countsGL;	
	if(LpostPerSpPerCat.size()==0  )
	{	
		LOGnOUT(4,<<endl<<"The required LpostPerSpPerCat is empty - run computeGain4Site to compute."<<endl);
		gainLoss4site gl4s(sc,tr,spVVec,gainDist,lossDist,outDir,unObservableData_p);
		gl4s.computeGain4Site();
		//gl4s.computeLoss4Site();	// No need to run both
		LpostPerSpPerCat = gl4s.getLpostPerSpPerCat();

	}
	countsGL = new computeCountsGL(sc,tr,spVVec,gainDist,lossDist,outDir,LpostPerSpPerCat,distanceFromNearestOTUForRecent);

	countsGL->run();	
	countsGL->printProbExp();						// Expectation and Probability PerPos
	countsGL->produceExpectationPerBranch();		// required before printExpectationPerBranch
	countsGL->printExpectationPerBranch();			// sum over all pos
	countsGL->updateTreeByGainLossExpectationPerBranch(_trGain,0,1);
	countsGL->updateTreeByGainLossExpectationPerBranch(_trLoss,1,0);
	countsGL->printProbabilityPerPosPerBranch();	// with probCutOff
	if(gainLossOptions::_isFewCutOffCounts)
		countsGL->printProbExpPerPosPerBranchFewCutOffs(gainLossOptions::_probCutOffPrintEvent);
	else
		countsGL->printProbExpPerPosPerBranch(gainLossOptions::_probCutOffPrintEvent,gainLossOptions::_probCutOffCounts);

	if(gainLossOptions::_printPropExpOfChangeFullData){
		MDOUBLE probCutOffPrintEvent = 0.0;	// if <0.05 results with a huge file
		countsGL->printProbExpPerPosPerBranch(probCutOffPrintEvent ,gainLossOptions::_probCutOffCounts);
	}
	if(gainLossOptions::_printExpPerPosPerBranchMatrix){
		countsGL->printExpPerPosPerBranchMatrix(0,1);
		countsGL->printExpPerPosPerBranchMatrix(1,0);
	}
	if(gainLossOptions::_printTreesWithExpectationValuesAsBP){
		countsGL->printTreesWithExpectationValuesAsBP();
	}
	if(gainLossOptions::_printTreesWithProbabilityValuesAsBP){
		countsGL->printTreesWithProbabilityValuesAsBP();
	}
	if(isUpdateMPPerPos)
		_SMPerPos = countsGL->get_expV();
	if(countsGL) delete countsGL;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}


/********************************************************************************************
A version used for simulated sequences
Good only for Gamma distribution
Unused!
*********************************************************************************************/
void gainLoss::computeCoEvolutionScoresBasedOnSimulatedData(sequenceContainer& scSimulated)
{
	LOGnOUT(4,<<endl<<" computeCoEvolutionScoresBasedOnSimulatedData..."<<endl);
	errorMsg::reportError("ERROR: computeCoEvolutionScoresBasedOnSimulatedData Not update. Check code or re-run");
	time_t t1,t2;
	time(&t1);
	VVdouble postProbPerCatPerPos;

	computeCountsGL* countsGL;	
	if(!gainLossOptions::_gainLossDist){
		if(_sp->categories()>1){	// to fill LpostPerCat - run computeRate4site()
			LOGnOUT(4,<<endl<<"The required LpostPerCat is empty - run Rate4Site to compute."<<endl);
			rate4siteGL r4s(scSimulated,_tr,_sp,gainLossOptions::_outDir, _unObservableData_p);
			r4s.run();
			postProbPerCatPerPos = r4s.getLpostPerCat();
		}
		else{
			postProbPerCatPerPos.resize(1);
			postProbPerCatPerPos[0].resize(scSimulated.seqLen());
			oneMatrix(postProbPerCatPerPos);
		}		
		countsGL = new computeCountsGL(scSimulated,_tr,_sp,gainLossOptions::_outDir,postProbPerCatPerPos, _distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
	}
	else{
		LOGnOUT(4,<<"ERROR - mixture model not supported in co-evolution by parametric bootstrap"<<endl);
	}
	countsGL->run();
	VVVVdouble expChanges_PosNodeXY_Sim = countsGL->getExpChanges();	// simulated data mapping


	////////// Correlations
	computeCorrelations* computeCorrel;
	Vint selectedSites;
	computeCorrel = new computeCorrelations(_tr,  gainLossOptions::_outDir, &_expChanges_PosNodeXY,  &expChanges_PosNodeXY_Sim);
	if(gainLossOptions::_printComputedCorrelationsAllSites || gainLossOptions::_selectedSitesForCorrelation==""){
		LOGnOUT(4,<<"Correlate all sites (all-against-all, STRING style print)"<<endl);
		readIntegersFromFileIntoVector(selectedSites,_sc.seqLen(), 0, NULL);	// all sites in range
	}else{
		LOGnOUT(4,<<" printComputedCorrelations... read positions from "<<gainLossOptions::_selectedSitesForCorrelation<<endl);
		readIntegersFromFileIntoVector(selectedSites,_sc.seqLen(), 0, &gainLossOptions::_selectedSitesForCorrelation);
	}
	LOGnOUT(4,<<"Warning: isNormalizeForBranch is by branch length. correlationForZscore false by Default. Both with and withour branch"<<endl);
	computeCorrel->runComputeCorrelations(selectedSites,_numOfGapsTillSite, gainLossOptions::_isNormalizeForBranchExpInCorrCompute);

	//computeCorrel->printComputedCorrelations(selectedSites, true/*, correlationForZscore*/); // DEB

	VVVdouble correlationsPerSitePerPosVecSim = computeCorrel->getcorrelationPerSitePerPosVec();	
	VVVdouble corPvalPerPos = _correlationsPerSitePerPosVec; //instead of resize
	computeCorrel->computedCorrelationsRankBasedOnSimulatedData(selectedSites, _correlationsPerSitePerPosVec,correlationsPerSitePerPosVecSim, corPvalPerPos);
	computeCorrel->produceSymeticMatrix(corPvalPerPos);
	bool correlationForZscore = false;
	computeCorrel->printComputedCorrelations(selectedSites,_evolvingSites, gainLossOptions::_isNormalizeForBranchExpInCorrCompute,correlationForZscore,&corPvalPerPos);

	if(countsGL) delete countsGL;
	if(computeCorrel) delete computeCorrel;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes (mapping+correlations with simulated data)"<<endl<<endl);
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::startParametricBootstapCorrelation(){
	LOGnOUT(4,<<"\n\n  startParametricBootstapCorrelation \n");

	gainLossAlphabet alph;
	sequenceContainer scSimulated;
	//_expChanges_PosNodeXYSampledData = _expChanges_PosNodeXY; // may be reduced, for BootStrap computation
	//_correlationsPerSitePerPosVecSampledData = _correlationsPerSitePerPosVec;

	
	if(Parameters::getInt("_isUpdateMinExpThresholdGivenRealDataQuantile") && Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair")>0){
		MDOUBLE gainQuantil = computeQuantileFrac(_gainPerPosCorr,gainLossOptions::_updateMinExpThresholdGivenRealDataQuantileVal);
		MDOUBLE lossQuantil = computeQuantileFrac(_lossPerPosCorr,gainLossOptions::_updateMinExpThresholdGivenRealDataQuantileVal);
		MDOUBLE qNminOfSimData = computeNminRforCorrelWithGainAndLoss(gainQuantil,lossQuantil);
		MDOUBLE minExpT = (double)Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
		if(minExpT < qNminOfSimData){
			Parameters::updateParameter("_minExpThresholdForPValComputationForCorrelatingPair",double2string(qNminOfSimData).c_str());
			LOGnOUT(4,<<"Update Nmin MinExpThreshold Given Read data quantile= "<<qNminOfSimData<<endl);
		}else{
			LOGnOUT(4,<<"No update Nmin Given Read data quantile= "<<qNminOfSimData<<" is smaller than current val="<<minExpT<<endl);
		}
	}

	MDOUBLE qNminOfRealData; // used for convergence
	if(_gainPerPosCorr.size()>3){
		MDOUBLE gainQuantil = computeQuantile(_gainPerPosCorr,gainLossOptions::_percentileOfNminWithCorr1RequiredForLastIteration);
		MDOUBLE lossQuantil = computeQuantile(_lossPerPosCorr,gainLossOptions::_percentileOfNminWithCorr1RequiredForLastIteration);
		qNminOfRealData = computeNminRforCorrelWithGainAndLoss(gainQuantil,lossQuantil);
	}else{
		qNminOfRealData = computeNminRforCorrelWithGainAndLoss(_meanGain,_meanLoss);
	}
	LOGnOUT(4,<<"\nStart Parametric bootstrap simulations. With up to "<<gainLossOptions::_numberOfIterations2simulate<<" iterations or till pair with Corr=1 simulated for Rate="<<qNminOfRealData<<"\n");

	MDOUBLE T_BH_prev = 0;
	string simCorrel = gainLossOptions::_outDir + "//" + "simCorrelationsFrequencies.txt";
	ofstream* simCorrelStream = new ofstream(simCorrel.c_str());
	int totalNumberOfSimulatedPairsAboveNmin= 0;
	bool isLastIteration = false;
	int numOfpairsWithRateAboveMinRequiredExp = (int)1E09; // temp init. the number of "hypothesis" tested, to be filled by CoMap
	
	////////////////////////////////////////////////////////////////////////////////
	int i =0;
	for(; i<gainLossOptions::_numberOfIterations2simulate && !isLastIteration ; ++i){
		LOGnOUT(4,<<"\n Parametric bootstrap iteration "<<i+1<<"\n");
		*simCorrelStream<<"iteration num  "<<i+1<<"\n";
		if(i==gainLossOptions::_numberOfIterations2simulate-1) // last, without convergence (median Bin with Corr=1)
			isLastIteration = true;		

		tree trSampled = _tr;
		sequenceContainer scSampled = _scEvolvingSites;
		if( gainLossOptions::_usePosSpecificSimulations){			
			startSimultePosteriorExpectationOfChange(gainLossOptions::_numberOfPositions2simulate,gainLossOptions::_numberOfSequences2simulate);
			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp1" + "//"+ "seqAll"+ "//" + "seq" + ".fa";
			ifstream in(strSeqNum.c_str());
			scSimulated = recognizeFormat::read(in,&alph);
		}
		else{
			scSimulated = simulateSequencesForParametricBootstrap(gainLossOptions::_numberOfPositions2simulate,scSampled,trSampled );
		}

		// sample leaves 2-remove. Not final. Need to sample from Real data and re-compute correlaion. Note: the _usePosSpecificSimulations is not updated
		bool isSample = false;
		if(isSample){
			MDOUBLE fractionOfSeq2Sample = 0.75;
			vector<int> seqIDs2remove;
			vector<tree::nodeP> nodes2remove;

			treeIterDownTopConst tIt(trSampled);
			for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
				if (mynode->isInternal()) 
					continue;
				MDOUBLE randV = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
				if(randV>fractionOfSeq2Sample && !mynode->father()->isRoot())
					nodes2remove.push_back(mynode);
			}
			LOGnOUT(3,<<"  In sampling, "<<nodes2remove.size()<<" sequences are removed"<<endl);

			for(int node=0; i<nodes2remove.size(); ++node){
				cout<<nodes2remove[node]->name()<<"\n";
				if(nodes2remove[node]->name()=="A")
					cout<<nodes2remove[node]->name()<<"\n";

				trSampled.removeLeaf(nodes2remove[node]);
			}
			//sequenceContainer::constTaxaIterator myseq=scSimulated.constTaxaBegin();
			//for (;myseq != scSimulated.constTaxaEnd(); ++myseq){
			//	if(talRandom::giveRandomNumberBetweenZeroAndEntry(1.0)<fractionOfSeq2Sample)
			//		seqIDs2remove.push_back(myseq->id());
			//}
			//for(int i=0; i<scSimulated.numberOfSeqs(); ++i){
			//	if(talRandom::giveRandomNumberBetweenZeroAndEntry(1.0)<fractionOfSeq2Sample)
			//		seqIDs2remove.push_back(i);	
			//}
			//for(int i=0; i<seqIDs2remove.size(); ++i){
			//	scSimulated.remove(seqIDs2remove[i]);
			//}
			intersectNamesInTreeAndSequenceContainer(trSampled,scSimulated);

			// Write seq and tree (required for re-labeling IDs
			string strSeqNum = gainLossOptions::_outDir + "//" + "seqSim."+ int2string(i) + ".fa";
			ofstream seq_out(strSeqNum.c_str());
			fastaFormat::  write(seq_out,scSimulated);
			string treeSampled = gainLossOptions::_outDir + "//" + "TheTree." + int2string(i) + ".ph"; 
			ofstream treeStream(treeSampled.c_str());
			trSampled.output(treeStream);		

			// re-Read
			ifstream in(strSeqNum.c_str());
			scSimulated = recognizeFormat::read(in,&alph);
			trSampled= tree(treeSampled);
		}

		if(Parameters::getInt("_isRemoveSimulatedPositionsWithExpectedLowNminBasedOnOccur")){
			LOGnOUT(3,<<"  Remove simulated position with too low/high occur to save later computation time (quick and (very) dirty)"<<endl);
			int minNumOfOnes = (int)Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
			int minNumOfZeros = (int)Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
			bool isRemovePosNotWithinMinMax = true;
			checkMinNumOfOnesOrZeros(scSimulated,minNumOfOnes,minNumOfZeros, isRemovePosNotWithinMinMax);
		}		
		if(Parameters::getInt("_accountForMissingData")){
			int minNumOfOnes = Parameters::getInt("_minNumOfOnes");
			int minNumOfZeros = Parameters::getInt("_minNumOfZeros");
			bool isRemovePosNotWithinMinMax = true;
			checkMinNumOfOnesOrZeros(scSimulated,minNumOfOnes,minNumOfZeros, isRemovePosNotWithinMinMax);
		}

		// version of startComputePosteriorExpectationOfChange with simulated seq. 1) get mapping vectors for each simulated seq 2) run correlation with both A,B
		if(totalNumberOfSimulatedPairsAboveNmin > numOfpairsWithRateAboveMinRequiredExp*10000.0 && totalNumberOfSimulatedPairsAboveNmin>1E09){
			isLastIteration = true; // in case there are 10000 more simulated pairs then tested pairs, last.
			LOGnOUT(4,<<"\n Last iteration of simulations, with sufficient simulated pairs "<<totalNumberOfSimulatedPairsAboveNmin<<" compared with tested pairs "<<numOfpairsWithRateAboveMinRequiredExp<<endl);
		}
		if((i+1) % gainLossOptions::_numberOfIterationsForPrintResults == 0 ) 	isLastIteration = true; //only for print
		totalNumberOfSimulatedPairsAboveNmin += computeCoEvolutionScoresBasedOnSimulatedDataCoMap(scSimulated,trSampled ,qNminOfRealData, isLastIteration ,numOfpairsWithRateAboveMinRequiredExp, T_BH_prev, simCorrelStream);
		if((i+1) % gainLossOptions::_numberOfIterationsForPrintResults == 0 ) 	isLastIteration = false; //revert back

		if(totalNumberOfSimulatedPairsAboveNmin < numOfpairsWithRateAboveMinRequiredExp*1000.0){
			isLastIteration = false; // revert to false, in case it was changed to 'true' due to simulation of Corr=1
			LOGnOUT(4,<<"More iterations of simulations required, with too few simulated pairs "<<totalNumberOfSimulatedPairsAboveNmin<<" compared with tested pairs "<<numOfpairsWithRateAboveMinRequiredExp<<" after "<<i+1<<" iterations\n");
		}		
	}
	LOGnOUT(4,<<"total NumberOf Pairs simulated and AboveNmin ="<<totalNumberOfSimulatedPairsAboveNmin<<" after "<<i+1<<" iterations\n");
}



/********************************************************************************************
A version used for simulated sequences
Good only for Gamma distribution
*********************************************************************************************/
int gainLoss::computeCoEvolutionScoresBasedOnSimulatedDataCoMap(sequenceContainer& scSimulated, tree& trSampled,  MDOUBLE qNminOfRealData, bool& isLastIteration, int& numOfpairsWithRateAboveMinRequiredExp, MDOUBLE& T_BH_prev, ofstream* simCorrelStream)
{
	LOGnOUT(4,<<endl<<"Compute: Mapping + Correlation. CoMap Algorithm"<<endl);
	time_t t1,t2;
	time(&t1);
	VVdouble postProbPerCatPerPos;
	VVVdouble postProbPerSpPerCatPerPos;
	
	Vdouble rate4siteSim;
	Vdouble rate4siteReal;
	Vdouble gainSim;
	Vdouble lossSim;

	computeCountsGL* countsGL = NULL;
	sankoffReconstructGL* sankoffReconstructMP = NULL;
	rate4siteGL* r4s = NULL;
	gainLoss4site* gl4s = NULL;

	if(gainLossOptions::_isRemoveSimulatedPositionsBasedOnMP){
		LOGnOUT(4,<<endl<<"MaxParsimonyChange for simulated data (remove positions with no events)..."<<endl);
		createDir(gainLossOptions::_outDir, "MPsimulations");
		string dirMP =  gainLossOptions::_outDir + "//" + "MPsimulations";
		sankoffReconstructGL sankoffReconstructMP(scSimulated, trSampled, dirMP,gainLossOptions::_costMatrixGainLossRatio, _distanceFromNearestOTUForRecent);
		VVVdouble MPPerPos = sankoffReconstructMP.getMPPerPos();
		vector<int> posToRemove(scSimulated.seqLen(),false);
		MDOUBLE minExpT_MP =Parameters::getFloat("_minNumOfMPEvent2RemoveSimulatedPositions")/2; //Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair")/2.0;
		MDOUBLE Nmin = 0;
		LOGnOUT(4,<<"min Number Of Max Parsimony Event to consider a Position is "<<minExpT_MP*2<<endl);
		int numOfRemovedPos=0;
		for (int pos = 0; pos<scSimulated.seqLen(); ++pos){
			MDOUBLE Gain = MPPerPos[pos][0][1];
			MDOUBLE Loss = MPPerPos[pos][1][0];		
			if(gainLossOptions::_isOnlyCorrelateWithBoth)
				Nmin = computeNminRforCorrelWithGainAndLoss(Gain,Loss);
			else
				Nmin = max(Gain,Loss); // thus, position are removed only if both their gain and loss values are below minT
			if(Nmin < minExpT_MP){
				posToRemove[pos] = true;
				numOfRemovedPos++;
			}
		}	
		scSimulated.removePositions(posToRemove);
		LOGnOUT(4,<<"removed="<<numOfRemovedPos<<endl);
		int numOfSimulatedPositionsAboveMinRate = scSimulated.seqLen();
		LOGnOUT(4,<<"After remove numOfPositions="<<numOfSimulatedPositionsAboveMinRate<<endl);
	}

	MDOUBLE meanGain, meanLoss, medianGain, medianLoss;
	VVVVdouble expChanges_PosNodeXY_Sim;
	if(gainLossOptions::_isCorrelationsBasedOnMaxParsimonyMapping){
		sankoffReconstructMP = new sankoffReconstructGL(scSimulated, trSampled, gainLossOptions::_outDir,gainLossOptions::_costMatrixGainLossRatio, _distanceFromNearestOTUForRecent);
		expChanges_PosNodeXY_Sim = sankoffReconstructMP->getMPPerPosPerNode();
		gainSim = sankoffReconstructMP-> getGainMPPerPos();
		lossSim = sankoffReconstructMP-> getLossMPPerPos();
	}else{
		if(!gainLossOptions::_gainLossDist){
			if(_sp->categories()>1){	// to fill LpostPerCat - run computeRate4site()
				LOGnOUT(4,<<endl<<"The required LpostPerCat is empty - run Rate4Site to compute..."<<endl);
				r4s = new rate4siteGL(scSimulated,trSampled,_sp,gainLossOptions::_outDir, _unObservableData_p);
				r4s->run();
				postProbPerCatPerPos = r4s->getLpostPerCat();
				if(gainLossOptions::_isUseRateForSiteAsNminForCorrelations){
					rate4siteSim = r4s->getRates();
					rate4siteReal = _rates;
				}
				if(r4s) delete r4s;
			}
			else{
				postProbPerCatPerPos.resize(1);
				postProbPerCatPerPos[0].resize(scSimulated.seqLen());
				oneMatrix(postProbPerCatPerPos);
			}		
			countsGL = new computeCountsGL(scSimulated,trSampled,_sp,gainLossOptions::_outDir,postProbPerCatPerPos, _distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
		}
		else{
			LOGnOUT(4,<<endl<<"The required LpostPerSpPerCat is empty - run computeGain4Site to compute..."<<endl);
			gl4s = new gainLoss4site(scSimulated,trSampled,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,_unObservableData_p);
			gl4s->computeGain4Site();
			postProbPerSpPerCatPerPos = gl4s->getLpostPerSpPerCat();
			if(gl4s) delete gl4s;
			countsGL = new computeCountsGL(scSimulated,trSampled,_spVVec,_gainDist,_lossDist,gainLossOptions::_outDir,postProbPerSpPerCatPerPos,_distanceFromNearestOTUForRecent);	//_distanceFromRootForRecent
		}
		countsGL->run();
		expChanges_PosNodeXY_Sim = countsGL->getExpChanges();	// simulated data mapping
		gainSim = countsGL-> get_expV01();
		lossSim = countsGL-> get_expV10();
	}
	meanGain = computeAverage(gainSim);
	meanLoss = computeAverage(lossSim);
	medianGain = computeMedian(gainSim);
	medianLoss = computeMedian(lossSim);
	LOGnOUT(4,<<"Mean   values Gain="<<meanGain<<"\tLoss="<<meanLoss<<endl);
	LOGnOUT(4,<<"Median values Gain="<<medianGain<<"\tLoss="<<medianLoss<<endl<<endl);
	
	if(Parameters::getInt("_isUpdateMinExpThresholdGivenSimulaitonsQuantile") && Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair")>0){
		MDOUBLE quantileVal = 0.1;
		MDOUBLE gainQuantil = computeQuantileFrac(gainSim,quantileVal);
		MDOUBLE lossQuantil = computeQuantileFrac(lossSim,quantileVal);
		MDOUBLE qNminOfSimData = computeNminRforCorrelWithGainAndLoss(gainQuantil,lossQuantil);
		MDOUBLE qNminOfSimDataPrev = qNminOfSimData;		
		while( qNminOfSimData-qNminOfSimDataPrev < 0.1 ){
			qNminOfSimDataPrev = qNminOfSimData;
			gainQuantil = computeQuantileFrac(gainSim,quantileVal);
			lossQuantil = computeQuantileFrac(lossSim,quantileVal);
			qNminOfSimData = computeNminRforCorrelWithGainAndLoss(gainQuantil,lossQuantil);			
			quantileVal += 0.1;
		}	

		MDOUBLE minExpT = (double)Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
		if(minExpT < qNminOfSimData){
			Parameters::updateParameter("_minExpThresholdForPValComputationForCorrelatingPair",double2string(qNminOfSimData).c_str());
			LOGnOUT(4,<<"Update MinExpThreshold GivenSimulaitonsQuantile= "<<qNminOfSimData <<" with respect to simulation quantile ="<<quantileVal<<endl);
		}else{
			LOGnOUT(4,<<"No update MinExpThreshold GivenSimulaitonsQuantile= "<<qNminOfSimData<<" is smaller than current val="<<minExpT<<endl);
		}
		Parameters::updateParameter("_isUpdateMinExpThresholdGivenSimulaitonsQuantile","0"); // Done once
	}
	if(countsGL) delete countsGL;
	if(sankoffReconstructMP) delete sankoffReconstructMP;

	// can remove only positions that are below Threshold in all correlation types
	MDOUBLE minExpT = Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
	LOGnOUT(4,<<"Remove simulated positions below Nmin="<<minExpT<<endl);
	LOGnOUT(4,<<"Before remove numOfPositions="<<scSimulated.seqLen()<<endl);
	vector<int> posToRemove(scSimulated.seqLen(),false);
	int numOfRemovedPos = 0;
	VVVVdouble expChanges_PosNodeXY_SimFinal;
	VVVdouble expChanges_PosXY_Sim;	
	computeRateValPerPos(expChanges_PosNodeXY_Sim, expChanges_PosXY_Sim);
	for (int pos = 0; pos<scSimulated.seqLen(); ++pos){
		MDOUBLE Nmin = 0;			
		MDOUBLE Gain = expChanges_PosXY_Sim[pos][0][1];
		MDOUBLE Loss = expChanges_PosXY_Sim[pos][1][0];
		
		if(gainLossOptions::_isOnlyCorrelateWithBoth)
			Nmin = computeNminRforCorrelWithGainAndLoss(Gain,Loss);
		else
			Nmin = max(Gain,Loss); // thus, position are removed only if both their gain and loss values are below minT

		if(Nmin < minExpT){
			posToRemove[pos] = true;
			numOfRemovedPos++;
		}
		else
			expChanges_PosNodeXY_SimFinal.push_back(expChanges_PosNodeXY_Sim[pos]);

	}	
	scSimulated.removePositions(posToRemove);
	LOGnOUT(4,<<"removed="<<numOfRemovedPos<<endl);
	int numOfSimulatedPositionsAboveMinRate = scSimulated.seqLen();
	LOGnOUT(4,<<"After remove numOfPositions="<<numOfSimulatedPositionsAboveMinRate<<endl);


	//// Correlations with simulated data
	computeCorrelations* computeCorrel =NULL;
	Vint selectedSitesSim;
	Vint numOfGapsTillSite;
	computeCorrel = new computeCorrelations(trSampled,  gainLossOptions::_outDir, &expChanges_PosNodeXY_SimFinal);
	readIntegersFromFileIntoVector(selectedSitesSim, numOfSimulatedPositionsAboveMinRate-1, 0, NULL,NULL);	// all sites in range
	numOfGapsTillSite.resize(selectedSitesSim.size(),0);

	//LOGnOUT(4,<<"Warning: isNormalizeForBranch is by branch length. correlationForZscore false by Default. Both with and without branch"<<endl);
	// Compute correlations of Sim data
	computeCorrel->runComputeCorrelations(selectedSitesSim,numOfGapsTillSite, gainLossOptions::_isNormalizeForBranchExpInCorrCompute);

	// sort Corr vector of Sim data
	computeCorrel->produceSortedVectorsOfAllCorrelations(rate4siteSim);	// maybe of size=0		
	
	// produce Bins of Sim data	
	int numberOfHighCorrInSimulationOfMedianNminBin = computeCorrel->produceSortedVectorsOfCorrelationsBinedByRate(qNminOfRealData, simCorrelStream);
	if(numberOfHighCorrInSimulationOfMedianNminBin >= 1   && gainLossOptions::_percentileOfNminWithCorr1RequiredForLastIteration<100){ // use 100 for no "Corr=1 based convergence"
		isLastIteration = true; // convergence (median Bin with Corr=1)
		LOGnOUT(4,<<"\n Last iteration of simulations, reached 'convergence' - simulated "<<numberOfHighCorrInSimulationOfMedianNminBin<<" pairs with Corr~=1 for bin with Rate in bin of "<<qNminOfRealData<<endl);
	}
			
	// compute correlations between real data (input _correlationsPerSitePerPosVec, _expChanges_PosNodeXY)
	//  and simulated data (according to Rate bins, already part of object)
	//   fill corPvalPerPos and _correlationsData
	resizeMatrix( _isComputePairWithRateAboveNim, _correlationsPerSitePerPosVecSampledData[0].size(),_correlationsPerSitePerPosVecSampledData[0][0].size()); // _isComputePairWithRateAboveNim - bool vector
	VVVdouble corPvalPerPos = _correlationsPerSitePerPosVecSampledData; //instead of resize - the vector of all-against-all Correlation coefficient determines the size of the vector PVals
	if(_correlationsPerSitePerPosVecSampledData.size() == 0 || _expChanges_PosNodeXYSampledData.size()==0)
		errorMsg::reportError("Real data correlation and expectation data is missing, can't compute simulation-based pVal");
	else
		numOfpairsWithRateAboveMinRequiredExp = computeCorrel->computedCorrelationsPValBasedOnSimulatedDataCoMapBins(_correlationsPerSitePerPosVecSampledData,_isComputePairWithRateAboveNim,_expChanges_PosNodeXYSampledData,corPvalPerPos
																	, _correlationsData, rate4siteReal ,_selectedSites,_numOfGapsTillSite,_evolvingSites, isLastIteration); // fill corPvalPerPos

	if(isLastIteration){ // compute FDR and print results
		bool correlationForZscore = false;
		//Vint selectedSites;
		//readIntegersFromFileIntoVector(selectedSites,_sc.seqLen(), 0, NULL);	// all sites in range

		string printType = "pVal";	
		if(gainLossOptions::_isPrintCorrelationsOfAllPairs_pVal)
			computeCorrel->printComputedCorrelations(_selectedSites,_evolvingSites, gainLossOptions::_isNormalizeForBranchExpInCorrCompute,correlationForZscore,&corPvalPerPos,&printType);

		if(gainLossOptions::_isFDRcorrectionForPValInCorrelation && _correlationsData.size()>0){
			Vdouble T_BH(corPvalPerPos.size()); // to be filled, for each corr type		
			// FDR			
			if(gainLossOptions::_isComputeQVals){
				VVVdouble corQvalPerPos = computeCorrel-> pVals2qVals (corPvalPerPos,_correlationsData,_isComputePairWithRateAboveNim, T_BH, _selectedSites,_evolvingSites);
				string printType = "qVal";
				computeCorrel->printComputedCorrelations(_selectedSites,_evolvingSites, gainLossOptions::_isNormalizeForBranchExpInCorrCompute,correlationForZscore,&corQvalPerPos,&printType);
			}else
				computeCorrel-> pVals2qVals (corPvalPerPos,_correlationsData,_isComputePairWithRateAboveNim, T_BH, _selectedSites,_evolvingSites);

			computeCorrel->printComputedCorrelationsData(gainLossOptions::_isNormalizeForBranchExpInCorrCompute,correlationForZscore,_correlationsData, T_BH);
			
			if(gainLossOptions::_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH){
				Vdouble minPValForPrint(corPvalPerPos.size(),gainLossOptions::_pValueCutOffForBootStrap); // same non-FDR min pVal for all correlation types	
				computeCorrel->printComputedCorrelationsData(gainLossOptions::_isNormalizeForBranchExpInCorrCompute,correlationForZscore,_correlationsData, minPValForPrint,gainLossOptions::_isPrintAllPairsOfCorrelatedSitesIncludingPValsAboveBH);
			}

			// Convergence of BH, less than 0.01% change
			MDOUBLE T_BH_currentMinuslast = T_BH[0]-T_BH_prev;
			LOGnOUT(4,<<" Convergence of BH: current Minus last="<<T_BH_currentMinuslast<<endl);
			T_BH_prev = T_BH[0]; // updated

		}else{
			LOGnOUT(4,<<" Note: No computation of FDR corrData (pVal significant) size="<<_correlationsData.size()<<endl);
		}
	}	

	if(computeCorrel) delete computeCorrel;
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes (mapping+correlations with simulated data)"<<endl);
	int numOfSimulatedPairsAboveMinRate = (numOfSimulatedPositionsAboveMinRate*(1+numOfSimulatedPositionsAboveMinRate))/2 - numOfSimulatedPositionsAboveMinRate;
	return numOfSimulatedPairsAboveMinRate;
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::startMaxParsimonyChange(bool isUpdateMPPerPos)
{
	LOGnOUT(4,<<endl<<"Starting MaxParsimonyChange..."<<endl);
	time_t t1,t2;
	time(&t1);
	sankoffReconstructGL sankoffReconstructMP(_sc, _tr, gainLossOptions::_outDir,gainLossOptions::_costMatrixGainLossRatio, _distanceFromNearestOTUForRecent);
	if(isUpdateMPPerPos){
		_MPPerPos = sankoffReconstructMP.getMPPerPos();
		_MP_PosNodeXY = sankoffReconstructMP.getMPPerPosPerNode();
		_gainMPPerPos = sankoffReconstructMP.getGainMPPerPos();
		_lossMPPerPos = sankoffReconstructMP.getLossMPPerPos();
		_CostOfTreeMP = sankoffReconstructMP.getNumOfGainEvnetsMP() + sankoffReconstructMP.getNumOfLossEvnetsMP() ;
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::startMaxParsimonyChange(sequenceContainer& sc, tree& tr, string& outDir,MDOUBLE costMatrixGainLossRatio, MDOUBLE distanceFromNearestOTUForRecent,bool isUpdateMPPerPos)
{
	LOGnOUT(4,<<endl<<"Starting MaxParsimonyChange..."<<endl);
	time_t t1,t2;
	time(&t1);
	sankoffReconstructGL sankoffReconstructMP(sc, tr, outDir,costMatrixGainLossRatio, distanceFromNearestOTUForRecent);
	if(isUpdateMPPerPos)
		_MPPerPos = sankoffReconstructMP.getMPPerPos();
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}



/********************************************************************************************
ancestralReconstructStates
*********************************************************************************************/
void gainLoss::ancestralReconstructor()
{
	LOGnOUT(4,<<endl<<"Starting ancestralReconstructor..."<<endl);
	time_t t1,t2;
	time(&t1);
	
	if(_alphVecDist.size()>_sc.getAlphabet()->size() && _alphVecDist[_sc.getAlphabet()->size()]>0)
		LOGnOUT(2,<<"\nWARNING !!! : ancestralReconstruct is not fully functional with missing data.\n Assume missing data indicates absence (indels)."<<endl<<endl);

	VVint statesV; // the vector with the states of the nodes, to be filled with reconstructed states (max over joint)
	ancestralReconstructStates ancestralReconst(_tr,_sc,_sp);	// Per POS,CAT

// compute joint reconstruction	(?)
	VVVdouble upL;
	VVVint backtrack;
	VVVint transitionTypeCount;
	ancestralReconst.traverseUpML(upL, backtrack);
	Vdouble LofJointV = ancestralReconst.traverseDownML(upL, backtrack, transitionTypeCount);	
	statesV = ancestralReconst.getStates();

	treeIterDownTopConst tIt(_tr);	// iterator used by following loops	

// sum over positions - joint
	Vint statesSum; // the vector with the Sum states of the nodes (joint)
	statesSum.resize(_tr.getNodesNum());	
	for (int pos = 0; pos <_sc.seqLen(); ++pos)	{
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()){			
			statesSum[mynode->id()]+=statesV[pos][mynode->id()];	// Sum over positions
		}
	}

	string AncestralReonstructSum = gainLossOptions::_outDir + "//" + "AncestralReconstructSumJoint.txt"; 
	ofstream AncestralReonstructSumStream(AncestralReonstructSum.c_str());
	AncestralReonstructSumStream<<"Node"<<"\t"<<"Sum"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		AncestralReonstructSumStream<<mynode->name()<<"\t"<<statesSum[mynode->id()]<<endl;
	}

// printAncestralReconstructFullData (joint)
	if(gainLossOptions::_printAncestralReconstructFullData){
		string AncestralReonstruct = gainLossOptions::_outDir + "//" + "AncestralReonstruct.txt"; 
		ofstream AncestralReonstructStream(AncestralReonstruct.c_str());
		AncestralReonstructStream<<"POS"<<"\t"<<"Node"<<"\t"<<"State"<<endl;
		for (int pos = 0; pos <_sc.seqLen(); ++pos)	{
			for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()){			
				AncestralReonstructStream<<pos+1<<"\t"<<mynode->name()<<"\t"<<statesV[pos][mynode->id()]<<endl;
			}
		}
	}
// print Trees
	if(gainLossOptions::_printTreesWithAncestralReconstructAsBP){
		createDir(gainLossOptions::_outDir, "TreesWithAncestralReonstruct");
		for (int pos = 0; pos <_sc.seqLen(); ++pos){
			string strTreeNum = gainLossOptions::_outDir + "//" + "TreesWithAncestralReonstruct" + "//" + "TreeAncRec" + int2string(pos+1) + ".ph";
			ofstream tree_out(strTreeNum.c_str());
			printTreeStatesAsBPValues(tree_out,statesV[pos],_tr);
		}
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}


/********************************************************************************************
ancestralReconstructorBasedOnJoint (Posterior)
//  interested in a the set of all the hypothetical taxonomic unit (HTU) sequences (joint reconstruction) 
// as oppose to a specific HTU whose sequence we would like to estimate (marginal reconstruction).
*********************************************************************************************/
void gainLoss::ancestralReconstructorBasedOnJoint()
{
	LOGnOUT(4,<<endl<<"Starting ancestralReconstructorBasedOnJoint..."<<endl);
	time_t t1,t2;
	time(&t1);

	VVint statesV; // the vector with the states of the nodes, to be filled with reconstructed states (max over joint)
	VVVdouble ancestralProbsPerPosNodeState;	// the vector with the probabilities of the nodes states, to be filled with reconstructed states (posterior)

	ancestralReconstructStates ancestralReconst(_tr,_sc,_sp);	// Per POS,CAT

	// compute posterior reconstruction
	// Prob(N=x|Data) = sum{fatherState}[P(N=x, father(N)=y|D)]}	
	if(_jointProb_PosNodeXY.size()==0){
		computePosteriorExpectationOfChangeRunOnly();	// this phase will also fill _jointProb_PosNodeXY
	}
	ancestralReconst.computeAncestralPosterior(_jointProb_PosNodeXY);
	ancestralProbsPerPosNodeState = ancestralReconst.getAncestralProbs(); // VVVdouble[pos][node][state] ancestralProbsPerPosNodeState

	treeIterDownTopConst tIt(_tr);	// iterator used by following loops
	// printAncestralReconstructFullData (posterior)
	if(gainLossOptions::_printAncestralReconstructPosterior){
		string AncestralReonstructPosterior = gainLossOptions::_outDir + "//" + "AncestralReconstructPosterior.txt"; 
		ofstream AncestralReonstructPosteriorStream(AncestralReonstructPosterior.c_str());
		AncestralReonstructPosteriorStream.precision(PRECISION);		

		AncestralReonstructPosteriorStream<<"POS"<<"\t"<<"Node"<<"\t"<<"State"<<"\t"<<"Prob"<<endl;
		for (int pos = 0; pos <_sc.seqLen(); ++pos)	{
			for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()){		
				int state=1;
				AncestralReonstructPosteriorStream<<pos+1<<"\t"<<mynode->name();
				//for (int state = 0; state <_sp->alphabetSize(); ++state){ // only state=1 is printed
				AncestralReonstructPosteriorStream<<"\t"<<state<<"\t"<<ancestralProbsPerPosNodeState[pos][mynode->id()][state];
				//}
				AncestralReonstructPosteriorStream<<endl;
			}
		}
	}
	// sum over positions - posterior
	VVdouble probStatesSum; // the vector with the Sum probes of the nodes (posterior)
	resizeMatrix(probStatesSum,_tr.getNodesNum(),_sp->alphabetSize());	
	Vdouble probOnesSum; // the vector with the Sum ones probes of the nodes (posterior) - good for {0,1}
	probOnesSum.resize(_tr.getNodesNum());
	for (int pos = 0; pos <_sc.seqLen(); ++pos)	{
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()){
			for (int state = 0; state <_sp->alphabetSize(); ++state){
				probStatesSum[mynode->id()][state]+=ancestralProbsPerPosNodeState[pos][mynode->id()][state];	// Sum over positions
				probOnesSum[mynode->id()] +=ancestralProbsPerPosNodeState[pos][mynode->id()][state]*state;	// if state==0 Nothing added 
			}			
		}
	}
	// print Sum: Table, Tree
	string AncestralReonstructPosteriorSum = gainLossOptions::_outDir + "//" + "AncestralReconstructPosteriorSum.txt"; 
	ofstream AncestralReonstructPosteriorSumStream(AncestralReonstructPosteriorSum.c_str());
	AncestralReonstructPosteriorSumStream.precision(PRECISION);

	AncestralReonstructPosteriorSumStream<<"Node"<<"\t"<<"State"<<"\t"<<"ProbSum"<<"\t"<<"Father"<<"\t"<<"StateFather"<<"\t"<<"ProbSumFather"<<endl;
	int state=1;
	//for (int state = 0; state <_sp->alphabetSize(); ++state){ // only state=1 is printed
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		AncestralReonstructPosteriorSumStream<<mynode->name()<<"\t"<<state<<"\t"<<probStatesSum[mynode->id()][state]<<"\t";
		if(!mynode->isRoot())
			AncestralReonstructPosteriorSumStream<<mynode->father()->name()<<"\t"<<state<<"\t"<<probStatesSum[mynode->father()->id()][state]<<endl;
		else
			AncestralReonstructPosteriorSumStream<<"NoFather2Root"<<"\t"<<"NA"<<"\t"<<"NA"<<endl;

	}
	//}

	// print Tree
	string TreeAncRecSum = gainLossOptions::_outDir + "//" + "TreeAncRecSum" + ".ph";
	ofstream TreeAncRecSumStr(TreeAncRecSum.c_str());
	TreeAncRecSumStr.precision(LOW_PRECISION);	
	printTreeStatesAsBPValues(TreeAncRecSumStr,probOnesSum,_tr);

	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}




/********************************************************************************************
*********************************************************************************************/
void gainLoss::computeBranchLegthDiffFactor(ostream& out){
	LOGnOUT(4,<<endl<<"Starting computeBranchLegthDiffFactor..."<<endl);
	LOGnOUT(4,<<" Likelihood reference (computed after BBL)="<<_logL<<endl);

	MDOUBLE percentOfLogLDiffTolerance = 0.01;
	MDOUBLE logLOrig;
	MDOUBLE branchLegthOrig;
	MDOUBLE branchLegthAfterBBL;
	tree treeComp = _tr;	// copy the tree
	treeIterTopDownConst tIt(treeComp);
	
	out<<"branch"<<"\t"<<"length@orginal"<<"\t"<<"lengthAfterBBL"<<"\t"<<"factor"<<"\t"<<"Diff"<<"\t"<<"logL@orginal"<<"\t"<<"logLAfterBBL"<<"\t"<<"logL_Diff"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot()) continue;

		branchLegthOrig = _trOrig.findNodeByName(mynode->name())->dis2father();
		branchLegthAfterBBL = _tr.findNodeByName(mynode->name())->dis2father();		
		
		treeComp.findNodeByName(mynode->name())->setDisToFather(branchLegthOrig);	// set BL to original
		
		if(_unObservableData_p){
			if(!gainLossOptions::_gainLossDist){_unObservableData_p->setLforMissingData(treeComp,_sp);}
			else{_unObservableData_p->setLforMissingData(treeComp,_spVVec,_gainDist,_lossDist);}
		}
		if(!gainLossOptions::_gainLossDist){logLOrig = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(treeComp,_scUniqPatterns,*_sp,_weightsUniqPatterns,_unObservableData_p);}
		else{logLOrig = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(treeComp,_scUniqPatterns,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);}
		
		treeComp.findNodeByName(mynode->name())->setDisToFather(branchLegthAfterBBL);	// set BL back

		if(logLOrig > _logL+((percentOfLogLDiffTolerance/100.0)*abs(_logL)) ){
			LOGnOUT(4,<<"WARN... logL with  estimated BL=" <<_logL<<" is lower than original BL="<<logLOrig<<endl);
		}
		LOGnOUT(6,<<"Likelihood for branch Diff="<<branchLegthAfterBBL-branchLegthOrig<<"\t previous logL ="<<logLOrig<<endl);
		out<<mynode->name()<<"\t"
			<<branchLegthOrig<<"\t"
			<<branchLegthAfterBBL<<"\t"
			<<branchLegthAfterBBL / branchLegthOrig <<"\t"
			<<branchLegthAfterBBL - branchLegthOrig<<"\t"
			<<logLOrig <<"\t"<<_logL<<"\t"<<_logL-logLOrig<<endl;
	}
}






/********************************************************************************************
*********************************************************************************************/
void gainLoss::startSimulateSequences(int numOfSequenceSets, int seqLengthInSeq)
{
	LOGnOUT(4,<<endl<< "simulating sequences with the same rate as _rates. numOfSequenceSets="<<numOfSequenceSets<<endl);	
	simulateSequences(numOfSequenceSets, seqLengthInSeq, gainLossOptions::_writeSeqSim,
		gainLossOptions::_useTheSameSpForSim, gainLossOptions::_isReversibleSim, gainLossOptions::_gainEQlossSim,gainLossOptions::_rateDistributionTypeSim);
	
	
	if(!gainLossOptions::_gainLossDist&&gainLossOptions::_calculateRate4siteSim){
		for(int i=0; i<numOfSequenceSets; ++i){
			//re-open seq
			gainLossAlphabet alph;
			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedSequences" + "//" + "seq" + int2string(i+1) + ".fa";
			ifstream in(strSeqNum.c_str());
			sequenceContainer seqReOpened = recognizeFormat::read(in,&alph);

			string outDirSeq = gainLossOptions::_outDir + "//" + "SimulatedSequences" + "//" + "seq" + int2string(i+1);
			createDir(gainLossOptions::_outDir + "//" + "SimulatedSequences", "seq" + int2string(i+1));
			rate4siteGL r4s(seqReOpened,_tr,_sp, outDirSeq,_unObservableData_p);
			r4s.run();
	
		}
	}
}

/********************************************************************************************
simulateSequences
*********************************************************************************************/
vector<sequenceContainer>  gainLoss::simulateSequences(int numOfSequenceSets, int seqLengthInSet, bool writeSeq,
													   bool useTheSame, bool isReversible, bool isGeqL, gainLossOptions::distributionType rateDistributionTypeSim) 
{
	int numOfSitesInSeq= seqLengthInSet;
	LOGnOUT(4,<< "simulating numOfSitesInSeq="<<numOfSitesInSeq<<endl);
	time_t t1,t2;
	time(&t1);

	gainLossAlphabet alph;
	vector<sequenceContainer> scV;
	scV.resize(numOfSequenceSets);

	tree trForSim;
	stochasticProcess* spForSim =NULL;

	if(useTheSame){
		LOGnOUT(4,<< "simulating sequences with the same stochastic proess"<<endl);	
		spForSim = _sp;
		trForSim = _tr;
		printModellValuesOfParams(spForSim, trForSim);		
	}
	else{
		LOGnOUT(4,<< "simulating sequences with the NEW stochastic proess"<<endl);
		LOGnOUT(4,<< "simulating sequences with a Reversible stochastic proess="<<isReversible<<endl);
		LOGnOUT(4,<< "simulating sequences with _gainEQlossSim="<<isGeqL<<endl);
		if(isGeqL){
			LOGnOUT(4,<< "WARNING: _gainLossRateAreFreq is overwritten with"<<isGeqL<<endl);	
			//Parameters::updateParameter("_gainEQloss","1");	// override to previous value assumed to be false
			gainLossOptions::_gainEQloss  =1;	// override to previous value assumed to be false
			//Parameters::updateParameter("_characterFreqEval","FiftyFifty");
			gainLossOptions::_characterFreqEval = gainLossOptions::FiftyFifty;
			//Parameters::updateParameter("_isReversible","1");
			gainLossOptions::_isReversible =1;
		}
		spForSim = startStochasticProcessGeneric(rateDistributionTypeSim, isReversible);
		//gainLossOptimizer glOpt(_tr,spForSim,_sc,
		//	gainLossOptions::_epsilonOptimizationIterationCycle,gainLossOptions::_maxNumOfIterations,
		//	gainLossOptions::_epsilonOptimizationModel,gainLossOptions::_maxNumOfIterationsModel,
		//	gainLossOptions::_epsilonOptimizationBBL,gainLossOptions::_maxNumOfIterationsBBL,_pLforMissingDataPerCat);
		bool isBBL = true; // in optimizer also check gainLossOptions::_isBBL. This if to differ from manyStarts
		gainLossOptimizer glOpt(_tr,spForSim,_scUniqPatterns,
			gainLossOptions::_epsilonOptimizationIterationCycle,gainLossOptions::_maxNumOfIterations,
			gainLossOptions::_epsilonOptimizationModel,gainLossOptions::_maxNumOfIterationsModel,
			gainLossOptions::_epsilonOptimizationBBL,gainLossOptions::_maxNumOfIterationsBBL,_weightsUniqPatterns,_unObservableData_p,
			isBBL ,gainLossOptions::_isbblLSWhenbblEMdontImprove);
		trForSim = glOpt.getOptTree();

	}

	if(_rates.size()==0 && gainLossOptions::_rateDistributionType==gainLossOptions::GAMMA){	// to fill _LpostPerCat - run computeRate4site()		
		rate4siteGL r4s(_sc,_tr,_sp, gainLossOptions::_outDir,_unObservableData_p);
		r4s.run();
		_postProbPerCatPerPos = r4s.getLpostPerCat();
		_rates  = r4s.getRates();
	}
	if(writeSeq)
		createDir(gainLossOptions::_outDir, "SimulatedSequences");
	for(int i=0; i<numOfSequenceSets; ++i){
		LOGnOUT(4,<< "simulating set="<<i<<endl);
		simulateTree st(trForSim, *spForSim, &alph);
		//st.generate_seqWithRateVector(_rates, _sc.seqLen());
		st.generate_seq(numOfSitesInSeq);
		scV[i] = st.toSeqDataWithoutInternalNodes();		
		if(writeSeq){			
			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedSequences" + "//" + "seq" + int2string(i+1) + ".fa";
			ofstream seq_out(strSeqNum.c_str());
			fastaFormat::  write(seq_out,scV[i]);
		}
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return scV;
}

/********************************************************************************************
simulateSequencesForParametricBootstrap
*********************************************************************************************/
sequenceContainer  gainLoss::simulateSequencesForParametricBootstrap(int seqLengthInSet, sequenceContainer& scSampled, tree& trSampled, bool writeSeq,
													   bool useTheSame) 
{
	LOGnOUT(4,<< "simulateSequencesForParametricBootstrap numOfSitesInSeq="<<seqLengthInSet<<endl);

	time_t t1,t2;
	time(&t1);

	MDOUBLE lowRateFactor = 0.1; // not clear it's working...
	MDOUBLE fractionOfPosForLowRate = 0.1;
	int numOfPos2SimulateLowRate = (int)(fractionOfPosForLowRate*seqLengthInSet);
	if(gainLossOptions::_isAddSimulationsWithLowRate){
		seqLengthInSet -= numOfPos2SimulateLowRate; // keep the entire number of positions fixed.
		LOGnOUT(4,<<"Simulate low rate for "<<numOfPos2SimulateLowRate<<" positions (with tree branches multiplied by "<<lowRateFactor<< ")"<<endl);
	}
	gainLossAlphabet alph;
	sequenceContainer sc;
	//if(gainLossOptions::_accountForMissingData)
		//sc.startZeroSequenceContainerGL(scSampled,alph,gainLossOptions::_minNumOfZeros,gainLossOptions::_minNumOfOnes); // reverse from the Zero sequence
		sc.startZeroSequenceContainerGL(scSampled,alph,0,0); // Just as precurasor for next concat.
	//fastaFormat::write(cout,sc); // DEBUG

	tree trForSim;
	tree trForSimForLowRate;
	stochasticProcess* spForSim =NULL;

	trForSim = trSampled;
	if(!gainLossOptions::_gainLossDist){
		spForSim = _sp;
		printModellValuesOfParams(spForSim, trForSim);

		simulateTree st(trForSim, *spForSim, &alph);
		//st.generate_seqWithRateVector(_rates, scSampled.seqLen());
		st.generate_seq(seqLengthInSet);
		sequenceContainer scTemp = st.toSeqDataWithoutInternalNodes();
		/*if(sc.seqLen()>0)*/
			sc.concatenate(scTemp);
		/*else
			sc = scTemp;*/
	}else{
		int numOfSpGain  = _spVVec.size();
		int numOfSpLoss  = _spVVec[0].size();
		int numOfSps = numOfSpGain*numOfSpLoss;
		int numOfPos2SimulatePerSp = seqLengthInSet/numOfSps;
		for (int gainCategor=0; gainCategor<numOfSpGain; gainCategor++){
			for (int lossCategor=0; lossCategor<numOfSpLoss; lossCategor++){
				spForSim = _spVVec[gainCategor][lossCategor];
				simulateTree st(trForSim, *spForSim, &alph);
				//st.generate_seqWithRateVector(_rates, scSampled.seqLen());
				st.generate_seq(numOfPos2SimulatePerSp);
				sequenceContainer scTemp = st.toSeqDataWithoutInternalNodes();
				/*if(sc.seqLen()>0)*/
					sc.concatenate(scTemp);
				/*else
					sc = scTemp;*/
			}
		}
		spForSim = _spSimple;
	}
	if(gainLossOptions::_isAddSimulationsWithLowRate){
		trForSimForLowRate = trSampled;
		trForSimForLowRate.multipleAllBranchesByFactor(lowRateFactor);
		simulateTree stLowRate(trForSimForLowRate, *spForSim, &alph);
		stLowRate.generate_seq(numOfPos2SimulateLowRate); // add 10% low rate simulations
		sequenceContainer scLowRate = stLowRate.toSeqDataWithoutInternalNodes();
		sc.concatenate(scLowRate);
	}
	if(writeSeq){			
		string strSeqNum = gainLossOptions::_outDir + "//" + "simulatedSeq" + ".fa";
		ofstream seq_out(strSeqNum.c_str());
		fastaFormat::  write(seq_out,sc);
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return sc;
}



/********************************************************************************************
Co-Evolution
*********************************************************************************************/
void gainLoss::findCoEvolvingSites(const int numberOfSequences2simulateForCoEvol) {
	// 1. get the observed Vi array	
	if(_postProbPerCatPerPos.size()==0){	// to fill LpostPerCat - run computeRate4site()
		rate4siteGL r4s(_sc,_tr,_sp,gainLossOptions::_outDir, _unObservableData_p);
		r4s.run();
		_rates = r4s.getRates();
		_postProbPerCatPerPos = r4s.getLpostPerCat();
	}	
	computeCountsGL countsGL(_sc,_tr,_sp,gainLossOptions::_outDir,_postProbPerCatPerPos,_distanceFromNearestOTUForRecent);
	countsGL.run();
	VVVVdouble posteriorsGivenTerminals;		// probChangesForBranch[pos][nodeID][x][y]
	resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),posteriorsGivenTerminals);
	posteriorsGivenTerminals = countsGL.getExpChanges();


	// 2. get the simulated Vi array*s*
	LOGnOUT(4,<<endl<< "simulating sequences with the same rate as _rates. numOfSequenceSets="<<numberOfSequences2simulateForCoEvol<<endl);	
	createDir(gainLossOptions::_outDir, "SimulatedSequences");

	simulateSequences(numberOfSequences2simulateForCoEvol,_sc.seqLen(), gainLossOptions::_writeSeqSim,
		gainLossOptions::_useTheSameSpForSim, gainLossOptions::_isReversibleSim, gainLossOptions::_gainEQlossSim,gainLossOptions::_rateDistributionTypeSim);

	VVVVVdouble posteriorsGivenTerminalsSim;		// posteriorsGivenTerminalsSim[Seq][pos][nodeID][x][y]
	posteriorsGivenTerminalsSim.resize(numberOfSequences2simulateForCoEvol);
	bool isSilent = true;
	for(int i=0; i<numberOfSequences2simulateForCoEvol; ++i){
		//re-open seq
		gainLossAlphabet alph;
		string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedSequences" + "//" + "seq" + int2string(i+1) + ".fa";
		ifstream in(strSeqNum.c_str());
		sequenceContainer seqReOpened = recognizeFormat::read(in,&alph);

		string outDirSeq = gainLossOptions::_outDir + "//" + "SimulatedSequences" + "//" + "seq" + int2string(i+1);
		createDir(gainLossOptions::_outDir + "//" + "SimulatedSequences", "seq" + int2string(i+1));		

		computeCountsGL countsGL(seqReOpened,_tr,_sp,outDirSeq,_postProbPerCatPerPos, isSilent);
		countsGL.run();
		//countsGL.printExpectationPerBranch();
		//countsGL.printProbabilityPerPosPerBranch();
		//countsGL.printProbExp();

		resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_sp->alphabetSize(),_sp->alphabetSize(),posteriorsGivenTerminalsSim[i]);
		posteriorsGivenTerminalsSim[i] = countsGL.getExpChanges();
	}

	// 3. Call a general class the finds co-evolving sites based on these VI arrays.
	VVdouble correlations; //[pos][pos]. The correlation between position i and position j.
	correlations.resize(_sc.seqLen());
	for (int k=0; k < correlations.size(); ++k) correlations[k].resize(_sc.seqLen());

	for (int i=0; i < posteriorsGivenTerminals.size() ; ++i) {
		for (int j=i+1; j < posteriorsGivenTerminals.size() ; ++j) {
			correlations[i][j] = computeCorrelationBetweenVis(posteriorsGivenTerminals[i],posteriorsGivenTerminals[j]);
		}
	}

	// computing the correlations between the simulated sequences
	VVVdouble correlationsSim; //[sim][pos][pos]
	resizeVVV(numberOfSequences2simulateForCoEvol,_sc.seqLen(),_sc.seqLen(),correlationsSim);
	for (int k=0; k < correlationsSim.size(); ++k) {
		for (int i=0; i < posteriorsGivenTerminals.size() ; ++i) {
			for (int j=i+1; j < posteriorsGivenTerminals.size() ; ++j) {
				correlationsSim[k][i][j] = computeCorrelationBetweenVis(posteriorsGivenTerminalsSim[k][i],posteriorsGivenTerminalsSim[k][j]);
			}
		}	
	}

	// sort and find where the actual corr is with respect to the simualted sequences.

	//	CoEvol glCoEvo(
	//	LOGnOUT(3,<<" starting to compute co evolving sites... "<<endl);
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLoss::computeCorrelationBetweenVis(const VVVdouble & VIpos_i, const VVVdouble & VIpos_j
											   //, char corrType
											   )
{
// corrType will be 0 for correlations of 0>1 in both (the two positions underand 
// the function gets as input two vectors of substitutions - one for position i and one for position j.
// it then computes the correlation between these two vectors by computing cov (vi, vj)/(sd(vi),sd(vj)).
// VIpos_i has the general structur [nodeId][char][char]
// 1. computing e(x,y)
	MDOUBLE corr = 0.0;
	MDOUBLE EXY = 0.0;
	MDOUBLE EX = 0.0;
	MDOUBLE EY = 0.0;
	for (int i=0; i < VIpos_i.size(); ++i) {// going over all nodes
		MDOUBLE tmp1 = VIpos_i[i][0][1]-VIpos_i[i][1][0];
		MDOUBLE tmp2 = VIpos_j[i][0][1]-VIpos_j[i][1][0];
		EX += tmp1;
		EY += tmp2;
		EXY += (tmp1*tmp2);
	}
	EXY /= VIpos_i.size();
	EX /= VIpos_i.size();
	EY /= VIpos_i.size();
	corr = EXY-EX*EY;
	return corr;
}

/********************************************************************************************
FlatSpBeforeOpt
*********************************************************************************************/
void gainLoss::FlatSpBeforeOpt(stochasticProcess& sp , unObservableData* unObservableData_p){
	LOGnOUT(4,<<"WARNING: FlatSpBeforeOpt.. "<<endl);
	bool isReversible = gainLossOptions::_isReversible;
	bool optimizeAlpha = isAlphaOptimization(sp.distr());
	bool optimizeBeta = isBetaOptimization(sp.distr());
	//bool optimizeMixture = isMixOptimization(sp.distr());
	bool probInvariant = isInvariantOptimization(sp.distr());
	bool evalTheta = isThetaOptimization();

	static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setMu1(1,isReversible);
	if (!isReversible){
		static_cast<gainLossModelNonReversible*>(sp.getPijAccelerator()->getReplacementModel())->setMu2(1);	}
	if(optimizeAlpha){
		setRateAlpha(sp.distr(),0.7);	}
	if(optimizeBeta){
		setRateBeta(sp.distr(),0.7);	}
	if(evalTheta){
		static_cast<gainLossModel*>(sp.getPijAccelerator()->getReplacementModel())->setTheta(0.5);}
	if(probInvariant){
		static_cast<generalGammaDistributionPlusInvariant*>(sp.distr())->setInvProb(0.01);}
	if(gainLossOptions::_isNormalizeQ)
		normalizeQ(&sp);
	if(unObservableData_p)
		unObservableData_p->setLforMissingData(_tr,&sp);
}
/********************************************************************************************
*********************************************************************************************/
void gainLoss::FlatSpBeforeOpt(vector<vector<stochasticProcess*> >& spVVec,distribution * gainDist, distribution * lossDist, unObservableData* unObservableData_p){
	LOGnOUT(4,<<"WARNING: FlatSpBeforeOpt.. "<<endl);
	bool isReversible = gainLossOptions::_isReversible;
	bool optimizeBetaGain = isBetaOptimization(gainDist);
	bool optimizeBetaLoss = isBetaOptimization(lossDist);
	bool optimizeGLProbInvariant = isInvariantOptimization(gainDist); // for both gain and loss
	bool evalTheta = isThetaOptimization();

	stochasticProcess sp = *spVVec[0][0];	
	bool optimizeRateAlpha = isAlphaOptimization((sp.distr()));
	bool optimizeRateProbInvariant = isInvariantOptimization((sp.distr())); 

	updateGainAlpha(1,spVVec,gainDist,lossDist,false);
	if(optimizeBetaGain) updateGainBeta(1,spVVec,gainDist,lossDist,false);
	if(optimizeGLProbInvariant) {
		updateGainProbInvariant(0.01,gainDist);
	}
	// Loss
	if (!isReversible){
		updateLossAlpha(1,spVVec,gainDist,lossDist,false);
		if(optimizeBetaLoss) updateLossBeta(1,spVVec,gainDist,lossDist,false);
		if(optimizeGLProbInvariant) {
			updateGainProbInvariant(0.01,lossDist);
		}
	}
	// overall rate
	if(optimizeRateAlpha)	updateRateAlpha(0.7,spVVec,gainDist,lossDist,false);
	if(optimizeRateProbInvariant)	updateRateProbInvariant(0.01,spVVec,gainDist,lossDist,false);

	if(evalTheta)		updateTheta(0.5,spVVec,gainDist,lossDist,false);
	normalizeQ(spVVec,gainDist,lossDist);
	if(unObservableData_p)
		_unObservableData_p->setLforMissingData(_tr,spVVec,gainDist,lossDist);
}






// prints
/********************************************************************************************
printOptionParameters
*********************************************************************************************/
void gainLoss::printOptionParameters(ostream & out) {
	LOGnOUT(4,<<"\n ---------------------- THE PARAMETERS ----------------------------"<<endl);
	if(gainLossOptions::_gainEQloss)
		LOGnOUT(4,<<"gain=loss model is used. =>freq(0)=freq(1)."<<endl);
	if(Parameters::getInt("_accountForMissingData")){
		LOGnOUT(4,<<"Likelihood computation is performed while acounting for un-oberved data"<<endl);
		LOGnOUT(4,<<"With min number of presences('1's)= "<<Parameters::getInt("_minNumOfOnes")<<endl);
		LOGnOUT(4,<<"With min number of absences('0's)= "<<Parameters::getInt("_minNumOfZeros")<<endl);
	}	
	if(!gainLossOptions::_isReversible && !gainLossOptions::_isRootFreqEQstationary){
		LOGnOUT(4,<<"Fixed-Root ('Non-Rev') model is used"<<endl);}
	else{
		LOGnOUT(4,<<"'Reversible'(Root.freq==stationary.freq) model is used"<<endl);}

	if(gainLossOptions::_isRootFreqEQstationary){
		LOGnOUT(4,<<"RootFreq EQ stationary (taken from each sp - gain/(gain+loss) )"<<endl);
	}
	else
	{
		switch (gainLossOptions::_characterFreqEval){
			case (gainLossOptions::FiftyFifty):  
				LOGnOUT(4,<<"frequencies were set to FiftyFifty "<<endl);
				break;
			case (gainLossOptions::LeavesAve):
				LOGnOUT(4,<<"frequencies are based on LeavesAve (-F option) "<<endl);
				break;
			case (gainLossOptions::optimizeOverTree):
				LOGnOUT(4,<<"frequencies (root '1'(Theta)/'0' freq) are model-based"<<endl);
				break;
		}
	}
	if (gainLossOptions::_treeFile.size()>0) LOGnOUT(4,<<"inTree file: "<< gainLossOptions::_treeFile<<endl);
	LOGnOUT(4,<<"inSeq file: "<<gainLossOptions::_seqFile<<endl);
	if 	(strcmp(gainLossOptions::_referenceSeq.c_str(),"non")!=0) LOGnOUT(4,<<"reference sequence is: "<<gainLossOptions::_referenceSeq<<endl);
	LOGnOUT(4,<<"log: "<<gainLossOptions::_logFile<<" with level= "<<gainLossOptions::_logValue<<endl);
	//LOGnOUT(4,<<"outDir: "<<gainLossOptions::_outDir<<endl);

	// _gainLossDist
	if 	(gainLossOptions::_gainLossDist) {
		if 	(gainLossOptions::_gainLossDistPlusInvariant) {
			LOGnOUT(4,<<"gain , loss ~ GammaPlusInvariant(Alpha,Beta) "<<endl);
			LOGnOUT(4,<<"gain - a Gamma prior distribution with: "<<gainLossOptions::_numberOfGainCategories+1<< " categories (1 Invariant)"<<endl);
			LOGnOUT(4,<<"loss - a Gamma prior distribution with: "<<gainLossOptions::_numberOfLossCategories+1<< " categories (1 Invariant)"<<endl);
		}
		else {
			LOGnOUT(4,<<"gain , loss ~ Gamma(Alpha,Beta) "<<endl);
			LOGnOUT(4,<<"gain - a Gamma prior distribution with: "<<gainLossOptions::_numberOfGainCategories<< " categories"<<endl);
			LOGnOUT(4,<<"loss - a Gamma prior distribution with: "<<gainLossOptions::_numberOfLossCategories<< " categories"<<endl);
		}
	}
	// _performOptimizations
	if(gainLossOptions::_performOptimizations){
		LOGnOUT(4,<< "Optimization of the model parmeters is performed"<<endl);
		if(gainLossOptions::_performOptimizationsManyStarts)
			LOGnOUT(4,<< "performOptimizationsManyStarts with numStarts= "<< gainLossOptions::_numberOfRandPointsInOptimization<<endl);
		switch (gainLossOptions::_rateEstimationMethod){
		case (gainLossOptions::ebExp):  
			{
				if(gainLossOptions::_rateDistributionType == gainLossOptions::GAMMA){
					LOGnOUT(4,<< "rate inference method is: empirical Bayesian estimate"<<endl);
					LOGnOUT(4,<< "using a Gamma prior distribution with: "<<gainLossOptions::_numberOfRateCategories<< " discrete categories"<<endl);
				}			
				else if(gainLossOptions::_rateDistributionType == gainLossOptions::GAMMA_MIXTURE){
					LOGnOUT(4,<< "rate inference method is: empirical Bayesian estimate"<<endl);
					LOGnOUT(4,<< "using a GAMMA_MIXTURE distribution with: "<<gainLossOptions::_numberOfRateComponents<< " components, "<< gainLossOptions::_numberOfRateCategories<< " categories"<<endl);
					if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::EM) LOGnOUT(4,<< "Optimize the Alpha and Beta parameters with EM algorithm"<<endl);
					if(gainLossOptions::_gammmaMixtureOptimizerAlg == gainLossOptions::ONE_DIM) LOGnOUT(4,<< "Optimize the Alpha and Beta parameters with ONE_DIM algorithm"<<endl);
				}break;
			}
		case (gainLossOptions::mlRate): 
			LOGnOUT(4,<< "rate inference method is: maximum likelihood (ML) "<<endl); break;
		}

		if(Parameters::getInt("_performOptimizationsBBL")){
			if(gainLossOptions::_isBblLS){
				LOGnOUT(4,<<"branch lengths optimization is performed using 'Line-Search'"<<endl);}
			else{
				LOGnOUT(4,<<"branch lengths optimization is performed using 'BBL-EM'"<<endl);}
		}
		else{
			LOGnOUT(4,<<"branch lengths are not optimized"<<endl);	}
	}
	else{
		LOGnOUT(4,<< "No optimization is performed"<<endl);
	}


	if(gainLossOptions::_isHGT_normal_Pij){
		//LOGnOUT(4,<<"'Normal' model used with: P01 = gain/(-eigenvalue)-exp(eigenvalue*d)*(1-loss/(-eigenvalue))"<<endl);
	}
	else {
		LOGnOUT(4,<<"The replacement model not allows HGT:  P01 = epsilon*d"<<endl);	}
	if(gainLossOptions::_isHGT_with_Q){
		//LOGnOUT(4,<<"'Normal' replacement model is used with gain => 0"<<endl);
	}
	else {
		LOGnOUT(4,<<"The replacement model with gain = 0"<<endl);	}

	if 	(gainLossOptions::_calculateRate4site) {
		LOGnOUT(4,<<"rate4site is calculated "<<endl);
	}
	if 	(gainLossOptions::_gainLossDist&&gainLossOptions::_calculeGainLoss4site) {
		LOGnOUT(4,<<"gain and loss 4site are calculated "<<endl);
	}
	if(gainLossOptions::_calculePosteriorExpectationOfChange){
		if(gainLossOptions::_isAnaliticComputeJumps){
			LOGnOUT(4,<<"calculePosteriorExpectationOfChange is done Analytically"<<endl);}
		else{
			LOGnOUT(4,<<"calculePosteriorExpectationOfChange is done with "<<gainLossOptions::_numOfSimulationsForPotExp<<" simulations"<<endl);}
	}
	LOGnOUT(4,<<"-------------------------------------------------------------------"<<endl);
}
/********************************************************************************************
printPij_t
*********************************************************************************************/
void gainLoss::printPij_t(MDOUBLE dist, ostream& out){
	out<<"-------------------------------"<<endl<<"The Pij("<<dist<<"): Matrix:"<<endl;
	if(gainLossOptions::_gainLossDist){
		MDOUBLE spPij_t00=0;
		MDOUBLE spPij_t01=0;
		MDOUBLE spPij_t10=0;
		MDOUBLE spPij_t11=0;		
		int numOfSPs = _gainDist->categories()*_lossDist->categories();
		for (int i=0; i < numOfSPs; ++i) {
			int gainIndex =fromIndex2gainIndex(i,_gainDist->categories(),_lossDist->categories());
			int lossIndex =fromIndex2lossIndex(i,_gainDist->categories(),_lossDist->categories());			
			spPij_t00 += _spVVec[gainIndex][lossIndex]->Pij_t(0,0,dist)* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t01 += _spVVec[gainIndex][lossIndex]->Pij_t(0,1,dist)* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t10 += _spVVec[gainIndex][lossIndex]->Pij_t(1,0,dist)* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t11 += _spVVec[gainIndex][lossIndex]->Pij_t(1,1,dist)* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
		}
		out<<"p0,0["<<dist<<"]: "<<spPij_t00<<endl;
		out<<"p0,1["<<dist<<"]: "<<spPij_t01<<endl;
		out<<"p1,0["<<dist<<"]: "<<spPij_t10<<endl;
		out<<"p1,1["<<dist<<"]: "<<spPij_t11<<endl;
		out<<endl;
	}
	else{
		out<<"p0,0["<<dist<<"]: "<<_sp->Pij_t(0,0,dist)<<endl;
		out<<"p0,1["<<dist<<"]: "<<_sp->Pij_t(0,1,dist)<<endl;
		out<<"p1,0["<<dist<<"]: "<<_sp->Pij_t(1,0,dist)<<endl;
		out<<"p1,1["<<dist<<"]: "<<_sp->Pij_t(1,1,dist)<<endl;
		out<<endl;
	}
}
/********************************************************************************************
printQ
*********************************************************************************************/
void gainLoss::printQ(ostream& out){
	VVdouble Q;
	out<<"-------------------------------"<<endl<<"The Q Matrix:"<<endl;
	if(gainLossOptions::_gainLossDist){
		MDOUBLE spPij_t00=0;
		MDOUBLE spPij_t01=0;
		MDOUBLE spPij_t10=0;
		MDOUBLE spPij_t11=0;		
		int numOfSPs = _gainDist->categories()*_lossDist->categories();
		for (int i=0; i < numOfSPs; ++i) {
			int gainIndex =fromIndex2gainIndex(i,_gainDist->categories(),_lossDist->categories());
			int lossIndex =fromIndex2lossIndex(i,_gainDist->categories(),_lossDist->categories());			
			Q = (static_cast<gainLossModel*>(_spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel())->getQ());
			spPij_t00 += Q[0][0]* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t01 += Q[0][1]* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t10 += Q[1][0]* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
			spPij_t11 += Q[1][1]* _gainDist->ratesProb(gainIndex)*_lossDist->ratesProb(lossIndex);
		}
		out<<"Q[0][0]= "<<spPij_t00<<endl;
		out<<"Q[0][1]= "<<spPij_t01<<endl;
		out<<"Q[1][0]= "<<spPij_t10<<endl;
		out<<"Q[1][1]= "<<spPij_t11<<endl;
		out<<endl;
	}
	else{
		VVdouble Q = (static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getQ());
		//out<<"freq[0]*Q[0][1]= "<<(1-static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getTheta())*Q[0][1]<<endl;	
		//out<<"freq[1]*Q[1][0]= "<<(static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getTheta())*Q[1][0]<<endl;
		out<<"Q[0][0]= "<<Q[0][0] <<endl;
		out<<"Q[0][1]= "<<Q[0][1] <<endl;
		out<<"Q[1][0]= "<<Q[1][0] <<endl;
		out<<"Q[1][1]= "<<Q[1][1] <<endl;
		out<<endl;
	}
}
/********************************************************************************************
printTreeLikelihoodAllPosAlphTheSame
*********************************************************************************************/
void gainLoss::printTreeLikelihoodAllPosAlphTheSame(bool isLOGnOUT ,ostream& out)
{
	MDOUBLE res;
	if(!gainLossOptions::_gainLossDist){
		res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,*_sp,_weightsUniqPatterns,_unObservableData_p);
		out.precision(9);
		if(isLOGnOUT){
			LOGnOUT(3,<<"The Tree Likelihood AllPosAlphTheSame is "<<res<<endl);}
		else{
			out<<"The Tree Likelihood AllPosAlphTheSame is "<<res<<endl;}
	}
	else{
		res = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,_spVVec,_gainDist,_lossDist,_weightsUniqPatterns,_unObservableData_p);
		out.precision(9);
		if(isLOGnOUT){
			LOGnOUT(3,<<"The Tree Likelihood AllPosAlphTheSame is "<<res<<endl);}
		else{
			out<<"The Tree Likelihood AllPosAlphTheSame is "<<res<<endl;}
		//res = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSameNoComputeUp(_tr,_sc,_spVVec,_gainDist,_lossDist);
		//out<<"The Tree Likelihood AllPosAlphTheSameNoComputeUp is "<<res<<endl;
	}
	_logL = res; // update the tree likelihood.
}


/********************************************************************************************
printLofPosBothModels
*********************************************************************************************/
void gainLoss::printLofPosBothModels(){
	LOGnOUT(4,<<"Starting printLofPosBothModels..."<<endl);
	string LofPosBothModels = gainLossOptions::_outDir + "//" + "printLofPosBothModels.txt"; 
	ofstream likeOfPosBothModelsStream(LofPosBothModels.c_str());
	likeOfPosBothModelsStream.precision(PRECISION);

	MDOUBLE treeL = printLofPosBothModels(likeOfPosBothModelsStream);
	LOGnOUT(4,<<"treeL= "<<treeL<<endl);
}

/********************************************************************************************
printLofPosBothModels
*********************************************************************************************/
MDOUBLE gainLoss::printLofPosBothModels(ostream& out){
	
	//MDOUBLE mu1 = static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getMu1();
	MDOUBLE res =0;
	unObservableData* unObservableData_p_0;
// single stochastic process	
	if(!gainLossOptions::_gainLossDist){
		out<<"The likelihood of each pos when gain and loss are scalars"<<endl;

		computePijGam piModel_0, piModel_1;		
		piModel_1.fillPij(_tr,*_sp);

		stochasticProcess* spModel_0 = _sp->clone();
		//static_cast<gainLossModel*>((*_sp).getPijAccelerator()->getReplacementModel())->setMu1(0.0,gainLossOptions::_isReversible);
		static_cast<gainLossModel*>((*spModel_0).getPijAccelerator()->getReplacementModel())->setMu1(0.0,gainLossOptions::_isReversible);	//NO NEED to update since the _sp is byRef

		piModel_0.fillPij(_tr,*spModel_0);
		if(_unObservableData_p){
			unObservableData_p_0 = new unObservableData(_sc, spModel_0, gainLossAlphabet(),Parameters::getInt("_minNumOfOnes"), Parameters::getInt("_minNumOfZeros"));
			unObservableData_p_0->setLforMissingData(_tr,spModel_0);
		}
		else
			unObservableData_p_0 = NULL;


		MDOUBLE LnofPos_Model_0, LnofPos_Model_1;
		int k;
		out<<"POS"<<"\t"<<"M_gain0"<<"\t"<<"M"<<"\t"<<"Diff"<<endl;
		for (k=0; k < _scWithFullLength.seqLen(); ++k) {
			LnofPos_Model_0 = log(likelihoodComputation::getLofPos(k,_tr,_scWithFullLength,piModel_0,*spModel_0));
			LnofPos_Model_1 = log(likelihoodComputation::getLofPos(k,_tr,_scWithFullLength,piModel_1,*_sp));
			if(_unObservableData_p){
				LnofPos_Model_1 = LnofPos_Model_1 - log(1- exp(_unObservableData_p->getlogLforMissingData()));
				LnofPos_Model_0 = LnofPos_Model_0 - log(1- exp(unObservableData_p_0->getlogLforMissingData()));
			}
			res += LnofPos_Model_0;
			out<<k+1<<"\t"<<LnofPos_Model_0<<"\t"<<LnofPos_Model_1<<"\t"<<LnofPos_Model_1-LnofPos_Model_0<<endl;
		}
		if(unObservableData_p_0) delete unObservableData_p_0;
		return res;
	}
// multiple stochastic processes	
	else{
		out<<"The likelihood of each pos when gain and loss ~Gamma(Alpha,Beta)"<<endl;
		out<<"...Not implemented (yet)"<<endl;
		return res;
	}
}

/********************************************************************************************
printLofPos
*********************************************************************************************/
void gainLoss::printLofPos(){
	LOGnOUT(4,<<"Starting printLofPos..."<<endl);
	//ofstream likeOfPosStream(gainLossOptions::_outFileLikeofPos.c_str());
	string g4s = gainLossOptions::_outDir + "//" + "likeOfPos.txt";
	ofstream likeOfPosStream(g4s.c_str());
	likeOfPosStream.precision(PRECISION);
	MDOUBLE treeL = printLofPos(likeOfPosStream);
	cout<<"Tree logL="<<treeL<<"\n";
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLoss::printLofPos(ostream& out){

	MDOUBLE res =0;
	out<<"# log Likelihood of tree (entire data)="<<"\t"<<_logL<<endl;
	out<<"POS"<<"\t"<<"logLofPos"<<endl;
	if(!gainLossOptions::_gainLossDist){
		//out<<"The likelihood of each pos when gain and loss are scalars"<<endl;
		computePijGam pi;
		pi.fillPij(_tr,*_sp);
		MDOUBLE LnofPos;
		for (int k=0; k < _scWithFullLength.seqLen(); ++k) {
			LnofPos = log(likelihoodComputation::getLofPos(k,_tr,_scWithFullLength,pi,*_sp));
			if(_unObservableData_p)
				LnofPos = LnofPos - log(1- exp(_unObservableData_p->getlogLforMissingData()));			
			res += LnofPos;
			out<<k+1<<"\t"<<LnofPos<<endl;
			// DEB
			if(LnofPos>0)
				out<<k+1<<"\t"<<LnofPos<<endl;
		}
		return res;
	}
	else{
		int numOfRateCategories = _spVVec[0][0]->categories();
		vector<computePijGam> pi_vec(numOfRateCategories);
		vector<suffStatGlobalGam> ssc_vec(numOfRateCategories);
		vector<computeUpAlg> cup_vec(numOfRateCategories);
		likelihoodComputationGL::fillPijAndUp(_tr,_sc, _spVVec,_gainDist,_lossDist,pi_vec,ssc_vec,cup_vec);
		Vdouble posLike;
		res = likelihoodComputationGL::getTreeLikelihoodFromUp2(_tr,_sc,_spVVec,ssc_vec,_gainDist, _lossDist,NULL,_unObservableData_p,&posLike);
		for (int k=0; k < _sc.seqLen(); ++k) {			
			out<<k+1<<"\t"<<posLike[k]<<endl;
		}
		return res;
	}
}
/********************************************************************************************
Util - printLikelihoodLandscape
*********************************************************************************************/
void gainLoss::printLikelihoodLandscape(stochasticProcess* sp){
	LOGnOUT(4,<<"start printLikelihoodLandscape for: ..."<<endl);
	if(gainLossOptions::_printLikelihoodLandscapeAlphaRate)
		LOGnOUT(4,<<" AlphaRate"<<endl);
	if(gainLossOptions::_printLikelihoodLandscapeGainLoss)
		LOGnOUT(4,<<"Gain and Loss"<<endl);
	if(gainLossOptions::_printLikelihoodLandscapeTheta)
		LOGnOUT(4,<<"Theta"<<endl);

	stochasticProcess* spTemp = sp->clone();
	string LikelihoodLandscape = gainLossOptions::_outDir + "//" + "LikelihoodLandscape.txt"; 
	ofstream LikelihoodLandscapeStream(LikelihoodLandscape.c_str());
	LikelihoodLandscapeStream.precision(PRECISION);
	LikelihoodLandscapeStream<<"Alpha"<<"\t"<<"Gain"<<"\t"<<"Loss"<<"\t"<<"Theta"<<"\t"<<"L"<<endl;
	cout<<"Alpha"<<"\t"<<"Gain"<<"\t"<<"Loss"<<"\t"<<"Theta"<<"\t"<<"L"<<endl;

	bool optimizeAlpha = isAlphaOptimization(_sp->distr());
	//bool optimizeBeta = isBetaOptimization(_sp->distr());
	//bool optimizeMixture = isMixOptimization(_sp->distr());
	//bool probInvariant = isInvariantOptimization(_sp->distr());
	bool evalTheta = isThetaOptimization();

	MDOUBLE AlphaRate,Gain,Loss,Theta;
	MDOUBLE Increment = 0.01;
	int BigEnoughToEndLoop = 100000000;
	MDOUBLE LL;

	// get all original values
	if(optimizeAlpha) 
		AlphaRate =static_cast<gammaDistribution*>(spTemp->distr())->getAlpha();
	Gain = static_cast<gainLossModel*>((*spTemp).getPijAccelerator()->getReplacementModel())->getMu1();
	if(!gainLossOptions::_isReversible)
		Loss= static_cast<gainLossModelNonReversible*>((*spTemp).getPijAccelerator()->getReplacementModel())->getMu2();
	if(evalTheta)
		Theta= static_cast<gainLossModel*>((*spTemp).getPijAccelerator()->getReplacementModel())->getTheta();

	// start the 1-3way loop for landscape	
	for (int i=1; i*Increment<=gainLossOptions::_userAlphaRateMax; i++){
		if(gainLossOptions::_printLikelihoodLandscapeAlphaRate){
			AlphaRate = i*Increment;
			if(optimizeAlpha)	setRateAlpha(spTemp->distr(),AlphaRate);
		}
		else
			i=BigEnoughToEndLoop;
		for (int j=1; j*Increment<=gainLossOptions::_userGainMax; j++){
			if(gainLossOptions::_printLikelihoodLandscapeGainLoss){
				Gain = j*Increment;
				static_cast<gainLossModel*>(spTemp->getPijAccelerator()->getReplacementModel())->setMu1(Gain, gainLossOptions::_isReversible);
			}
			else
				j=BigEnoughToEndLoop;			
			for (int k=1; k*Increment<=gainLossOptions::_userLossMax; k++){
				if(gainLossOptions::_printLikelihoodLandscapeGainLoss){
					Loss = k*Increment;
					if (!gainLossOptions::_isReversible) static_cast<gainLossModelNonReversible*>(spTemp->getPijAccelerator()->getReplacementModel())->setMu2(Loss);
				}
				else
					k=BigEnoughToEndLoop;
				for (int l=1; l*Increment<=gainLossOptions::_userThetaMax; l++){
					if(gainLossOptions::_printLikelihoodLandscapeTheta){
						Theta = l*Increment;
						if(evalTheta) static_cast<gainLossModel*>(spTemp->getPijAccelerator()->getReplacementModel())->setTheta(Theta);
					}
					else
						l=BigEnoughToEndLoop;							
					LL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*spTemp,_weightsUniqPatterns,_unObservableData_p);
					LikelihoodLandscapeStream<<AlphaRate<<"\t"<<Gain<<"\t"<<Loss<<"\t"<<Theta<<"\t"<<LL<<endl;
					cout<<AlphaRate<<"\t"<<Gain<<"\t"<<Loss<<"\t"<<Theta<<"\t"<<LL<<endl;
				}
			}
		}
	}
	delete spTemp;
}
/********************************************************************************************
printLikelihoodLandscapeStatFreqRatioAndRootFreqRatio
*********************************************************************************************/
void gainLoss::printLikelihoodLandscapeStatFreqRatioAndRootFreqRatio(){
	bool isCloneEveryIteration = false;	// otherwise same model in all iteration (other than the update)
	stochasticProcess* spTemp=NULL;
	vector<vector<stochasticProcess*> > spVVec;
	unObservableData* unObservableData_p=NULL;
	distribution* gainDist=NULL;
	distribution* lossDist=NULL;
	tree tempTree = _tr;
	if(_unObservableData_p)
		unObservableData_p = _unObservableData_p->clone();

	
	if(gainLossOptions::_gainLossDist){
		LOGnOUT(4,<<"start printLikelihoodLandscape for: gainLossDist (mixture) with gainLoss ratio and Theta (Root'1'Freq)"<<endl);
		LOGnOUT(4,<<"increment="<<gainLossOptions::_likelihoodLandscapeIncrement<<endl);
		//LOGnOUT(4,<<"WARNING: the _spVVec,_gainDist,_lossDist are overwritten"<<endl);
		cloneSpVVec(_spVVec,spVVec);
		gainDist = _gainDist->clone();
		lossDist = _lossDist->clone();		
		if(gainLossOptions::_optBBL_LS_InIteration || gainLossOptions::_optBBL_EM_InIteration){
			errorMsg::reportError("Error: BBL not implemented with gainLossDist for printLikelihoodLandscape\n");
		}
	}
	else{
		if(!gainLossOptions::_gainLossRateAreFreq){
			LOGnOUT(4,<<"WARNING:: choose _gainLossRateAreFreq for printLikelihoodLandscapeStatFreqRatioAndRootFreqRatio\n");
		}
		LOGnOUT(4,<<"start printLikelihoodLandscape for: Gain (Stationary'1'Freq) and Theta (Root'1'Freq)"<<endl);
		LOGnOUT(4,<<"increment="<<gainLossOptions::_likelihoodLandscapeIncrement<<endl);
		LOGnOUT(4,<<"gainLossOptions::_optAlphaInIteration="<<gainLossOptions::_optAlphaInIteration<<" optBBL_LS_InIteration="<<gainLossOptions::_optBBL_LS_InIteration<<" optBBL_EM_InIteration="<<gainLossOptions::_optBBL_EM_InIteration<<endl);
		spTemp = _sp->clone();
	}

	bool optimizeAlpha=false;
	if(!gainLossOptions::_gainLossDist)
		optimizeAlpha= isAlphaOptimization((*spTemp).distr());
	string LikelihoodLandscape = gainLossOptions::_outDir + "//" + "LikelihoodLandscape.txt"; 
	ofstream LikelihoodLandscapeStream(LikelihoodLandscape.c_str());
	LikelihoodLandscapeStream<<"StationaryFreq"<<"\t"<<"Theta"<<"\t"<<"Like"<<endl;
	LOGnOUT(4,<<"StationaryFreq"<<"\t"<<"Theta"<<"\t"<<"Like"<<endl);
	MDOUBLE AlphaGain,AlphaLoss, BetaGain,BetaLoss,   gainLossRatioToCompleteByBeta, ratio  =1;
	if(gainLossOptions::_gainLossDist){
		AlphaGain = getRateAlpha(_gainDist);
		AlphaLoss = getRateAlpha(_lossDist);
	}

	MDOUBLE Gain,Theta=1;
	MDOUBLE Increment = gainLossOptions::_likelihoodLandscapeIncrement;
	MDOUBLE LL=1;

	MDOUBLE optLike=1;
	MDOUBLE AlphaRate=1;
	MDOUBLE currAlpha=1;

	MDOUBLE tollForPairwiseDist=0.01; // the BBL default, epsilon per branch (brent's value)
	MDOUBLE bblEMfactor = 10;
	int numberOfBranchs = _tr.getNodesNum();
	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/5; // (is 1.5) for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationBBLIter = gainLossOptions::_epsilonOptimizationBBL*epsilonOptimizationIterFactor/bblEMfactor;	// The next iteration epsilon, multiply per-branch value

	// get all original values
	if(gainLossOptions::_optAlphaInIteration){
		if(optimizeAlpha) 
			AlphaRate =static_cast<gammaDistribution*>(spTemp->distr())->getAlpha();
	}


	//////////////////////////////////////////////////////////////////////////
	for (int j=1; j*Increment<=0.99999; j++){
		Gain = j*Increment;
        if(gainLossOptions::_gainLossDist){
			ratio = Gain/(1-Gain);
			gainLossRatioToCompleteByBeta = ratio*(AlphaLoss/AlphaGain);
			BetaGain =sqrt(1/gainLossRatioToCompleteByBeta);			// AlphaGain = 0.35
			BetaLoss =sqrt(gainLossRatioToCompleteByBeta);				// AlphaLoss = 0.9
			updateGainBeta(BetaGain,spVVec,_gainDist,_lossDist);
			updateLossBeta(BetaLoss,spVVec,_gainDist,_lossDist);
		}
		else{
			if(gainLossOptions::_gainLossRateAreFreq)
				static_cast<gainLossModel*>(spTemp->getPijAccelerator()->getReplacementModel())->setMu1(Gain, gainLossOptions::_isReversible);
			else{
				static_cast<gainLossModel*>(spTemp->getPijAccelerator()->getReplacementModel())->setMu1(Gain, gainLossOptions::_isReversible);
				static_cast<gainLossModelNonReversible*>(spTemp->getPijAccelerator()->getReplacementModel())->setMu2((1-Gain));
			}
		}		
		
		//////////////////////////////////////////////////////////////////////////
		for (int l=1; l*Increment<=0.99999; l++){
			tree tempTree = _tr;
			Theta = l*Increment;
			if(gainLossOptions::_gainLossDist){
				updateTheta(Theta,spVVec,_gainDist,_lossDist);
				if(unObservableData_p) unObservableData_p->setLforMissingData(_tr,spVVec,_gainDist,_lossDist); // No need?
			}
			else{
				static_cast<gainLossModel*>(spTemp->getPijAccelerator()->getReplacementModel())->setTheta(Theta);
				if(unObservableData_p)	unObservableData_p->setLforMissingData(_tr,_sp);
			}
			
			if(gainLossOptions::_gainLossDist){
				LL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(_tr,_scUniqPatterns,spVVec,_gainDist,_lossDist,_weightsUniqPatterns,unObservableData_p);
			}
			else{
				LL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*spTemp,_weightsUniqPatterns,unObservableData_p);
			}
// BBL-LS
			if(gainLossOptions::_optBBL_LS_InIteration){
				bblLS bbl;
				LL = bbl.optimizeBranches(tempTree,spTemp,_sc,_weightsUniqPatterns,unObservableData_p,1, gainLossOptions::_epsilonOptimizationBBL, gainLossOptions::_maxNumOfIterationsBBL,LL);
			}
// BBL-EM
			if(gainLossOptions::_optBBL_EM_InIteration){
				bblEM bblEM1(tempTree, _sc, *spTemp, NULL, (int)(gainLossOptions::_maxNumOfIterationsBBL*bblEMfactor), epsilonOptimizationBBLIter,tollForPairwiseDist,unObservableData_p,&LL);
				LL = bblEM1.getTreeLikelihood();
			}
// optAlpha
			if(optimizeAlpha && gainLossOptions::_optAlphaInIteration){
				LL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_tr,_sc,*spTemp,_weightsUniqPatterns,unObservableData_p);
				optLike = -brent(MINIMUM_ALPHA_PARAM,AlphaRate,MAXIMUM_ALPHA_PARAM
					,C_evalParam(_tr,*spTemp,_sc,C_evalParam::rateAlpha,gainLossOptions::_isReversible,_weightsUniqPatterns,unObservableData_p),gainLossOptions::_epsilonOptimizationModel,&currAlpha);			
				if (optLike>LL)
					setRateAlpha(spTemp->distr(),currAlpha);
			}
			LikelihoodLandscapeStream<<Gain<<"\t"<<Theta<<"\t"<<LL<<endl;
			LOGnOUT(4,<<Gain<<"\t"<<Theta<<"\t"<<LL<<endl);
			// clone every iteration (is needed?)
			tempTree = _tr;
			if(isCloneEveryIteration){
				if(gainLossOptions::_gainLossDist){
					deleteSpVVec(&spVVec);
					cloneSpVVec(_spVVec,spVVec);
					delete(gainDist);
					gainDist = _gainDist->clone();
					delete(lossDist);
					lossDist = _lossDist->clone();
					if(unObservableData_p){
						delete unObservableData_p;
						unObservableData_p = _unObservableData_p->clone();
					}

				}
				else{
					delete(spTemp);
					spTemp = _sp->clone();
					if(unObservableData_p){
						delete unObservableData_p;
						unObservableData_p = _unObservableData_p->clone();
					}
				}
			}
		}			
	}
	// final deletions
	if(spTemp)
		delete(spTemp);
	if(spVVec.size()>0)
		deleteSpVVec(&spVVec);
	if(gainDist)
		delete gainDist;
	if(lossDist)
		delete lossDist;
	if(unObservableData_p)
		delete unObservableData_p;

}








/********************************************************************************************
*********************************************************************************************/
void gainLoss::initMixtureParams(Vdouble& initAlphaRates, Vdouble& initBetaRates, Vdouble& initCompProbRates, int numOfGammaComp,
								 MDOUBLE initAlphaRate, MDOUBLE initBetaRate, MDOUBLE initCompProbRate) 
{
	initAlphaRates.resize(numOfGammaComp);
	initBetaRates.resize(numOfGammaComp);
	initCompProbRates.resize(numOfGammaComp);
	for (int i = 0; i < numOfGammaComp; ++i)
	{
		initAlphaRates[i] = initAlphaRate*(numOfGammaComp-i)/numOfGammaComp;
		initBetaRates[i] = initBetaRate*(i+1)/numOfGammaComp;
		initCompProbRates[i] = initCompProbRate/numOfGammaComp;
	}
}


/********************************************************************************************
*********************************************************************************************/
void gainLoss::convertGainLossRatesToFreq(){
	LOGnOUT(4,<<"Starting convertGainLossRatesToFreq..."<<endl);
	MDOUBLE gainLossSum = 0.0;
	if(!gainLossOptions::_gainLossDist){
		MDOUBLE gain = static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->getMu1();
		MDOUBLE loss = static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->getMu2();
		gainLossSum	= gain+loss;		
		static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->setMu1(gain/gainLossSum,gainLossOptions::_isReversible);
		static_cast<gainLossModelNonReversible*>(_sp->getPijAccelerator()->getReplacementModel())->setMu2(loss/gainLossSum);
	}
	else{
		//gainLossSum = normalizeQ(_spVVec, _gainDist, _lossDist);
	}
	_tr.multipleAllBranchesByFactor(gainLossSum); //Needed in order to maintain the overall expected number of event
	printTreeLikelihoodAllPosAlphTheSame();
}


/********************************************************************************************
Normalize the rates by setting the expected number of substitutions per site (per unit time) to 1:
setting Sum over i q_ii*freq_i = 1
*********************************************************************************************/
void gainLoss::normalizeQandTree(bool isComputeLikelihood, bool isMultipleAllBranchesByNormFactor){
	LOGnOUT(4,<<"Starting normalizeQandTree...(so that sumQii=1 (or weighted ave. of sunOii's for many Qs))"<<endl);
	MDOUBLE norm_factor = 0.0;
	if(!gainLossOptions::_gainLossDist){
		norm_factor = normalizeQ(_sp);
	}
	else{
		norm_factor = normalizeQ(_spVVec, _gainDist, _lossDist);
	}
	LOGnOUT(4,<<"Q were multiplied by "<<1.0/norm_factor<<endl);
	if(isMultipleAllBranchesByNormFactor){
		_tr.multipleAllBranchesByFactor(norm_factor); ////Needed in order to maintain the overall expected number of event, Q was multi in 1/norm_factor	
		LOGnOUT(4,<<"Tree branches multi by "<<norm_factor<<endl);
	}
	if(isComputeLikelihood)
		printTreeLikelihoodAllPosAlphTheSame();
}


/********************************************************************************************
This manipulation produces an un normalized Q matrices 
*********************************************************************************************/
void gainLoss::AlphaEqBetaManipulation(){
	LOGnOUT(4,<<"Starting AlphaEqBetaManipulation..."<<endl);
	MDOUBLE lossAlpha = getRateAlpha(_lossDist);
	MDOUBLE lossBeta = getRateBeta(_lossDist);
	MDOUBLE factor2MultiplyBy = lossAlpha/lossBeta;
	bool isNormalizeQ = false;

	updateLossBeta(getRateBeta(_lossDist)*factor2MultiplyBy,_spVVec,_gainDist,_lossDist,isNormalizeQ);
	updateGainBeta(getRateBeta(_gainDist)*factor2MultiplyBy,_spVVec,_gainDist,_lossDist,isNormalizeQ);	
	_tr.multipleAllBranchesByFactor(factor2MultiplyBy);
	if(_unObservableData_p) _unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist); // No need?
	LOGnOUT(4,<<"Finish AlphaEqBetaManipulation.");
	printTreeLikelihoodAllPosAlphTheSame();
}

/********************************************************************************************
Aiming to classify branch specific event as either Recent or Ancient,
compute the distance from root cut-off
This is a basic method to compute the cut-off while finding a balance in total branch lengths
so to minimize "totalBranchLengthAncient - totalBranchLengthRecent"
This method don't consider "distance to OTU" - i.e., that some nodes will be recent by 
*********************************************************************************************/
MDOUBLE gainLoss::computeDistanceFromRootForRecent(tree& tr)
{
	MDOUBLE distanceFromRootForRecentCutOff;
	MDOUBLE MeanDistanceFromRoot;
	//MDOUBLE MeanDistanceFromNearestOTU;

	MDOUBLE totalBranchLengthRecent = 0;
	MDOUBLE totalBranchLengthAncient = 0;
	MDOUBLE diffTotalBranchLengthRecentAncient = 0;
	
	int numberOfNodes = tr.getNodesNum();
	Vdouble DistanceFromRoot(numberOfNodes-1);	// -1 because of Root
	//Vdouble DistanceFromNearestOTU(numberOfNodes);
	Vdouble Distance2father(numberOfNodes-1);

	treeIterDownTopConst tIt(tr);
	int i = 0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())
			break;
		//DistanceFromRoot.push_back(getDistance2ROOT(mynode));
		//Distance2father.push_back(mynode->dis2father());		
		Distance2father[i] = mynode->dis2father();
		DistanceFromRoot[i] = mynode->getDistance2ROOT();
		//DistanceFromNearestOTU[i] =  getMinimalDistance2OTU(mynode);
		//cout<<mynode->name()<<" "<<DistanceFromRoot[i]<<endl;	// DEBUG
		++i;
	}
	MeanDistanceFromRoot = computeAverage(DistanceFromRoot, &Distance2father);
	
	distanceFromRootForRecentCutOff = MeanDistanceFromRoot; // Starting point as Mean
	for(i = 0; i<numberOfNodes; ++i){
		if(DistanceFromRoot[i] < distanceFromRootForRecentCutOff)
			totalBranchLengthAncient+= Distance2father[i];
		else
			totalBranchLengthRecent+= Distance2father[i];
	}	
	diffTotalBranchLengthRecentAncient = totalBranchLengthAncient - totalBranchLengthRecent;
	bool isRecentBiggerAncient = true;
	if(totalBranchLengthAncient>totalBranchLengthRecent)
		isRecentBiggerAncient = false;

	bool isImprovedRecentEstimation = true;
	int numberOfIterations = 0;
	while(isImprovedRecentEstimation && (numberOfIterations<10000)){
		MDOUBLE prevDiffTotalBranchLengthRecentAncient = diffTotalBranchLengthRecentAncient; // init
		MDOUBLE prevDistanceFromRootForRecent = distanceFromRootForRecentCutOff;
		MDOUBLE prevtotalBranchLengthAncient = totalBranchLengthAncient;
		MDOUBLE prevtotalBranchLengthRecent = totalBranchLengthRecent;

		distanceFromRootForRecentCutOff = distanceFromRootForRecentCutOff - diffTotalBranchLengthRecentAncient/(numberOfNodes*100); // cont. correction

		for(i=0, totalBranchLengthAncient=0, totalBranchLengthRecent=0; i<numberOfNodes; ++i){
			if(DistanceFromRoot[i] < distanceFromRootForRecentCutOff)
				totalBranchLengthAncient+= Distance2father[i];
			else
				totalBranchLengthRecent+= Distance2father[i];
		}	
		diffTotalBranchLengthRecentAncient = totalBranchLengthAncient - totalBranchLengthRecent;
		if(abs(diffTotalBranchLengthRecentAncient) > abs(prevDiffTotalBranchLengthRecentAncient) 
			//&& 	((totalBranchLengthAncient>totalBranchLengthRecent)*isRecentBiggerAncient)	// to make sure that Ancient is not more than Recent, wait for "flip"
			)
		{
			isImprovedRecentEstimation = false;
			distanceFromRootForRecentCutOff = prevDistanceFromRootForRecent; // go back to last estimation.
			totalBranchLengthAncient = prevtotalBranchLengthAncient;
			totalBranchLengthRecent = prevtotalBranchLengthRecent;
		}
		//cout<<diffTotalBranchLengthRecentAncient<<" "<<distanceFromRootForRecentCutOff<<"\n";	// DEBUG
		numberOfIterations++;
	}
	LOGnOUT(4,<<"The computed distanceFromRootForRecentCutOff="<<distanceFromRootForRecentCutOff<<" with TotalBranchLength Ancient="<<totalBranchLengthAncient<<" Recent="<<totalBranchLengthRecent<<" Converged "<<!isImprovedRecentEstimation<<endl);
	return distanceFromRootForRecentCutOff;
}



/********************************************************************************************
Aiming to classify branch specific event as either Recent or Ancient,
compute the distance from Leaf
This is a basic method to compute the cut-off while finding a balance in total branch lengths
so to minimize "totalBranchLengthAncient - totalBranchLengthRecent"
*********************************************************************************************/
MDOUBLE gainLoss::computeDistanceNearestOTUforRecent(tree& tr)
{
	MDOUBLE distance2NearestOTUForRecent;
	MDOUBLE MeanDistanceFromNearestOTU;

	MDOUBLE totalBranchLengthRecent = 0;
	MDOUBLE totalBranchLengthAncient = 0;
	MDOUBLE diffTotalBranchLengthRecentAncient = 0;

	int numberOfNodes = tr.getNodesNum();
	Vdouble DistanceFromNearestOTU(numberOfNodes-1);
	Vdouble Distance2father(numberOfNodes-1);

	treeIterDownTopConst tIt(tr);
	int i = 0;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if(mynode->isRoot())
			break;
		Distance2father[i] = mynode->dis2father();
		DistanceFromNearestOTU[i] =  mynode->getMinimalDistance2OTU();
		//cout<<mynode->name()<<" "<<DistanceFromNearestOTU[i]<<"\n";	// DEBUG
		++i;
	}

	MeanDistanceFromNearestOTU = computeAverage(DistanceFromNearestOTU, &Distance2father);
	distance2NearestOTUForRecent = MeanDistanceFromNearestOTU; // Starting point
	for(i = 0, totalBranchLengthAncient=0, totalBranchLengthRecent=0; i<numberOfNodes-1; ++i){
		if(DistanceFromNearestOTU[i] > distance2NearestOTUForRecent)
			totalBranchLengthAncient+= Distance2father[i];
		else
			totalBranchLengthRecent+= Distance2father[i];
	}	
	diffTotalBranchLengthRecentAncient = totalBranchLengthAncient - totalBranchLengthRecent;
	bool isRecentBiggerAncient = true;
	if(totalBranchLengthAncient>totalBranchLengthRecent)
		isRecentBiggerAncient = false;

	bool isImprovedRecentEstimation = true;
	int numberOfIterations = 0;
	while(isImprovedRecentEstimation && (numberOfIterations<100000)){
		MDOUBLE prevDiffTotalBranchLengthRecentAncient = diffTotalBranchLengthRecentAncient; // init
		MDOUBLE prevDistance2NearestOTUForRecent = distance2NearestOTUForRecent;
		MDOUBLE prevtotalBranchLengthAncient = totalBranchLengthAncient;
		MDOUBLE prevtotalBranchLengthRecent = totalBranchLengthRecent;
		distance2NearestOTUForRecent = distance2NearestOTUForRecent + diffTotalBranchLengthRecentAncient/(numberOfNodes*100000); // cont. correction

		for(i = 0, totalBranchLengthAncient=0, totalBranchLengthRecent=0; i<numberOfNodes-1; ++i){
			if(DistanceFromNearestOTU[i] > distance2NearestOTUForRecent)
				totalBranchLengthAncient+= Distance2father[i];
			else
				totalBranchLengthRecent+= Distance2father[i];
		}	
		diffTotalBranchLengthRecentAncient = totalBranchLengthAncient - totalBranchLengthRecent;

		if(abs(diffTotalBranchLengthRecentAncient) > abs(prevDiffTotalBranchLengthRecentAncient) 
			//&& ((totalBranchLengthAncient>totalBranchLengthRecent)*isRecentBiggerAncient)	// to make sure that Ancient is not more than Recent, wait for "flip"
			)
		{
			isImprovedRecentEstimation = false;
			distance2NearestOTUForRecent = prevDistance2NearestOTUForRecent; // go back to last estimation.
			diffTotalBranchLengthRecentAncient = prevDiffTotalBranchLengthRecentAncient;
			totalBranchLengthAncient = prevtotalBranchLengthAncient;
			totalBranchLengthRecent = prevtotalBranchLengthRecent;
		}
		//cout<<diffTotalBranchLengthRecentAncient<<" "<<distance2NearestOTUForRecent<<"\n";	// DEBUG
		numberOfIterations++;
	}
	LOGnOUT(4,<<"The computed distance2NearestOTUForRecent="<<distance2NearestOTUForRecent<<" with TotalBranchLength Ancient="<<totalBranchLengthAncient<<" Recent="<<totalBranchLengthRecent<<" Converged "<<!isImprovedRecentEstimation<<endl);
	return distance2NearestOTUForRecent;
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss::updateSetLofMissingData(){
	if(!gainLossOptions::_gainLossDist)
		_unObservableData_p->setLforMissingData(_tr,_sp);
	else
		_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
}


void gainLoss::multipleAllBranchesByFactorAtStartByMaxParsimonyCost(int costOfTreeMP){
	MDOUBLE branchLengthSum = _tr.getAllBranchesLengthSum();
	MDOUBLE requiredBranchLengthSumByMaxParsimonyCost = (double)costOfTreeMP/_sc.seqLen();
	MDOUBLE factorBL = requiredBranchLengthSumByMaxParsimonyCost / branchLengthSum;
	_tr.multipleAllBranchesByFactor(factorBL);
	MDOUBLE updatedBranchLengthSum = _tr.getAllBranchesLengthSum();
	LOGnOUT(4,<<"  multipleAllBranchesByFactorAtStartByMaxParsimonyCost Total branch lengths: "<<updatedBranchLengthSum <<" with respect to costOfTreeMP "<<costOfTreeMP<<endl);
}

/********************************************************************************************
brent doesn't work if the limits are "non-average-able".
e.g., if min=0.01 and max=100 ave =~50 and brent will not go towards the lower values
Solution: the exponent of log10, thus,
1 = 10^0,
0.1 = 10-1
10 = 10^1
*********************************************************************************************/
void gainLoss::multipleAllBranchesByFactorAtStart(MDOUBLE epsilonOptimization){

	//printTreeLikelihoodAllPosAlphTheSame(); // updates _logL
	MDOUBLE branchLengthSum = _tr.getAllBranchesLengthSum();
	MDOUBLE factorBL = 1;
	LOGnOUT(4,<<" Start multipleAllBranchesByFactorAtStart use epsilonOptimization "<<epsilonOptimization<<" Total branch lengths:"<<branchLengthSum <<endl);
	MDOUBLE minBranchProportionExponent = -8; 
	MDOUBLE maxBranchProportionExponent = 8; 
	MDOUBLE bestBranchProportionExponent = 0;
	MDOUBLE currBranchProportionExponent = 0;
	MDOUBLE currBestL = VERYSMALL;
	MDOUBLE logLimprovement = 0;
	bool isStopAfterNoImprovment = false;

	//while(maxBranchProportionExponent>=0){	// allow up to 8 orders of magnitude change
		LOGnOUT(4,<<"Allow proportion: "<<pow(10,minBranchProportionExponent) <<" to "<<pow(10,maxBranchProportionExponent)<<endl);

		if(gainLossOptions::_gainLossDist)
			currBestL = -brent(minBranchProportionExponent,bestBranchProportionExponent,maxBranchProportionExponent,evalBranchProportionExponentSPvv(&_tr, _sc, _spVVec,_gainDist,_lossDist,NULL,_unObservableData_p),epsilonOptimization,&currBranchProportionExponent);
		else
			currBestL = -brent(minBranchProportionExponent,bestBranchProportionExponent,maxBranchProportionExponent,evalBranchProportionExponent(&_tr, _sc, _sp,NULL,_unObservableData_p),epsilonOptimization,&currBranchProportionExponent);
		factorBL = pow(10,currBranchProportionExponent);
		logLimprovement = currBestL-_logL;
		if(logLimprovement > 0){
			_tr.multipleAllBranchesByFactor(factorBL);
			if(_unObservableData_p){
				if(gainLossOptions::_gainLossDist)
					_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);
				else
					_unObservableData_p->setLforMissingData(_tr,_sp);
			}
			printTreeLikelihoodAllPosAlphTheSame(); // updates _logL
			LOGnOUT(4,<<"Tree multiplied by "<<factorBL<< " Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
			printTree(_tr);
			if(! currBestL > _logL+epsilonOptimization && isStopAfterNoImprovment){
				LOGnOUT(4,<<"Last iteration with maxBranchProportionExponent "<<maxBranchProportionExponent<<" and Likelihood improvement of "<<logLimprovement<< " Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
				//break;
			}
		}
		else{
			LOGnOUT(4,<<" Branch length LengthSum was not changed. factor="<<factorBL <<endl);
			LOGnOUT(4,<<" Total branch lengths:"<<_tr.getAllBranchesLengthSum() <<endl);
			return;
		}
	//	minBranchProportionExponent += 2;
	//	maxBranchProportionExponent -= 2;
	//}
}



/********************************************************************************************
numOfSequenceSets = original seq. length (~4,800)
numOfRepeats = replications of simulation (~10)
The simulations if by simulateOnePos (not efficient) 
- for each positions we randoms sample gain and loss rates and create new sp

There are few options for data simulation:
gain,loss rates are sampled from (for each position):
1. Gamma distributions [with empirical parameters]					//_initParamsFromTrueEstimation
!!! Need to run with previously found user parameters: Tree, _userTheta,_userAlphaGain,_userBetaGain,_userAlphaLoss,_userBetaLoss
2. MP estimated rates (the empirical _MPPerPos) +minimalRate added	//initParamsFromMPEstimation
3. uniform distributions											//_initParamsAtRandPointsInSimPostExp, Default

Theta ("1" freq):
1. taken from the empirical theta for all positions					//_initParamsFromTrueEstimation
2. The observed frequencies for all positions+perturbation			//isTheataFromObservedFreq
2. sampled from uniform distribution
a. once for all positions										//_initParamsAtRandPointsInSimPostExp, Default
b. for each position											//_initRootFreqAtRandPointsInSimPostExpEachPos
*********************************************************************************************/
void gainLoss::startSimultePosteriorExpectationOfChange(int numOfSequenceSets, const int numOfRepeats)
{
	bool isNormalizeQAfterRatesSample = true;			// After Norm, multi by gain+loss
	bool isNormalizeQwithEmpricialQ = true;	// applicable for MPemt and SMest
	bool isMultBy2_normQ = false;					// False, old and wrong correction, no need.
	bool isComputeEmpiricalCorrection = false;	// Failed trial. MP and SM, sampling rates
	bool isGammaRatioAdjusted = false;	// thus, the user data is override
	bool isThetaSampledForGamma = false;
	bool isThetaFromObservedForEmpiricalSimulations = true; // for MPest and SMest
	bool isMPcostEmpirical = false; //if false, the loss2gainRatioToSim
	bool isRateEQnumOfEvents = true;	// gain=#gainEvents
	bool isUsePeudoCountForEmpirical = true;	// keep it true, it's actually minRate as before, but with correct sampling
	MDOUBLE minPeudoCountForEmpirical = 0.01;
	
	MDOUBLE glRatioTieBreakerInCostMatrix = 0.0;	//is positive, losses are favored (gain cost is higher)
	MDOUBLE epsilonForgainLossRatio = 0.01;		// was 0.1
	MDOUBLE loss2gainRatioToSim = gainLossOptions::_loss2gainRatioToSim;
	MDOUBLE MaxGLratio = 2.1; //was 37, here, only 1,2 is accounted for in multiple MP costs

	// Evolutionary model based, with separate variables (possibly "flat" each replication)
	stochasticProcess* spSim=NULL;
	stochasticProcess* spSimpleSim=NULL;
	vector<vector<stochasticProcess*> > spVVecSim;
	unObservableData* unObservableDataSim=NULL;
	distribution* gainDistSim=NULL;
	distribution* lossDistSim=NULL;		
	tree trSim;
	VVdouble LpostPerCatSim; // the posterior probability for each position for each category
	VVVdouble LpostPerSpPerCatSim;

	MDOUBLE minThetaRandSample = 0.1;	// was 0.01. change all from 0.01 to 0.05, and later to 0.1
	MDOUBLE maxThetaRandSample = 0.9;	// was 0.09
	MDOUBLE minGainRandSample = 0.1;	// was 0.01
	MDOUBLE maxGainRandSample = 2.0;	// was 2.5, now E(val) = 1	
	MDOUBLE minLossRandSample = 0.1;	// was 0.01
	MDOUBLE maxLossRandSample = loss2gainRatioToSim*2;
	MDOUBLE meanGaussianGain = 1.0;
	MDOUBLE varianceGaussianGain = 1.0;
	MDOUBLE minAllowedRate = 0.01;	// 0.01, An important parameter. used to avoid too low or high rates in Gamma and MP
	MDOUBLE maxAllowedRate = 100;

	MDOUBLE meanGainFromEMP=1;
	MDOUBLE meanLossFromEMP=1;
	MDOUBLE meanQrateFromEMP=1;
	// these parameter need to be part of gainLossOptions
	bool printTreeForEachReplication = true;
	//bool isUseMeanEventFromEMP = true;	// sum of gain and loss, if T: meanGainFromMP=meanLossFromMP=Events (for computation)
	MDOUBLE meanEventsFromEMP=1;
	MDOUBLE expectedQvalEmpirical=1;

	MDOUBLE meanGaussianLoss = loss2gainRatioToSim;
	MDOUBLE varianceGaussianLoss = loss2gainRatioToSim;
	MDOUBLE Theta = gainLossOptions::_userTheta; // 
	MDOUBLE AlphaGain = gainLossOptions::_userAlphaGain; //
	MDOUBLE BetaGain = gainLossOptions::_userBetaGain; // 
	MDOUBLE AlphaLoss = gainLossOptions::_userAlphaLoss; // 
	MDOUBLE BetaLoss = gainLossOptions::_userBetaLoss; // 
	MDOUBLE AlphaRate =  gainLossOptions::_userAlphaRate; //
	if(gainLossOptions::_performParametricBootstapCorrelation){
		isNormalizeQAfterRatesSample = false;
		if(gainLossOptions::_gainLossDist){
			Theta =static_cast<gainLossModel*>((*_spVVec[0][0]).getPijAccelerator()->getReplacementModel())->getTheta(); //  gainLossOptions::_userTheta
			AlphaGain = getRateAlpha(_gainDist); // gainLossOptions::_userAlphaGain
			BetaGain = getRateBeta(_gainDist); // gainLossOptions::_userBetaGain
			AlphaLoss = getRateAlpha(_lossDist); // gainLossOptions::_userAlphaLoss
			BetaLoss = getRateBeta(_lossDist); // gainLossOptions::_userBetaLoss
		}else{
			Theta =static_cast<gainLossModel*>(_sp->getPijAccelerator()->getReplacementModel())->getTheta();
			AlphaRate =  getRateAlpha(_sp->distr());
		}
		minGainRandSample = 0.01;	
		maxGainRandSample = VERYBIG;		
		minLossRandSample = 0.01;	
		maxLossRandSample = VERYBIG;
		minAllowedRate = 0.01;
		maxAllowedRate = VERYBIG;
	}

	MDOUBLE gainPlusLossExpectancyGamma = (AlphaGain/BetaGain)+(AlphaLoss/BetaLoss);

	MDOUBLE costMatrixGainLossRatio = gainLossOptions::_costMatrixGainLossRatio; // to be updated according to simulation
	MDOUBLE costMatrixGainLossRatioCorrectionFactor =1;
	MDOUBLE minAllowedMeanEMP = 0.01;
	bool normalizationFactorForLoss1AsInTree = false;
	MDOUBLE randomNoise =0;

	Vdouble freq(2,0.0);
	MDOUBLE init_gain = 0.5; //gainLossOptions::_userGain taken from original runs of COG data, get it from params file
	MDOUBLE init_loss = 0.5; //gainLossOptions::_userLoss
	MDOUBLE rateSample = 1;
	MDOUBLE lossGainRatioSample = 1;
	bool _isHGT_normal_Pij = gainLossOptions::_isHGT_normal_Pij;
	bool _isHGT_with_Q = gainLossOptions::_isHGT_with_Q;

	// DEBUG Test for gain events in Eq sequences (change isTestForGainEventsInEqSeq=true)
	bool isTestForGainEventsInEqSeq =false;

	LOGnOUT(4,<<endl<<"****************************************************\n startSimultePosteriorExpectationOfChange... "<<endl);
	LOGnOUT(4,<<"Replicates="<<numOfRepeats<<" Positions="<<numOfSequenceSets<<endl);
	LOGnOUT(4,<<" simulationType {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise;..."<<endl);
	LOGnOUT(4,<<" EQ_gEql,EQ_gVrl,Gam_gEql,GamgVrl}="<<gainLossOptions::_simulationType<<endl);

	LOGnOUT(4,<<" loss/gain = ");
	if(!(gainLossOptions::_simulationType == gainLossOptions::Gamma))
		LOGnOUT(4,<<loss2gainRatioToSim<<endl)
	else
	LOGnOUT(4,<<(gainLossOptions::_userAlphaLoss/gainLossOptions::_userBetaLoss) /  (gainLossOptions::_userAlphaGain/gainLossOptions::_userBetaGain)<<endl);
	if(!gainLossOptions::_isRootFreqEQstationaryInSimulations){
		if(gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos)
			LOGnOUT(4,<<"Root(1) freq is sampled seperatly for each pos"<<endl)
		else
		LOGnOUT(4,<<"Root(1) freq is sampled once for the entire replication (sim all positions)"<<endl)
		if(isComputeEmpiricalCorrection 
			&& (gainLossOptions::_simulationType == gainLossOptions::MPestEmp || gainLossOptions::_simulationType == gainLossOptions::SMestEmp  ))
			LOGnOUT(3,<<"!!! WARN !!! the ComputeEmpiricalCorrection for SM and MP assumes RootFreq=Stationary"<<endl);
	}
	else
		LOGnOUT(4,<<"Root(1) freq = gain/(gain+loss)"<<endl);	

	if(isGammaRatioAdjusted){
		LOGnOUT(4,<<" LossBeta is assigned to maintain  loss2gainRatioToSim="<<loss2gainRatioToSim<<endl);
		BetaLoss = BetaGain*AlphaLoss/(AlphaGain*loss2gainRatioToSim);
	}
	LOGnOUT(4,<<" All rate simulations are done with minAllowedRate "<<minAllowedRate<<" and maxAllowedRate "<<maxAllowedRate<<endl);

	if(gainLossOptions::_isRootFreqEQstationaryInSimulations & !gainLossOptions::_isRootFreqEQstationary)
		LOGnOUT(3,<<"\n\n WARNING!!! _isRootFreqEQstationaryInSimulations " <<gainLossOptions::_isRootFreqEQstationaryInSimulations<<" and model "<<gainLossOptions::_isRootFreqEQstationary<<endl<<endl<<endl);

	if(gainLossOptions::_simulationType == gainLossOptions::MPestEmp){
		if(isNormalizeQwithEmpricialQ)
			isNormalizeQAfterRatesSample = false;	// in this case, no other normalization
		if(_MPPerPos.size()==0){
			LOGnOUT(4,<<" _MPPerPos size="<<_MPPerPos.size()<<endl);
			startMaxParsimonyChange();
		}
		meanGainFromEMP = max(getVMatrixJK(_MPPerPos,0,1)/_MPPerPos.size(), minAllowedMeanEMP);	// used for initParamsFromMPEstimation
		meanLossFromEMP = max(getVMatrixJK(_MPPerPos,1,0)/_MPPerPos.size(), minAllowedMeanEMP);	// used for initParamsFromMPEstimation
		MDOUBLE sumQrateFromEMP = 0;
		for(int i=0; i<_MPPerPos.size(); ++i){
			MDOUBLE gain = _MPPerPos[i][0][1];
			MDOUBLE loss = _MPPerPos[i][1][0];
			if(gain+loss>0)
				sumQrateFromEMP += 2*gain*loss/(gain+loss);
		}
		meanQrateFromEMP = sumQrateFromEMP/_MPPerPos.size();
		if(isUsePeudoCountForEmpirical){
			LOGnOUT(4,<<"To avoid zero rates, Use pseudo counts=  "<<minPeudoCountForEmpirical<<endl);
			meanGainFromEMP += minPeudoCountForEmpirical;
			meanLossFromEMP += minPeudoCountForEmpirical;
		}
		meanEventsFromEMP = meanGainFromEMP+meanLossFromEMP;		
		LOGnOUT(4,<<" The mean number by MP of gain= "<<meanGainFromEMP<<", loss= "<<meanLossFromEMP<<endl);
		if(normalizationFactorForLoss1AsInTree){ // not is use
			meanGainFromEMP = meanGainFromEMP/meanLossFromEMP;
			meanLossFromEMP = meanLossFromEMP/meanLossFromEMP;	// results in 1
			LOGnOUT(4,<<" The mean number by MP after normalization of Loss=1 (compatible with tree). gain= "<<meanGainFromEMP<<", loss= "<<meanLossFromEMP<<endl);
		}
		if(isComputeEmpiricalCorrection){
			expectedQvalEmpirical=ComputeEmpiricalExpectedQforStationaryProcess(_MPPerPos,minPeudoCountForEmpirical);
			cout<<"expectedQvalEmpirical="<<expectedQvalEmpirical<<endl;
		}
	}
	if(gainLossOptions::_simulationType == gainLossOptions::SMestEmp){
		if(isNormalizeQwithEmpricialQ)
			isNormalizeQAfterRatesSample = false;	// in this case, no other normalization

		if(_SMPerPos.size()==0){
			LOGnOUT(4,<<" _SMPerPos size="<<_SMPerPos.size()<<endl);
			startComputePosteriorExpectationOfChange();
		}
		meanGainFromEMP = max(getVMatrixJK(_SMPerPos,0,1)/_SMPerPos.size(), minAllowedMeanEMP);	// used for initParamsFromMPEstimation
		meanLossFromEMP = max(getVMatrixJK(_SMPerPos,1,0)/_SMPerPos.size(), minAllowedMeanEMP);	// used for initParamsFromMPEstimation
		MDOUBLE sumQrateFromEMP = 0;
		for(int i=0; i<_SMPerPos.size(); ++i){
			MDOUBLE gain = _SMPerPos[i][0][1];
			MDOUBLE loss = _SMPerPos[i][1][0];
			if(gain+loss>0)
				sumQrateFromEMP += 2*gain*loss/(gain+loss);
		}
		meanQrateFromEMP = sumQrateFromEMP/_SMPerPos.size();
		if(isUsePeudoCountForEmpirical){
			LOGnOUT(4,<<"To avoid zero rates, Use pseudo counts=  "<<minPeudoCountForEmpirical<<endl);
			meanGainFromEMP += minPeudoCountForEmpirical;
			meanLossFromEMP += minPeudoCountForEmpirical;
		}
		meanEventsFromEMP = meanGainFromEMP+meanLossFromEMP;
		LOGnOUT(4,<<" The mean number by SM of gain= "<<meanGainFromEMP<<", loss= "<<meanLossFromEMP<<endl);		
		if(normalizationFactorForLoss1AsInTree){ // not is use
			meanGainFromEMP = meanGainFromEMP/meanLossFromEMP;
			meanLossFromEMP = meanLossFromEMP/meanLossFromEMP;	// results in 1
			LOGnOUT(4,<<" The mean number by SM after normalization of Loss=1 (compatible with tree). gain= "<<meanGainFromEMP<<", loss= "<<meanLossFromEMP<<endl);
		}
		if(isComputeEmpiricalCorrection){
			expectedQvalEmpirical=ComputeEmpiricalExpectedQforStationaryProcess(_SMPerPos,minPeudoCountForEmpirical);
			cout<<"expectedQvalEmpirical="<<expectedQvalEmpirical<<endl;
		}
	}
	if(gainLossOptions::_isRootFreqEQstationaryInSimulations){
		LOGnOUT(4,<<" Theta (Root freq of '1's) is equal to the stationary one"<<endl);
	}
	else if(gainLossOptions::_isTheataFromObservedFreq  /*&& !(gainLossOptions::_simulationType==gainLossOptions::Gamma)*/){
		Vdouble observedFreq = evaluateCharacterFreq(_scWithFullLength);
		Theta = observedFreq[1];
		MDOUBLE maxDiffFromObservedFreq = 0.0; // 0.2
		maxThetaRandSample = min(maxThetaRandSample,Theta+maxDiffFromObservedFreq) ;
		minThetaRandSample = max(minThetaRandSample, Theta-maxDiffFromObservedFreq);
		LOGnOUT(4,<<" Theta taken from 'counting' freq= "<<observedFreq[1]<<" with random perturbation of "<<maxDiffFromObservedFreq<<endl);
	}// else it's from _userTheta

	if(isTestForGainEventsInEqSeq){
		numOfSequenceSets = 1000;
		LOGnOUT(4,<<endl<<"Using TestForGainEventsInEqSeq with numOfSequenceSets= "<<numOfSequenceSets<<endl);
	}

	//string treeFile = gainLossOptions::_treeFile; // input tree - same for all iterations
	if(gainLossOptions::_isRootFreqEQstationaryInSimulations)
		LOGnOUT(4,<<"\tIn statrionary model - Theta=Root(1) is driven from the gain/gain+loss"<<endl);
	if(gainLossOptions::_isMPratio)
		LOGnOUT(4,<<"\tMPratio simulations: gain is sampled, and loss if multiplied by the MPcost= "<<costMatrixGainLossRatio<<endl);
	switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
	{
	case gainLossOptions::Uniform:
		LOGnOUT(4,<<"\n\n(*) UNIFORM simulations: Using RandSample parameters:"<<endl);
		LOGnOUT(4,<<"\tGain Min="<<minGainRandSample<<" Max="<<maxGainRandSample<<" "<<endl);
		if(!gainLossOptions::_isMPratio)
			LOGnOUT(4,<<"\tLoss Min="<<minLossRandSample<<" Max="<<maxLossRandSample<<" "<<endl);
		if(!gainLossOptions::_isRootFreqEQstationaryInSimulations)
			LOGnOUT(4,<<"\tTheta Min="<<minThetaRandSample<<" Max="<<maxThetaRandSample<<endl);
		break;
	case gainLossOptions::Normal:
		LOGnOUT(4,<<"\n\n(*) Gaussian simulations: Using RandSample parameters:"<<endl);
		LOGnOUT(4,<<"\tGain mean="<<meanGaussianGain<<" var="<<varianceGaussianGain<<" "<<endl);
		if(!gainLossOptions::_isMPratio)
			LOGnOUT(4,<<"\tLoss mean="<<meanGaussianLoss<<" var="<<varianceGaussianLoss<<" "<<endl);
		break;
	case gainLossOptions::Gamma:
		LOGnOUT(4,<<"\n\n(*) GAMMA simulations: Sample from Gamma with parameters:\n");
		if(isThetaSampledForGamma){
			LOGnOUT(4,<<"Note: Theta is sampled for Gamma, not taken from userTheta\n"<<endl);
			LOGnOUT(4,<<"Theta Min="<<minThetaRandSample<<" Max="<<maxThetaRandSample); 		}
		else if (!gainLossOptions::_isRootFreqEQstationaryInSimulations)
			LOGnOUT(4,<<"\n_userTheta="<<gainLossOptions::_userTheta);
		LOGnOUT(4,<<"\n_userAlphaGain=\t"<<AlphaGain
			<<"\n_userBetaGain=\t"<<BetaGain<<endl);
		if(!gainLossOptions::_isMPratio)
			LOGnOUT(4,<<"_userAlphaLoss=\t"<<AlphaLoss
			<<"\n_BetaLoss=\t"<<BetaLoss<<endl);
		if(isGammaRatioAdjusted)
			LOGnOUT(4,<<"Note: Gamma is RatioAdjusted, the BetaLoss is determined by the required gain:loss ratio "<<loss2gainRatioToSim<<endl);
		loss2gainRatioToSim = (AlphaLoss/BetaLoss)/(AlphaGain/BetaGain);
		break;
	case gainLossOptions::MPestEmp:
		LOGnOUT(4,<<"\n\n(*) MP simulations: Using Maximum parsimony empirical estimated gain and loss rates \n sampled from "<<_MPPerPos.size()<<" positions "<<endl);
		loss2gainRatioToSim = meanLossFromEMP/meanGainFromEMP;
		//if(isUseMeanEventFromEMP)
		//	LOGnOUT(4,<<" Note: meanGainFromMP = meanLossFromMP = Events (Thus, sampleGain/all and sampleLoss/all) "<<endl);			
		if(isThetaFromObservedForEmpiricalSimulations && !gainLossOptions::_isRootFreqEQstationaryInSimulations)
			LOGnOUT(4,<<"Note: for the Empirical simulation - theta is taken from Observed freq"<<endl);
		break;
	case gainLossOptions::SMestEmp:
		LOGnOUT(4,<<"\n\n(*) SM simulations: Using stochstic mapping empirical estimated gain and loss rates \n sampled from "<<_SMPerPos.size()<<" positions "<<endl);
		loss2gainRatioToSim = meanLossFromEMP/meanGainFromEMP;
		//if(isUseMeanEventFromEMP)
		//	LOGnOUT(4,<<" Note: meanGainFromEMP = meanLossFromEMP = Events (Thus, sampleGain/all and sampleLoss/all) "<<endl);
		if(isThetaFromObservedForEmpiricalSimulations && !gainLossOptions::_isRootFreqEQstationaryInSimulations)
			LOGnOUT(4,<<"Note: for the Empirical simulation - theta is taken from Observed freq"<<endl);
		break;
	case gainLossOptions::GammaNoise:
		LOGnOUT(4,<<"\n\n(*) GAMMA simulations with noise level "<< gainLossOptions::_noiseLevelInGammaSimulation<<" parameter before noise:");
		LOGnOUT(4,<<"\n_userTheta="<<gainLossOptions::_userTheta);
		LOGnOUT(4,<<"\n_userAlphaGain="<<AlphaGain
			<<"\n_userBetaGain="<<BetaGain<<endl);
		break;
	case gainLossOptions::EQ_gEql:
		LOGnOUT(4,<<" Simulation - EQ rate, gain=loss (via stationary freq)"<<endl);
		loss2gainRatioToSim = 1.0;
		break;
	case gainLossOptions::EQ_gVrl:
		LOGnOUT(4,<<" Simulation - EQ rate, gain/loss ratio unif variable. Mean loss/gain="<<loss2gainRatioToSim<<" epsilon="<< epsilonForgainLossRatio<<endl);	
		break;
	case gainLossOptions::Gam_gEql:
		LOGnOUT(4,<<" Simulation - Gamma rate, alpha="<<AlphaRate<<", gain=loss (via stationary freq)"<<endl);
		loss2gainRatioToSim = 1.0;
		break;
	case gainLossOptions::Gam_gVrl:
		LOGnOUT(4,<<" Simulation - Gamma rate, alpha="<<AlphaRate<<", gain/loss ratio unif variable. Mean loss/gain="<<loss2gainRatioToSim<<" epsilon="<< epsilonForgainLossRatio<<endl);			
		break;
	default:
		errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp GammaNoise, MPratio}");
	}
	if(isNormalizeQwithEmpricialQ
		&& (gainLossOptions::_simulationType == gainLossOptions::MPestEmp || gainLossOptions::_simulationType == gainLossOptions::SMestEmp) )
		LOGnOUT(4,<<" Q matrix is normalized given empirical expected Q values."<<endl);
	if(isNormalizeQAfterRatesSample)
		LOGnOUT(4,<<" Q matrix is normalized after sampling."<<endl);			

	if(isMultBy2_normQ 
		&& !(gainLossOptions::_simulationType == gainLossOptions::Gamma) 
		&& !(gainLossOptions::_simulationType == gainLossOptions::SMestEmp) 
		&& !(gainLossOptions::_simulationType == gainLossOptions::MPestEmp) )	// with mult=2, Q matrix is normalized with respect to the tree (after multiplied by freq=0.5)
		LOGnOUT(4,<<" Gain and loss rates multiplied by 2. Mainitaining normalized Q matrix."<<endl);

	
	////////////////////////////////////////////////////////////////////////// Replicates	
	for(int replicat=1; replicat<=numOfRepeats; ++replicat){
		LOGnOUT(4,<<endl<<".......................................Replicate= "<<replicat<<endl);
		time_t t1,t2;
		time(&t1);

		createDir(gainLossOptions::_outDir, "SimulatedPostExp"+ int2string(replicat));
		string outDirSeq = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seqAll" ;
		createDir(gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat), "seqAll" );
		string simulatedEventsSimString = outDirSeq + "//" + "simulatedEvents.txt";
		ofstream* simulatedEventsFile = new ofstream(simulatedEventsSimString.c_str());
		string posSim = outDirSeq + "//" + "nodesContentSim" + ".txt";
		ofstream* posSim_out = new ofstream(posSim.c_str());



		string perPosStat = outDirSeq + "//" + "statPos.txt";
		ofstream perPosStatStream(perPosStat.c_str());
		perPosStatStream<<"pos"<<"\t"<<"rate"<<"\t"<<"theta"<<"\t"<<"occur"<<"\n";

		string perBranchStat = outDirSeq + "//" + "statBranch.txt";
		ofstream perBranchStatStream(perBranchStat.c_str());
		perBranchStatStream<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<endl;
		treeIterTopDownConst tit(_tr);
		for (tree::nodeP myN = tit.first();myN!=tit.end(); myN = tit.next()) {
			if(myN->isRoot())
				continue;
			perBranchStatStream<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()<<endl;
		}
		perBranchStatStream.close();
		MDOUBLE init_gainsForCostMatrix = 0.0; // sum all position
		MDOUBLE init_lossesForCostMatrix = 0.0; // sum all position
		MDOUBLE init_losses2gainRatioForCostMatrixSum = 0.0;
		//MDOUBLE QnormTest = 0.0;

		// produce random noise
		if(gainLossOptions::_simulationType == gainLossOptions::GammaNoise){
			randomNoise = talRandom::giveRandomNumberBetweenTwoPoints(-gainLossOptions::_noiseLevelInGammaSimulation, gainLossOptions::_noiseLevelInGammaSimulation);
			// if noiseLevel=200% than param may be up to x3 or down to x0.33 its value
			if(randomNoise>=0)
				randomNoise = 1+randomNoise;
			else
				randomNoise = 1/(1-randomNoise);
			LOGnOUT(4,<<"Noise over all parameters="<< randomNoise<<endl);
		}		
		// Theta for all positions, (not relevant to stationary models)
		if(!gainLossOptions::_isRootFreqEQstationaryInSimulations){ // else Theta is driven from gain/gain+loss
			switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
			{
			case gainLossOptions::Uniform:
			case gainLossOptions::Normal:
				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
				break;
			case gainLossOptions::MPestEmp:
			case gainLossOptions::SMestEmp:
				if(isThetaFromObservedForEmpiricalSimulations)
					freq[1]=Theta;
				else
					freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
				break;
			case gainLossOptions::Gamma:
				if(isThetaSampledForGamma)
					freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
				else
					freq[1]= gainLossOptions::_userTheta;
				break;
			case gainLossOptions::GammaNoise:
				freq[1] = gainLossOptions::_userTheta*randomNoise;
				freq[1] = max(freq[1],minThetaRandSample);	// added to avoid too small or too big theta
				freq[1] = min(freq[1],maxThetaRandSample);
				break;
			case gainLossOptions::EQ_gEql:
			case gainLossOptions::EQ_gVrl:
			case gainLossOptions::Gam_gEql:
			case gainLossOptions::Gam_gVrl:
				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
				break;
			default:
				errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise}");
			}
			if(!gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos)
				LOGnOUT(4,<<" For all positions, Root(1)= "<<freq[1]<<endl);
		}
		else{
			LOGnOUT(4,<<" Stationary model - Theta=Root(1) is driven from the gain/gain+loss"<<endl);
		}

		vector<bool> isGainEventInAnode;	// DEBUG
		if(isTestForGainEventsInEqSeq)
			isGainEventInAnode.resize(numOfSequenceSets+1);

		
		////////////////////////////////////////////////////////////////////////// Positions
		MDOUBLE ratePerPosSum = 0;
		int randomPosition; //used in MPestEmp or SMestEmp
		sequenceContainer seqSimulated;
		gainLossAlphabet alph;


		for(int i=0; i<numOfSequenceSets; ++i)
		{			
			rateSample = 1;
			lossGainRatioSample = 1;
			switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
			{
			case gainLossOptions::Uniform:
				init_gain = talRandom::giveRandomNumberBetweenTwoPoints(minGainRandSample, maxGainRandSample);
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else
					init_loss = talRandom::giveRandomNumberBetweenTwoPoints(minLossRandSample, maxLossRandSample);
				break;
			case gainLossOptions::Normal:
				init_gain = talRandom::rand_gaussian(meanGaussianGain, varianceGaussianGain);
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else
					init_loss = talRandom::rand_gaussian(meanGaussianLoss, varianceGaussianLoss);					
				break;
			case gainLossOptions::Gamma:
				init_gain = talRandom::SampleGamma(AlphaGain,BetaGain);
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else					
					init_loss = talRandom::SampleGamma(AlphaLoss,BetaLoss);					
				break;
			case gainLossOptions::MPestEmp:
				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _MPPerPos.size());
				if(isRateEQnumOfEvents)
					init_gain = _MPPerPos[randomPosition][0][1];
				else{
					if(isUsePeudoCountForEmpirical)
						init_gain = (_MPPerPos[randomPosition][0][1]+minPeudoCountForEmpirical)/(meanGainFromEMP*2) *sqrt(meanGainFromEMP/meanLossFromEMP) ; // was /meanEventsFromEMP
					else
						init_gain = _MPPerPos[randomPosition][0][1]/meanEventsFromEMP; // was /meanEventsFromEMP
				}				
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else{
					if(isRateEQnumOfEvents)
						init_loss = _MPPerPos[randomPosition][1][0];
					else{
						if(isUsePeudoCountForEmpirical)
							init_loss = (_MPPerPos[randomPosition][1][0]+minPeudoCountForEmpirical)/(meanLossFromEMP*2)*sqrt(meanLossFromEMP/meanGainFromEMP) ; // was /meanEventsFromEMP
						else
							init_loss = _MPPerPos[randomPosition][1][0]/meanEventsFromEMP ; // was /meanEventsFromEMP
					}
				}
				break;
			case gainLossOptions::SMestEmp:
				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _SMPerPos.size());
				if(isRateEQnumOfEvents)
					init_gain = _SMPerPos[randomPosition][0][1];
				else{
					if(isUsePeudoCountForEmpirical)
						init_gain = (_SMPerPos[randomPosition][0][1]+minPeudoCountForEmpirical)/(meanGainFromEMP*2) *sqrt(meanGainFromEMP/meanLossFromEMP) ; // was /meanEventsFromEMP
					else
						init_gain = _SMPerPos[randomPosition][0][1]/meanEventsFromEMP; // was /meanEventsFromEMP
				}				
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else{
					if(isRateEQnumOfEvents)
						init_loss = _SMPerPos[randomPosition][1][0];
					else{
						if(isUsePeudoCountForEmpirical)
							init_loss = (_SMPerPos[randomPosition][1][0]+minPeudoCountForEmpirical)/(meanLossFromEMP*2)*sqrt(meanLossFromEMP/meanGainFromEMP) ; // was /meanEventsFromEMP
						else
							init_loss = _SMPerPos[randomPosition][1][0]/meanEventsFromEMP ; // was /meanEventsFromEMP
					}
				}
				break;
			case gainLossOptions::GammaNoise:
				init_gain = talRandom::SampleGamma( (AlphaGain*randomNoise)
					,(BetaGain*randomNoise));
				if(gainLossOptions::_isMPratio)
					init_loss = init_gain * costMatrixGainLossRatio;
				else{					
					init_loss = talRandom::SampleGamma((AlphaLoss*randomNoise)
						,(BetaLoss*randomNoise));
				}
				break;
			case gainLossOptions::EQ_gEql:
				init_gain = -(rateSample/(-1-lossGainRatioSample));			//init_gain = init_loss =0.5;
				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
				break;
			case gainLossOptions::Gam_gEql:
				rateSample = talRandom::SampleGamma(AlphaRate);	//init_gain = init_loss = 0.5*rateSample;
				init_gain = -(rateSample/(-1-lossGainRatioSample));
				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
				break;
			case gainLossOptions::EQ_gVrl:
				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
				init_gain = -(rateSample/(-1-lossGainRatioSample));
				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
				break;
			case gainLossOptions::Gam_gVrl:
				rateSample = talRandom::SampleGamma(AlphaRate);
				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
				init_gain = -(rateSample/(-1-lossGainRatioSample));
				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
				//init_gain = -(1/(-1-lossGainRatioSample))*rateSample;
				//init_loss = 1+(1/(-1-lossGainRatioSample))*rateSample;
				break;				
			default:
				errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp GammaNoise}");
			}
			init_gain = min(maxAllowedRate, max(init_gain,minAllowedRate));	// added to avoid too small gain rate
			init_loss = min(maxAllowedRate,max(init_loss,minAllowedRate));
			if(isMultBy2_normQ 
				&& !(gainLossOptions::_simulationType == gainLossOptions::Gamma) 
				&& !(gainLossOptions::_simulationType == gainLossOptions::SMestEmp) 
				&& !(gainLossOptions::_simulationType == gainLossOptions::MPestEmp)){	// with mult=2, Q matrix is normalized with respect to the tree (after multiplied by freq=0.5)
				init_gain *=2;
				init_loss *=2;
			}
			///////////// Theta random per pos 
			if(gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos && !gainLossOptions::_isRootFreqEQstationaryInSimulations){	
				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
			}			
			if(gainLossOptions::_isRootFreqEQstationaryInSimulations){	//Theta=Root(1) is driven from the gain/gain+loss
				freq[1]= init_gain/(init_gain+init_loss);
			}
			freq[0]= 1 - freq[1];
			gainLossModelNonReversible  glm(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,_isHGT_normal_Pij,_isHGT_with_Q);
			trivialAccelerator pijAcc(&glm);
			uniDistribution uniDistr;
			stochasticProcess *spSimSingle = NULL;
			spSimSingle = new stochasticProcess(&uniDistr,&pijAcc,false);
			MDOUBLE sumQii = 1.0;
			
			if(isNormalizeQAfterRatesSample ){	// if normalizeQ for each position, there is no rate variability in practice
				sumQii = normalizeQ(spSimSingle);	// (1) Normalize
				MDOUBLE scalingParameterExpectancyOfOne = init_gain+init_loss; //init_gain+init_loss
				if(gainLossOptions::_simulationType == gainLossOptions::Gamma)
					scalingParameterExpectancyOfOne /=gainPlusLossExpectancyGamma;
				static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->norm((scalingParameterExpectancyOfOne));	// (2) multiply by g+l
			}
			if(isNormalizeQwithEmpricialQ
				&& (gainLossOptions::_simulationType == gainLossOptions::MPestEmp || gainLossOptions::_simulationType == gainLossOptions::SMestEmp) )
				static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->norm((1/meanQrateFromEMP));
			if(isComputeEmpiricalCorrection 
				&& (gainLossOptions::_simulationType == gainLossOptions::MPestEmp || gainLossOptions::_simulationType == gainLossOptions::SMestEmp) ){
					static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->norm(1/expectedQvalEmpirical);					
			}
//////////////////////////////////////////////////////////////////////////
			//MDOUBLE gGLM = static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->getMu1();
			//MDOUBLE lGLM = static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->getMu2();
			//MDOUBLE freq1 = static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel())->getTheta();
			//MDOUBLE sumPijQijGLM=(static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel()))->sumPijQij();
			//MDOUBLE rateGLM = gGLM*(1-freq1)+lGLM*(freq1);

			//MDOUBLE gFormula = (1+lossGainRatioSample)/(2*lossGainRatioSample);
			//MDOUBLE lFormula = gFormula*lossGainRatioSample;
			//MDOUBLE rateFormula = (2*gFormula*lFormula)/(gFormula+lFormula);
			//cout<<gGLM<<"\t"<<gFormula<<"\t"<<gGLM-gFormula<<endl;
			//cout<<rateGLM<<"\t"<<rateFormula<<"\t"<<rateGLM-rateFormula<<endl;			
//////////////////////////////////////////////////////////////////////////
			//QnormTest +=  init_gain*freq[0]+init_loss*freq[1];
			init_losses2gainRatioForCostMatrixSum += init_loss/init_gain;
			init_gainsForCostMatrix += init_gain;
			init_lossesForCostMatrix += init_loss;

			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
			//string resFile = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "resSim" + int2string(i+1) + ".sim";

			simulateOnePos *simulateOnePosObj = NULL;
			MDOUBLE ratePerPos=0;
			if(gainLossOptions::_is3states){
				Vdouble init_cpN_vals(4);
				init_cpN_vals[0]=gainLossOptions::_3statesGain; //gain (0->1)
				init_cpN_vals[1]=gainLossOptions::_3statesMore; //more (1->more)
				init_cpN_vals[2]=gainLossOptions::_3statesLess; // less (more->1) 
				init_cpN_vals[3]=gainLossOptions::_3statesLoss; // loss (1->0)
				Vdouble freq_cpN(3);
				freq_cpN[0]=gainLossOptions::_3states0;
				freq_cpN[1]=gainLossOptions::_3states1;
				freq_cpN[2]=1 - (freq_cpN[0] + freq_cpN[1]);
				simulateOnePosObj = new simulateOnePos(strSeqNum, posSim_out, simulatedEventsFile, i,gainLossOptions::_treeFile,init_cpN_vals[0]+init_cpN_vals[3],freq[1],gainLossOptions::_is3states,NULL,&_tr,&init_cpN_vals,&freq_cpN);
			}
			else{
				ratePerPos=(static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel()))->sumPijQij();
				simulateOnePosObj = new simulateOnePos(strSeqNum, posSim_out, simulatedEventsFile, i,gainLossOptions::_treeFile,ratePerPos,freq[1],gainLossOptions::_is3states,spSimSingle,&_tr);
			}
			ratePerPosSum+=ratePerPos;
			perPosStatStream<<i+1<<"\t"<<ratePerPos<<"\t"<<freq[1]<<"\t"<<simulateOnePosObj->getOccurFraction()<<"\n";

			if(spSimSingle) delete spSimSingle;
			if(simulateOnePosObj) 	delete simulateOnePosObj;
			if(isTestForGainEventsInEqSeq){	// DEBUG
				if(simulateOnePosObj->getChangesForBranch(2)[0][1]>0)	// "A" == 2
					isGainEventInAnode[i+1] = true;
			}			
			//if(i==0){	
			//	seqSimulated = sequenceContainer(simulateOnePosObj->getSequenceContainer(),&alph);
			//}
			//else{
			//	sequenceContainer tempSeq = sequenceContainer(simulateOnePosObj->getSequenceContainer(),&alph);
			//	seqSimulated.concatenate(tempSeq);
			//	fastaFormat::write(cout,seqSimulated);
			//}

		}
		if(gainLossOptions::_isMatrixGainLossFromRatioInSimulations) // e.g., val=2, loss rate is double that of loss
			costMatrixGainLossRatio = init_lossesForCostMatrix/init_gainsForCostMatrix;						

		//LOGnOUT(5,<<"QnormTest=\t"<<QnormTest/numOfSequenceSets<<"\n");		
		LOGnOUT(5,<<"AveLoss/AveGain=\t"<<costMatrixGainLossRatio<<"\n");
		LOGnOUT(5,<<"Ave (loss/gain)=\t"<<init_losses2gainRatioForCostMatrixSum/numOfSequenceSets<<"\n");
		LOGnOUT(4,<<"All positions Q Ave="<<ratePerPosSum/(double)numOfSequenceSets<<"\n");
		time(&t2);
		LOGnOUT(4,<<"End simulations.\nTIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
		////////////////////////////////////////////////////////////////////////// end of per-position simulations

		//fastaFormat::write(cout,seqSimulated);
		//vector<int> posToRemove(seqSimulated.seqLen(),false);
		//posToRemove[0] = true;
		//seqSimulated.removePositions(posToRemove);
		//fastaFormat::write(cout,seqSimulated);

		//re-open seq
		string strSeqFirst = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(1) + ".fa";
		ifstream in(strSeqFirst.c_str());
		sequenceContainer seqReOpened = recognizeFormat::read(in,&alph);
		in.close();
		remove( strSeqFirst.c_str() ); // remove seq

		// Test for gain events in Eq sequences
		int totalNumberOfEqSeqs = 0;
		int totalNumberOfGainsInEqSeqs = 0;

		for(int i=1; i<numOfSequenceSets; i++){
			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
			ifstream in(strSeqNum.c_str());
			sequenceContainer seqSeqNum = recognizeFormat::read(in,&alph);
			in.close();
			seqReOpened.concatenate(seqSeqNum);
			if(isTestForGainEventsInEqSeq){
				if(_sc==seqSeqNum){
					++totalNumberOfEqSeqs;
					if(isGainEventInAnode[i+1]){			// to be consistent with previous calculations
						++totalNumberOfGainsInEqSeqs;
						//LOGnOUT(4,<<i+1<<" gain event\n");
					}
					//LOGnOUT(4,<<i+1<<" same seq\n");
				}
				else{
					//LOGnOUT(4,<<i<<" Diff seq\n");
				}
			}
			remove( strSeqNum.c_str() ); // remove seq
		}

		if(isTestForGainEventsInEqSeq){
			LOGnOUT(3,<<totalNumberOfEqSeqs<<" total same seqs\n");
			LOGnOUT(3,<<totalNumberOfGainsInEqSeqs<<" with gain event\n");
			LOGnOUT(3,<<(float)totalNumberOfGainsInEqSeqs/totalNumberOfEqSeqs<<" posteriorProb empirical\n");
		}
		LOGnOUT(5,<<"seqReOpened length "<<seqReOpened.seqLen()<<endl);
		string treeSimString = outDirSeq + "//" + "TreeSim.ph";
		string seqSim = outDirSeq + "//" + "seq" + ".fa";
		ofstream seq_out(seqSim.c_str());
		fastaFormat::  write(seq_out,seqReOpened);

		if(gainLossOptions::_isOnlySimulateSeq)
			continue;

		// Parsimony
		if(gainLossOptions::_calculeMaxParsimonyChangeSeveralGainLossRatios){			
			MDOUBLE GLratioMulti = 1;
			for(MDOUBLE glRatio = 1+glRatioTieBreakerInCostMatrix; glRatio <=MaxGLratio; glRatio+=GLratioMulti){
				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
					,glRatio,_distanceFromNearestOTUForRecent,false);
				if(glRatio>2)
					GLratioMulti*=2;
			}
			if(!gainLossOptions::_isMPratio && isMPcostEmpirical)
				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
				,costMatrixGainLossRatio*costMatrixGainLossRatioCorrectionFactor,_distanceFromNearestOTUForRecent,false);
			startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
				,loss2gainRatioToSim+glRatioTieBreakerInCostMatrix,_distanceFromNearestOTUForRecent,false);
		}
		else{
			if(isMPcostEmpirical)
				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
				,costMatrixGainLossRatio*costMatrixGainLossRatioCorrectionFactor,_distanceFromNearestOTUForRecent,false);
			startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
				,loss2gainRatioToSim+glRatioTieBreakerInCostMatrix,_distanceFromNearestOTUForRecent,false);
		}

		// Estimation of model paramers + Stochastic mapping
		tree trSim = _tr;
		if(gainLossOptions::_gainLossDist){
			cloneSpVVec(_spVVec,spVVecSim);
			gainDistSim = _gainDist->clone();
			lossDistSim = _lossDist->clone();
		}
		else{
			spSim = _sp->clone();
		}
		if(_unObservableData_p){
			unObservableDataSim = _unObservableData_p->clone();
		}
		if(gainLossOptions::_isFlatTreeBeforOpt){
			FlatTree(trSim);
		}		

		if(!gainLossOptions::_gainLossDist){// a single Stochastic processes (M)
			if(Parameters::getInt("_isFlatSpBeforeOpt")){
				FlatSpBeforeOpt(*spSim,unObservableDataSim);
			}
			if(Parameters::getInt("_isInitGainLossByEmpiricalFreqSimulatePostExp")){
				Vdouble freqSim = evaluateCharacterFreq(seqReOpened);
				LOGnOUT(4,<<"\nBefore optimization - init sp with simulated freq(1)= "<<freqSim[1]<<endl);
				MDOUBLE init_gain = freqSim[1];
				MDOUBLE init_loss = freqSim[0];
				static_cast<gainLossModel*>(spSim->getPijAccelerator()->getReplacementModel())->setMu1(init_gain, gainLossOptions::_isReversible);
				static_cast<gainLossModelNonReversible*>(spSim->getPijAccelerator()->getReplacementModel())->setMu2(init_loss);
				if(isThetaOptimization())
					static_cast<gainLossModel*>(spSim->getPijAccelerator()->getReplacementModel())->setTheta(freqSim[1]);
				printModellValuesOfParams(spSim,trSim);
				_logL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(trSim,seqReOpened,*spSim,_weightsUniqPatterns,unObservableDataSim);

				spSimpleSim =  startStochasticProcessSimpleGamma(freqSim[1],freqSim[0],freqSim); // simple initialization, based on empiricalCounting of '1' and '0'
				if(gainLossOptions::_isFlatTreeBeforOpt || gainLossOptions::_isbBLEMwithSimpleSpSimulatePostExp){
					bBLEMwithSimpleSpBeforeFullOptimization(trSim,seqReOpened,spSimpleSim,spSim,spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
				}
			}
			if(gainLossOptions::_modelOptimizationSimPostExp){
				gainLossOptimizer glOpt(trSim,spSim,seqReOpened,
					gainLossOptions::_epsilonOptimizationIterationCycle*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterations*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					gainLossOptions::_epsilonOptimizationModel*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterationsModel*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					gainLossOptions::_epsilonOptimizationBBL*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterationsBBL*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					NULL,unObservableDataSim, gainLossOptions::_BBLOptimizationSimPostExp, gainLossOptions::_isbblLSWhenbblEMdontImprove);
				if(gainLossOptions::_BBLOptimizationSimPostExp && printTreeForEachReplication){
					trSim = glOpt.getOptTree();
					printTree(trSim, treeSimString);
				}				
			}			
		}

		else{// Mixture of Stochastic processes (GLM)
			if(Parameters::getInt("_isFlatSpBeforeOpt")){
				FlatSpBeforeOpt(spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
			}
			if(Parameters::getInt("_isInitGainLossByEmpiricalFreqSimulatePostExp")){
				Vdouble freqSim = evaluateCharacterFreq(seqReOpened);
				LOGnOUT(4,<<"\nBefore optimization - init sp with simulated freq(1)= "<<freqSim[1]<<endl);
				MDOUBLE init_gain = freqSim[1];
				MDOUBLE init_loss = freqSim[0];
				MDOUBLE AlphasGainLossRatio = getRateAlpha(gainDistSim)/getRateAlpha(lossDistSim);
	
				updateGainBeta((1/init_gain)*(1/AlphasGainLossRatio), spVVecSim,gainDistSim,lossDistSim,false);
				updateLossBeta((1/init_loss)*AlphasGainLossRatio, spVVecSim,gainDistSim,lossDistSim,false);
				if(isThetaOptimization())
					updateTheta(freqSim[1], spVVecSim,gainDistSim, lossDistSim);
				normalizeQ(spVVecSim,gainDistSim,lossDistSim);
				printModellValuesOfParams(trSim,spVVecSim,gainDistSim,lossDistSim);
				_logL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(trSim,seqReOpened,spVVecSim,gainDistSim,lossDistSim,_weightsUniqPatterns,unObservableDataSim);

				spSimpleSim =  startStochasticProcessSimpleGamma(freqSim[1],freqSim[0],freqSim); // simple initialization, based on empiricalCounting of '1' and '0'
				if(gainLossOptions::_isFlatTreeBeforOpt || gainLossOptions::_isbBLEMwithSimpleSpSimulatePostExp){
					bBLEMwithSimpleSpBeforeFullOptimization(trSim,seqReOpened,spSimpleSim,spSim,spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
				}
			}
			if(gainLossOptions::_modelOptimizationSimPostExp){
				gainLossOptimizer glOpt(trSim,spVVecSim,gainDistSim,lossDistSim,seqReOpened,
					gainLossOptions::_epsilonOptimizationIterationCycle*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterations*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					gainLossOptions::_epsilonOptimizationModel*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterationsModel*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					gainLossOptions::_epsilonOptimizationBBL*gainLossOptions::_epsilonOptForPostExpSimFactor,
					(int)ceil(gainLossOptions::_maxNumOfIterationsBBL*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
					NULL, _unObservableData_p ,gainLossOptions::_BBLOptimizationSimPostExp, gainLossOptions::_isbblLSWhenbblEMdontImprove);
				if(gainLossOptions::_BBLOptimizationSimPostExp && printTreeForEachReplication){
					trSim = glOpt.getOptTree();
					printTree(trSim, treeSimString);
				}		
			}				
		}
		//////////////////////////////////////// compute Stochastic Mapping
		MDOUBLE distanceFromNearestOTUForRecent = computeDistanceNearestOTUforRecent(trSim);
		if(!gainLossOptions::_gainLossDist){
			startComputePosteriorExpectationOfChange(seqReOpened,trSim,spSim,LpostPerCatSim,unObservableDataSim,outDirSeq,distanceFromNearestOTUForRecent,false);
			LpostPerCatSim.clear();	// when cleared - each replicate will recompute the _LpostPerCat
		}
		else{
			startComputePosteriorExpectationOfChange(seqReOpened,trSim,spVVecSim,gainDistSim,lossDistSim,LpostPerSpPerCatSim,unObservableDataSim,outDirSeq,distanceFromNearestOTUForRecent,false);
			LpostPerSpPerCatSim.clear();	// when cleared - each replicate will recompute the _LpostPerSpPerCat
		}
		time(&t2);
		LOGnOUT(4,<<"Replicate SimultePosteriorExpectationOfChange RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	}	
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLoss::ComputeEmpiricalExpectedQforStationaryProcess(VVVdouble& EmpPerPos, MDOUBLE minRate){
	MDOUBLE expectedQvalEmpirical = 0;
	for(int pos=0; pos<EmpPerPos.size();++pos){
		MDOUBLE rateForPos=0;
		MDOUBLE gain =EmpPerPos[pos][0][1];
		MDOUBLE loss =EmpPerPos[pos][1][0];
		gain += minRate; 
		loss += minRate; 
		MDOUBLE Freq_0 = loss/(gain+loss);
		MDOUBLE Freq_1 = gain/(gain+loss);
		rateForPos = gain*Freq_0 + loss*Freq_1;
		expectedQvalEmpirical += rateForPos;
	}
	expectedQvalEmpirical /= EmpPerPos.size();
	return expectedQvalEmpirical;
}



/********************************************************************************************
*********************************************************************************************/
void gainLoss::RemoveSeqWithUnknownForSelectedSiteForCorrelation(sequenceContainer&  sc, tree& tr){

	gainLossAlphabet alph;
	int pos_2_remove = _sc.seqLen()-1;
	char char_2_match = alph.unknown();
	vector<int> seqIDs2remove;
	vector<string> SeqNamesThatMatchPos = _sc.getSeqNamesThatMatchPos(pos_2_remove,char_2_match);

	sequenceContainer::constTaxaIterator myseq=sc.constTaxaBegin();
	for (;myseq != sc.constTaxaEnd(); ++myseq){
		bool bFound = false;
		for (int i=0; i<SeqNamesThatMatchPos.size(); ++i) {			
	
			if (myseq->name() == SeqNamesThatMatchPos[i]) 
			{
				bFound = true;
				break;
			}
		}
		if (bFound == true) 
		{
			string errMsg = "The taxID name:\t";
			errMsg += myseq->name();
			errMsg += "\twas found in with missing data. Removed.";
			LOGnOUT(4,<<errMsg<<endl);
			seqIDs2remove.push_back(myseq->id());			
		}
	}
	for(int i=0; i<seqIDs2remove.size(); ++i){
		sc.remove(seqIDs2remove[i]);
	}
	intersectNamesInTreeAndSequenceContainer(tr,sc);
	// Write seq and tree (required for re-labeling IDs
	string strSeqNum = gainLossOptions::_outDir + "//" + "seq.noUnknown.fa";
	ofstream seq_out(strSeqNum.c_str());
	fastaFormat::  write(seq_out,sc);
	string treeSampled = gainLossOptions::_outDir + "//" + "TheTree.noUnknonwn.ph"; 
	ofstream treeStream(treeSampled.c_str());
	tr.output(treeStream);		

	// re-Read
	ifstream in(strSeqNum.c_str());
	sc = recognizeFormat::read(in,&alph);
	tr= tree(treeSampled);

}





/********************************************************************************************
*********************************************************************************************/
//void gainLoss::simultePhyleticData(const int numOfSequenceSets, string strSeqFirst,MDOUBLE loss2gainRatioToSim, gainLossOptions::simulationType simulationType
//								   , MDOUBLE AlphaGain, MDOUBLE BetaGain, MDOUBLE AlphaLoss, MDOUBLE BetaLoss, MDOUBLE AlphaRate)
//{
//	MDOUBLE minThetaRandSample = 0.1;	// was 0.01. change all from 0.01 to 0.05, and later to 0.1
//	MDOUBLE maxThetaRandSample = 0.9;	// was 0.09
//	MDOUBLE observedTheta = 0.5;
//	MDOUBLE minGainRandSample = 0.1;	// was 0.01
//	MDOUBLE maxGainRandSample = 2.0;	// was 2.5, now E(val) = 1	
//	MDOUBLE minLossRandSample = 0.1;	// was 0.01
//	MDOUBLE maxLossRandSample = loss2gainRatioToSim*2;
//	MDOUBLE meanGaussianGain = 1.0;
//	MDOUBLE varianceGaussianGain = 1.0;
//	MDOUBLE meanGaussianLoss = loss2gainRatioToSim;
//	MDOUBLE varianceGaussianLoss = loss2gainRatioToSim;
//	//MDOUBLE AlphaGain = gainLossOptions::_userAlphaGain;
//	//MDOUBLE BetaGain = gainLossOptions::_userBetaGain;
//	//MDOUBLE AlphaLoss = gainLossOptions::_userAlphaLoss;
//	//MDOUBLE BetaLoss = gainLossOptions::_userBetaLoss;
//	//MDOUBLE AlphaRate = gainLossOptions::_userAlphaRate;
//
//	Vdouble freq(2,0.0);
//	MDOUBLE init_gain = 1.0; //gainLossOptions::_userGain taken from original runs of COG data, get it from params file
//	MDOUBLE init_loss = 1.0; //gainLossOptions::_userLoss
//
//	
//	MDOUBLE init_gainsForCostMatrix = 0.0; // sum all position
//	MDOUBLE init_lossesForCostMatrix = 0.0; // sum all position
//	MDOUBLE init_losses2gainRatioForCostMatrixSum = 0.0;
//	MDOUBLE randomNoise = 0.0;
//	MDOUBLE costMatrixGainLossRatio = gainLossOptions::_costMatrixGainLossRatio; // to be updated according to simulation
//
//	// produce random noise
//	if(gainLossOptions::_simulationType == gainLossOptions::GammaNoise){
//		randomNoise = talRandom::giveRandomNumberBetweenTwoPoints(-gainLossOptions::_noiseLevelInGammaSimulation, gainLossOptions::_noiseLevelInGammaSimulation);
//		// if noiseLevel=200% than param may be up to x3 or down to x0.33 its value
//		if(randomNoise>=0)
//			randomNoise = 1+randomNoise;
//		else
//			randomNoise = 1/(1-randomNoise);
//		LOGnOUT(4,<<"Noise over all parameters="<< randomNoise<<endl);
//	}		
//	// Theta for all positions, (not relevant to stationary models)
//	if(!gainLossOptions::_isStationaryModelForSim){ // else Theta is driven from gain/gain+loss
//		switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
//		{
//		case gainLossOptions::Uniform:
//		case gainLossOptions::Normal:
//			freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			break;
//		case gainLossOptions::MPestEmp:
//		case gainLossOptions::SMestEmp:
//			if(isThetaFromObservedForEmpiricalSimulations)
//				freq[1]=observedTheta;
//			else
//				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			break;
//		case gainLossOptions::Gamma:
//			if(isThetaSampledForGamma)
//				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			else
//				freq[1]= gainLossOptions::_userTheta;
//			break;
//		case gainLossOptions::GammaNoise:
//			freq[1] = gainLossOptions::_userTheta*randomNoise;
//			freq[1] = max(freq[1],minThetaRandSample);	// added to avoid too small or too big theta
//			freq[1] = min(freq[1],maxThetaRandSample);
//			break;
//		case gainLossOptions::EQ_gEql:
//		case gainLossOptions::EQ_gVrl:
//		case gainLossOptions::Gam_gEql:
//		case gainLossOptions::Gam_gVrl:
//			freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			break;
//		default:
//			errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise}");
//		}
//		freq[0]= 1 - freq[1];
//		if(!gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos)
//			LOGnOUT(4,<<" For all positions, Root(1)= "<<freq[1]<<endl);
//	}
//	else{
//		LOGnOUT(4,<<" Statrionary model - Theta=Root(1) is driven from the gain/gain+loss"<<endl);
//	}
//
//	vector<bool> isGainEventInAnode;	// DEBUG
//	if(isTestForGainEventsInEqSeq)
//		isGainEventInAnode.resize(numOfSequenceSets+1);
//
//	////////////////////////////////////////////////////////////////////////// Positions
//	int randomPosition; //used in MPestEmp or SMestEmp
//	for(int i=0; i<numOfSequenceSets; ++i)
//	{
//		rateSample = 1;
//		lossGainRatioSample = 1;
//		switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
//		{
//			case gainLossOptions::Uniform:
//				init_gain = talRandom::giveRandomNumberBetweenTwoPoints(minGainRandSample, maxGainRandSample);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = talRandom::giveRandomNumberBetweenTwoPoints(minLossRandSample, maxLossRandSample);
//				break;
//			case gainLossOptions::Normal:
//				init_gain = talRandom::rand_gaussian(meanGaussianGain, varianceGaussianGain);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = talRandom::rand_gaussian(meanGaussianLoss, varianceGaussianLoss);					
//				break;
//			case gainLossOptions::Gamma:
//				init_gain = talRandom::SampleGamma(AlphaGain,BetaGain);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else					
//					init_loss = talRandom::SampleGamma(AlphaLoss,BetaLoss);					
//				break;
//			case gainLossOptions::MPestEmp:
//				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _MPPerPos.size());
//				init_gain = _MPPerPos[randomPosition][0][1]/meanEventsFromEMP;
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = _MPPerPos[randomPosition][1][0]/meanEventsFromEMP;
//				break;
//			case gainLossOptions::SMestEmp:
//				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _SMPerPos.size());
//				init_gain = _SMPerPos[randomPosition][0][1]/meanEventsFromEMP;
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = _SMPerPos[randomPosition][1][0]/meanEventsFromEMP;
//				break;
//			case gainLossOptions::GammaNoise:
//				init_gain = talRandom::SampleGamma( (AlphaGain*randomNoise)
//					,(BetaGain*randomNoise));
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else{					
//					init_loss = talRandom::SampleGamma((AlphaLoss*randomNoise)
//						,(BetaLoss*randomNoise));
//				}
//				break;
//			case gainLossOptions::EQ_gEql:
//				init_gain = -(rateSample/(-1-lossGainRatioSample));			//init_gain = init_loss =0.5;
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::Gam_gEql:
//				rateSample = talRandom::SampleGamma(AlphaRate);	//init_gain = init_loss = 0.5*rateSample;
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::EQ_gVrl:
//				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::Gam_gVrl:
//				rateSample = talRandom::SampleGamma(AlphaRate);
//				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				//init_gain = -(1/(-1-lossGainRatioSample))*rateSample;
//				//init_loss = 1+(1/(-1-lossGainRatioSample))*rateSample;
//				break;				
//			default:
//				errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp GammaNoise}");
//		}
//		//cout<<init_loss/init_gain<<"\n";
//
//		init_gain = max(init_gain,minAllowedRate);	// added to avoid too small gain rate
//		init_gain = min(init_gain,maxAllowedRate);
//		init_loss = max(init_loss,minAllowedRate);
//		init_loss = min(init_loss,maxAllowedRate);
//
//		init_losses2gainRatioForCostMatrixSum += init_loss/init_gain;
//		init_gainsForCostMatrix += init_gain;
//		init_lossesForCostMatrix += init_loss;			
//
//		///////////// Theta random per pos 
//		if(gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos){	
//			freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			freq[0]= 1 - freq[1];
//		}			
//		if(gainLossOptions::_isStationaryModelForSim){	//Theta=Root(1) is driven from the gain/gain+loss
//			freq[1]= init_gain/(init_gain+init_loss);
//			freq[0]= 1 - freq[1];
//		}
//		gainLossModelNonReversible  glm(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,_isHGT_normal_Pij,_isHGT_with_Q);
//		trivialAccelerator pijAcc(&glm);
//		uniDistribution uniDistr;
//		stochasticProcess *spSimSingle;
//		spSimSingle = new stochasticProcess(&uniDistr,&pijAcc,false);
//		if(isnormalizeQ){
//			MDOUBLE sumQii = normalizeQ(spSimSingle);	// added.
//			LOG(6,<<" Pos= "<<i+1<<
//				"\tfreq1="<<freq[1]<<"\tgain="<<init_gain/sumQii<<"\tloss="<<init_loss/sumQii<<endl);
//		}
//		string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
//		string resFile = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "resSim" + int2string(i+1) + ".sim";		
//		simulateOnePos simulateOnePosObj(strSeqNum, resFile, i,gainLossOptions::_treeFile,spSimSingle);			
//		if(spSimSingle) delete spSimSingle;
//		if(isTestForGainEventsInEqSeq){	// DEBUG
//			if(simulateOnePosObj.getChangesForBranch(2)[0][1]>0)	// "A" == 2
//				isGainEventInAnode[i+1] = true;
//		}			
//	}
//	if(gainLossOptions::_isMatrixGainLossFromRatioInSimulations) // e.g., val=2, loss rate is double that of loss
//		costMatrixGainLossRatio = init_lossesForCostMatrix/init_gainsForCostMatrix;						
//
//	cout<<"AveLoss/AveGain"<<costMatrixGainLossRatio<<"\n";
//	cout<<"Ave(loss/gain)"<<init_losses2gainRatioForCostMatrixSum/numOfSequenceSets<<"\n";		
//	////////////////////////////////////////////////////////////////////////// end of per-position simulations
//
//	//re-open seq
//	gainLossAlphabet alph;
//	//string strSeqFirst = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(1) + ".fa";
//	ifstream in(strSeqFirst.c_str());
//	sequenceContainer seqReOpened = recognizeFormat::read(in,&alph);
//	in.close();
//	remove( strSeqFirst.c_str() ); // remove seq
//
//	// Test for gain events in Eq sequences
//	int totalNumberOfEqSeqs = 0;
//	int totalNumberOfGainsInEqSeqs = 0;
//
//	for(int i=1; i<numOfSequenceSets; i++){
//		string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
//		ifstream in(strSeqNum.c_str());
//		sequenceContainer seqSeqNum = recognizeFormat::read(in,&alph);
//		in.close();
//		seqReOpened.concatenate(seqSeqNum);
//		remove( strSeqNum.c_str() ); // remove seq
//	}
//}


//void gainLoss::startSimultePosteriorExpectationOfChange(int numOfSequenceSets, const int numOfRepeats)
//{
//	LOGnOUT(4,<<endl<<"****************************************************\n startSimultePosteriorExpectationOfChange... "<<endl);
//	LOGnOUT(4,<<"Replicates="<<numOfRepeats<<" Positions="<<numOfSequenceSets<<endl);
//	LOGnOUT(4,<<" simulationType {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise;..."<<endl);
//	LOGnOUT(4,<<" EQ_gEql,EQ_gVrl,Gam_gEql,GamgVrl}="<<gainLossOptions::_simulationType<<endl);
//	
//	if(gainLossOptions::_simulationType == gainLossOptions::SMestEmp && _SMPerPos.size()==0){
//		LOGnOUT(4,<<" WARN!!! _SMPerPos size="<<_SMPerPos.size()<<endl);
//		startComputePosteriorExpectationOfChange();
//	}
//	
//	if(gainLossOptions::_simulationType == gainLossOptions::MPestEmp && _MPPerPos.size()==0){
//		LOGnOUT(4,<<" WARN!!! _MPPerPos size="<<_MPPerPos.size()<<endl);
//		startMaxParsimonyChange();
//	}
//	gainLossAlphabet alph;
//	simulatePhyleticPatternsAndPredictEvents simulateObj(_tr,_sp, alph);
//	
//
//
//	////////////////////////////////////////////////////////////////////////// Replicates	
//	for(int replicat=1; replicat<=numOfRepeats; ++replicat){
//		LOGnOUT(4,<<endl<<".......................................Replicate= "<<replicat<<endl);
//		time_t t1,t2;
//		time(&t1);
//
//		createDir(gainLossOptions::_outDir, "SimulatedPostExp"+ int2string(replicat));
//		string outDirSeq = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seqAll" ;
//		createDir(gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat), "seqAll" );
//		string simulatedEventsSimString = outDirSeq + "//" + "simulatedEvents.txt";
//		ofstream* simulatedEventsFile = new ofstream(simulatedEventsSimString.c_str());
//
//		string perPosStat = outDirSeq + "//" + "statPos.txt";
//		ofstream perPosStatStream(perPosStat.c_str());
//		perPosStatStream<<"pos"<<"\t"<<"rate"<<"\t"<<"theta"<<"\t"<<"occur"<<"\n";
//
//		string perBranchStat = outDirSeq + "//" + "statBranch.txt";
//		ofstream perBranchStatStream(perBranchStat.c_str());
//		perBranchStatStream<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<endl;
//		treeIterTopDownConst tit(_tr);
//		for (tree::nodeP myN = tit.first();myN!=tit.end(); myN = tit.next()) {
//			if(myN->isRoot())
//				continue;
//			perBranchStatStream<<myN->name()<<"\t"<<myN->dis2father()<<"\t"<<myN->getDistance2ROOT()<<"\t"<<myN->getMinimalDistance2OTU()<<"\t"<<myN->getMinimalNumOfNodes2OTU()<<endl;
//		}
//		perBranchStatStream.close();
//		MDOUBLE init_gainsForCostMatrix = 0.0; // sum all position
//		MDOUBLE init_lossesForCostMatrix = 0.0; // sum all position
//		MDOUBLE init_losses2gainRatioForCostMatrixSum = 0.0;
//		MDOUBLE QnormTest = 0.0;
//
//		// produce random noise
//		if(gainLossOptions::_simulationType == gainLossOptions::GammaNoise){
//			randomNoise = talRandom::giveRandomNumberBetweenTwoPoints(-gainLossOptions::_noiseLevelInGammaSimulation, gainLossOptions::_noiseLevelInGammaSimulation);
//			// if noiseLevel=200% than param may be up to x3 or down to x0.33 its value
//			if(randomNoise>=0)
//				randomNoise = 1+randomNoise;
//			else
//				randomNoise = 1/(1-randomNoise);
//			LOGnOUT(4,<<"Noise over all parameters="<< randomNoise<<endl);
//		}		
//		// Theta for all positions, (not relevant to stationary models)
//		if(!gainLossOptions::_isStationaryModelForSim){ // else Theta is driven from gain/gain+loss
//			switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
//			{
//			case gainLossOptions::Uniform:
//			case gainLossOptions::Normal:
//				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//				break;
//			case gainLossOptions::MPestEmp:
//			case gainLossOptions::SMestEmp:
//				if(isThetaFromObservedForEmpiricalSimulations)
//					freq[1]=observedTheta;
//				else
//					freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//				break;
//			case gainLossOptions::Gamma:
//				if(isThetaSampledForGamma)
//					freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//				else
//					freq[1]= gainLossOptions::_userTheta;
//				break;
//			case gainLossOptions::GammaNoise:
//				freq[1] = gainLossOptions::_userTheta*randomNoise;
//				freq[1] = max(freq[1],minThetaRandSample);	// added to avoid too small or too big theta
//				freq[1] = min(freq[1],maxThetaRandSample);
//				break;
//			case gainLossOptions::EQ_gEql:
//			case gainLossOptions::EQ_gVrl:
//			case gainLossOptions::Gam_gEql:
//			case gainLossOptions::Gam_gVrl:
//				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//				break;
//			default:
//				errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp,SMestEmp, GammaNoise}");
//			}
//			if(!gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos)
//				LOGnOUT(4,<<" For all positions, Root(1)= "<<freq[1]<<endl);
//		}
//		else{
//			LOGnOUT(4,<<" Stationary model - Theta=Root(1) is driven from the gain/gain+loss"<<endl);
//		}
//
//		vector<bool> isGainEventInAnode;	// DEBUG
//		if(isTestForGainEventsInEqSeq)
//			isGainEventInAnode.resize(numOfSequenceSets+1);
//
//		////////////////////////////////////////////////////////////////////////// Positions
//		int randomPosition; //used in MPestEmp or SMestEmp
//		sequenceContainer seqSimulated;
//		gainLossAlphabet alph;
//
//
//		for(int i=0; i<numOfSequenceSets; ++i)
//		{
//			rateSample = 1;
//			lossGainRatioSample = 1;
//			switch (gainLossOptions::_simulationType) //{Uniform, Normal, Gamma, MPestEmp, GammaNoise}
//			{
//			case gainLossOptions::Uniform:
//				init_gain = talRandom::giveRandomNumberBetweenTwoPoints(minGainRandSample, maxGainRandSample);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = talRandom::giveRandomNumberBetweenTwoPoints(minLossRandSample, maxLossRandSample);
//				break;
//			case gainLossOptions::Normal:
//				init_gain = talRandom::rand_gaussian(meanGaussianGain, varianceGaussianGain);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = talRandom::rand_gaussian(meanGaussianLoss, varianceGaussianLoss);					
//				break;
//			case gainLossOptions::Gamma:
//				init_gain = talRandom::SampleGamma(AlphaGain,BetaGain);
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else					
//					init_loss = talRandom::SampleGamma(AlphaLoss,BetaLoss);					
//				break;
//			case gainLossOptions::MPestEmp:
//				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _MPPerPos.size());
//				init_gain = _MPPerPos[randomPosition][0][1]/meanEventsFromEMP;
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = _MPPerPos[randomPosition][1][0]/meanEventsFromEMP;
//				break;
//			case gainLossOptions::SMestEmp:
//				randomPosition = (int)talRandom::giveRandomNumberBetweenTwoPoints(0, _SMPerPos.size());
//				init_gain = _SMPerPos[randomPosition][0][1]/meanEventsFromEMP;
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else
//					init_loss = _SMPerPos[randomPosition][1][0]/meanEventsFromEMP;
//				break;
//			case gainLossOptions::GammaNoise:
//				init_gain = talRandom::SampleGamma( (AlphaGain*randomNoise)
//					,(BetaGain*randomNoise));
//				if(gainLossOptions::_isMPratio)
//					init_loss = init_gain * costMatrixGainLossRatio;
//				else{					
//					init_loss = talRandom::SampleGamma((AlphaLoss*randomNoise)
//						,(BetaLoss*randomNoise));
//				}
//				break;
//			case gainLossOptions::EQ_gEql:
//				init_gain = -(rateSample/(-1-lossGainRatioSample));			//init_gain = init_loss =0.5;
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::Gam_gEql:
//				rateSample = talRandom::SampleGamma(AlphaRate);	//init_gain = init_loss = 0.5*rateSample;
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::EQ_gVrl:
//				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				break;
//			case gainLossOptions::Gam_gVrl:
//				rateSample = talRandom::SampleGamma(AlphaRate);
//				lossGainRatioSample = talRandom::giveRandomNumberBetweenTwoPoints(epsilonForgainLossRatio, loss2gainRatioToSim*2-epsilonForgainLossRatio);
//				init_gain = -(rateSample/(-1-lossGainRatioSample));
//				init_loss = rateSample+(rateSample/(-1-lossGainRatioSample));
//				//init_gain = -(1/(-1-lossGainRatioSample))*rateSample;
//				//init_loss = 1+(1/(-1-lossGainRatioSample))*rateSample;
//				break;				
//			default:
//				errorMsg::reportError("unknown type in optimizationLevel - {Uniform, Normal, Gamma, MPestEmp GammaNoise}");
//			}
//			if(isMultBy2_normQ){	// with mult=2, Q matrix is normalized with respect to the tree (after multiplied by freq=0.5)
//				init_gain *=2;
//				init_loss *=2;
//			}
//			init_gain = max(init_gain,minAllowedRate);	// added to avoid too small gain rate
//			init_gain = min(init_gain,maxAllowedRate);
//			init_loss = max(init_loss,minAllowedRate);
//			init_loss = min(init_loss,maxAllowedRate);
//
//			///////////// Theta random per pos 
//			if(gainLossOptions::_initRootFreqAtRandPointsInSimPostExpEachPos){	
//				freq[1]= talRandom::giveRandomNumberBetweenTwoPoints(minThetaRandSample, maxThetaRandSample);
//			}			
//			if(gainLossOptions::_isStationaryModelForSim){	//Theta=Root(1) is driven from the gain/gain+loss
//				freq[1]= init_gain/(init_gain+init_loss);
//			}
//			freq[0]= 1 - freq[1];
//			gainLossModelNonReversible  glm(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,_isHGT_normal_Pij,_isHGT_with_Q);
//			trivialAccelerator pijAcc(&glm);
//			uniDistribution uniDistr;
//			stochasticProcess *spSimSingle;
//			spSimSingle = new stochasticProcess(&uniDistr,&pijAcc,false);
//			MDOUBLE sumQii = 1.0;
//			if(isNormalizeQ){	// if normalizeQ for each position, there is no rate variability in practice
//				sumQii = normalizeQ(spSimSingle);					
//				LOG(6,<<" Pos= "<<i+1<<
//					"\tfreq1="<<freq[1]<<"\tgain="<<init_gain/sumQii<<"\tloss="<<init_loss/sumQii<<endl);
//			}			
//			//QnormTest +=  init_gain*freq[0]+init_loss*freq[1];
//			init_losses2gainRatioForCostMatrixSum += init_loss/init_gain;
//			init_gainsForCostMatrix += init_gain;
//			init_lossesForCostMatrix += init_loss;
//
//			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
//			string resFile = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "resSim" + int2string(i+1) + ".sim";
//
//			simulateOnePos *simulateOnePosObj = NULL;
//			MDOUBLE ratePerPos=0;
//			if(gainLossOptions::_is3states){
//				Vdouble init_cpN_vals(4);
//				init_cpN_vals[0]=gainLossOptions::_3statesGain; //gain (0->1)
//				init_cpN_vals[1]=gainLossOptions::_3statesMore; //more (1->more)
//				init_cpN_vals[2]=gainLossOptions::_3statesLess; // less (more->1) 
//				init_cpN_vals[3]=gainLossOptions::_3statesLoss; // loss (1->0)
//				Vdouble freq_cpN(3);
//				freq_cpN[0]=gainLossOptions::_3states0;
//				freq_cpN[1]=gainLossOptions::_3states1;
//				freq_cpN[2]=1 - (freq_cpN[0] + freq_cpN[1]);
//				simulateOnePosObj = new simulateOnePos(strSeqNum, resFile, simulatedEventsFile, i,gainLossOptions::_treeFile,init_cpN_vals[0]+init_cpN_vals[3],freq[1],gainLossOptions::_is3states,NULL,&_tr,&init_cpN_vals,&freq_cpN);
//			}
//			else{
//				ratePerPos=(static_cast<gainLossModel*>(spSimSingle->getPijAccelerator()->getReplacementModel()))->sumPijQij();
//				simulateOnePosObj = new simulateOnePos(strSeqNum, resFile, simulatedEventsFile, i,gainLossOptions::_treeFile,ratePerPos,freq[1],gainLossOptions::_is3states,spSimSingle,&_tr);
//			}
//			perPosStatStream<<i+1<<"\t"<<ratePerPos<<"\t"<<freq[1]<<"\t"<<simulateOnePosObj->getOccurFraction()<<"\n";
//
//
//			if(spSimSingle) delete spSimSingle;
//			if(isTestForGainEventsInEqSeq){	// DEBUG
//				if(simulateOnePosObj->getChangesForBranch(2)[0][1]>0)	// "A" == 2
//					isGainEventInAnode[i+1] = true;
//			}			
//			//if(i==0){	
//			//	seqSimulated = sequenceContainer(simulateOnePosObj->getSequenceContainer(),&alph);
//			//}
//			//else{
//			//	sequenceContainer tempSeq = sequenceContainer(simulateOnePosObj->getSequenceContainer(),&alph);
//			//	seqSimulated.concatenate(tempSeq);
//			//	fastaFormat::write(cout,seqSimulated);
//			//}
//
//		}
//		if(gainLossOptions::_isMatrixGainLossFromRatioInSimulations) // e.g., val=2, loss rate is double that of loss
//			costMatrixGainLossRatio = init_lossesForCostMatrix/init_gainsForCostMatrix;						
//
//		//LOGnOUT(5,<<"QnormTest=\t"<<QnormTest/numOfSequenceSets<<"\n");		
//		LOGnOUT(5,<<"AveLoss/AveGain=\t"<<costMatrixGainLossRatio<<"\n");
//		LOGnOUT(5,<<"Ave (loss/gain)=\t"<<init_losses2gainRatioForCostMatrixSum/numOfSequenceSets<<"\n");
//
//		time(&t2);
//		LOGnOUT(4,<<"End simulations.\nTIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
//		////////////////////////////////////////////////////////////////////////// end of per-position simulations
//
//		//fastaFormat::write(cout,seqSimulated);
//		//vector<int> posToRemove(seqSimulated.seqLen(),false);
//		//posToRemove[0] = true;
//		//seqSimulated.removePositions(posToRemove);
//		//fastaFormat::write(cout,seqSimulated);
//
//		//re-open seq
//		string strSeqFirst = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(1) + ".fa";
//		ifstream in(strSeqFirst.c_str());
//		sequenceContainer seqReOpened = recognizeFormat::read(in,&alph);
//		in.close();
//		remove( strSeqFirst.c_str() ); // remove seq
//
//		// Test for gain events in Eq sequences
//		int totalNumberOfEqSeqs = 0;
//		int totalNumberOfGainsInEqSeqs = 0;
//
//		for(int i=1; i<numOfSequenceSets; i++){
//			string strSeqNum = gainLossOptions::_outDir + "//" + "SimulatedPostExp"+ int2string(replicat) + "//" + "seq" + int2string(i+1) + ".fa";
//			ifstream in(strSeqNum.c_str());
//			sequenceContainer seqSeqNum = recognizeFormat::read(in,&alph);
//			in.close();
//			seqReOpened.concatenate(seqSeqNum);
//			if(isTestForGainEventsInEqSeq){
//				if(_sc==seqSeqNum){
//					++totalNumberOfEqSeqs;
//					if(isGainEventInAnode[i+1]){			// to be consistent with previous calculations
//						++totalNumberOfGainsInEqSeqs;
//						//LOGnOUT(4,<<i+1<<" gain event\n");
//					}
//					//LOGnOUT(4,<<i+1<<" same seq\n");
//				}
//				else{
//					//LOGnOUT(4,<<i<<" Diff seq\n");
//				}
//			}
//			remove( strSeqNum.c_str() ); // remove seq
//		}
//
//		if(isTestForGainEventsInEqSeq){
//			LOGnOUT(3,<<totalNumberOfEqSeqs<<" total same seqs\n");
//			LOGnOUT(3,<<totalNumberOfGainsInEqSeqs<<" with gain event\n");
//			LOGnOUT(3,<<(float)totalNumberOfGainsInEqSeqs/totalNumberOfEqSeqs<<" posteriorProb empirical\n");
//		}
//		LOGnOUT(5,<<"seqReOpened length "<<seqReOpened.seqLen()<<endl);
//		string treeSimString = outDirSeq + "//" + "TreeSim.ph";
//		string seqSim = outDirSeq + "//" + "seq" + ".fa";
//		ofstream seq_out(seqSim.c_str());
//		fastaFormat::  write(seq_out,seqReOpened);
//
//		// Parsimony
//		if(gainLossOptions::_calculeMaxParsimonyChangeSeveralGainLossRatios){			
//			MDOUBLE GLratioMulti = 1;
//			for(MDOUBLE glRatio = 1+glRatioTieBreakerInCostMatrix; glRatio <=MaxGLratio; glRatio+=GLratioMulti){
//				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
//					,glRatio,_distanceFromNearestOTUForRecent);
//				if(glRatio>2)
//					GLratioMulti*=2;
//			}
//			if(!gainLossOptions::_isMPratio && isMPcostEmpirical)
//				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
//				,costMatrixGainLossRatio*costMatrixGainLossRatioCorrectionFactor,_distanceFromNearestOTUForRecent);
//			startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
//				,loss2gainRatioToSim+glRatioTieBreakerInCostMatrix,_distanceFromNearestOTUForRecent);
//		}
//		else{
//			if(isMPcostEmpirical)
//				startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
//				,costMatrixGainLossRatio*costMatrixGainLossRatioCorrectionFactor,_distanceFromNearestOTUForRecent);
//			startMaxParsimonyChange(seqReOpened,_tr,outDirSeq
//				,loss2gainRatioToSim+glRatioTieBreakerInCostMatrix,_distanceFromNearestOTUForRecent);
//		}
//
//		// Estimation of model paramers + Stochastic mapping
//		tree trSim = _tr;
//		if(gainLossOptions::_gainLossDist){
//			cloneSpVVec(_spVVec,spVVecSim);
//			gainDistSim = _gainDist->clone();
//			lossDistSim = _lossDist->clone();
//		}
//		else{
//			spSim = _sp->clone();
//		}
//		if(_unObservableData_p){
//			unObservableDataSim = _unObservableData_p->clone();
//		}
//		if(gainLossOptions::_isFlatTreeBeforOpt){
//			FlatTree(trSim);
//		}		
//
//		if(!gainLossOptions::_gainLossDist){// a single Stochastic processes (M)
//			if(Parameters::getInt("_isFlatSpBeforeOpt")){
//				FlatSpBeforeOpt(*spSim,unObservableDataSim);
//			}
//			if(Parameters::getInt("_isInitGainLossByEmpiricalFreqSimulatePostExp")){
//				Vdouble freqSim = evaluateCharacterFreq(seqReOpened);
//				LOGnOUT(4,<<"\nBefore optimization - init sp with simulated freq(1)= "<<freqSim[1]<<endl);
//				MDOUBLE init_gain = freqSim[1];
//				MDOUBLE init_loss = freqSim[0];
//				static_cast<gainLossModel*>(spSim->getPijAccelerator()->getReplacementModel())->setMu1(init_gain, gainLossOptions::_isReversible);
//				static_cast<gainLossModelNonReversible*>(spSim->getPijAccelerator()->getReplacementModel())->setMu2(init_loss);
//				static_cast<gainLossModel*>(spSim->getPijAccelerator()->getReplacementModel())->setTheta(freqSim[1]);
//
//				spSimpleSim =  startStochasticProcessSimpleGamma(freqSim[1],freqSim[0],freqSim); // simple initialization, based on empiricalCounting of '1' and '0'
//				_logL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(trSim,seqReOpened,*spSim,_weightsUniqPatterns,unObservableDataSim);
//
//				if(gainLossOptions::_isFlatTreeBeforOpt || gainLossOptions::_isbBLEMwithSimpleSpSimulatePostExp){
//					bBLEMwithSimpleSpBeforeFullOptimization(trSim,seqReOpened,spSimpleSim,spSim,spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
//				}
//			}
//			if(gainLossOptions::_modelOptimizationSimPostExp){
//				gainLossOptimizer glOpt(trSim,spSim,seqReOpened,
//					gainLossOptions::_epsilonOptimizationIterationCycle*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterations*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					gainLossOptions::_epsilonOptimizationModel*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterationsModel*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					gainLossOptions::_epsilonOptimizationBBL*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterationsBBL*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					NULL,unObservableDataSim, gainLossOptions::_BBLOptimizationSimPostExp, gainLossOptions::_isbblLSWhenbblEMdontImprove);
//				if(gainLossOptions::_BBLOptimizationSimPostExp && printTreeForEachReplication){
//					trSim = glOpt.getOptTree();
//					printTree(trSim, treeSimString);
//				}				
//			}			
//		}
//
//		else{// Mixture of Stochastic processes (GLM)
//			if(Parameters::getInt("_isFlatSpBeforeOpt")){
//				FlatSpBeforeOpt(spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
//			}
//			if(Parameters::getInt("_isInitGainLossByEmpiricalFreqSimulatePostExp")){
//				Vdouble freqSim = evaluateCharacterFreq(seqReOpened);
//				LOGnOUT(4,<<"\nBefore optimization - init ssp with simulated freq(1)= "<<freqSim[1]<<endl);
//				MDOUBLE init_gain = freqSim[1];
//				MDOUBLE init_loss = freqSim[0];
//				MDOUBLE gainLossRatioToCompleteByBeta = (init_gain/init_loss)*(gainLossOptions::_userAlphaLoss/gainLossOptions::_userAlphaGain);
//				MDOUBLE initBetaGain =sqrt(1/gainLossRatioToCompleteByBeta);			// AlphaGain = 0.35
//				MDOUBLE initBetaLoss =sqrt(gainLossRatioToCompleteByBeta);				// AlphaLoss = 0.9
//				updateGainBeta(initBetaGain, spVVecSim,gainDistSim,lossDistSim);
//				updateLossBeta(initBetaLoss, spVVecSim,gainDistSim,lossDistSim);
//				updateTheta(freqSim[1], spVVecSim,gainDistSim, lossDistSim);
//
//				spSimpleSim =  startStochasticProcessSimpleGamma(freqSim[1],freqSim[0],freqSim); // simple initialization, based on empiricalCounting of '1' and '0'
//				_logL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(trSim,seqReOpened,spVVecSim,gainDistSim,lossDistSim,_weightsUniqPatterns,unObservableDataSim);
//
//				if(gainLossOptions::_isFlatTreeBeforOpt || gainLossOptions::_isbBLEMwithSimpleSpSimulatePostExp){
//					bBLEMwithSimpleSpBeforeFullOptimization(trSim,seqReOpened,spSimpleSim,spSim,spVVecSim,gainDistSim,lossDistSim,unObservableDataSim);
//				}
//			}
//			if(gainLossOptions::_modelOptimizationSimPostExp){
//				gainLossOptimizer glOpt(trSim,spVVecSim,gainDistSim,lossDistSim,seqReOpened,_spSimple,
//					gainLossOptions::_epsilonOptimizationIterationCycle*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterations*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					gainLossOptions::_epsilonOptimizationModel*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterationsModel*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					gainLossOptions::_epsilonOptimizationBBL*gainLossOptions::_epsilonOptForPostExpSimFactor,
//					(int)ceil(gainLossOptions::_maxNumOfIterationsBBL*gainLossOptions::_numOfIterationsOptForPostExpSimFactor),
//					NULL, _unObservableData_p ,gainLossOptions::_BBLOptimizationSimPostExp, gainLossOptions::_isbblLSWhenbblEMdontImprove);
//				if(gainLossOptions::_BBLOptimizationSimPostExp && printTreeForEachReplication){
//					trSim = glOpt.getOptTree();
//					printTree(trSim, treeSimString);
//				}		
//			}				
//		}
//		//////////////////////////////////////// compute Stochastic Mapping
//		if(!gainLossOptions::_gainLossDist){
//			startComputePosteriorExpectationOfChange(seqReOpened,trSim,spSim,LpostPerCatSim,unObservableDataSim,outDirSeq);
//			LpostPerCatSim.clear();	// when cleared - each replicate will recompute the _LpostPerCat
//		}
//		else{
//			startComputePosteriorExpectationOfChange(seqReOpened,trSim,spVVecSim,gainDistSim,lossDistSim,LpostPerSpPerCat,unObservableDataSim,outDirSeq);
//			LpostPerSpPerCat.clear();	// when cleared - each replicate will recompute the _LpostPerSpPerCat
//		}
//		time(&t2);
//		LOGnOUT(4,<<"Replicate SimultePosteriorExpectationOfChange RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
//	}	
//}


