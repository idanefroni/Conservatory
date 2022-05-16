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
#include "computeCorrelations.h"
#include "gainLossUtils.h"
#include "gainLossAlphabet.h"


/********************************************************************************************
computeCorrelations
Input: _expChanges_PosNodeXY - required,
if _expChanges_PosNodeXY_B not NULL, compute correlation against this data

1. Compute correlation before simulations, based only on real dataset 
(R is computed for each pair in real data)
startComputeAmongSitesCorrelations()
 correl->runComputeCorrelations() // with Real data

Perform several iteration of simulations:
startParametricBootstapCorrelation()
 Foreach iteration of simulations:
	1.1. simulated data with same model

	2. Compute correlation of simulated data
	computeCoEvolutionScoresBasedOnSimulatedDataCoMap()
	 2.1 fill LpostPerCat using rate4site or GL4site
	 2.2 fill expChanges_PosNodeXY_Sim stochastic mapping using computeCountsGL
	 2.3 new computeCorrel object with both real and simulated data used:
	  2.3.1 runComputeCorrelations
	  2.3.2 sort - produceSortedVectorsOfAllCorrelations
	  2.3.3 bins - produceSortedVectorsOfCorrelationsBinedByRate
	  2.3.4 pVal - computedCorrelationsPValBasedOnSimulatedDataCoMapBins
	  2.3.5 FDR pVals2qVals
	  2.3.6 printComputedCorrelationsData (smart print of map values) 

*********************************************************************************************/
computeCorrelations::computeCorrelations(tree& tr,  string& outDir, VVVVdouble* expChanges_PosNodeXY, VVVVdouble* expChanges_PosNodeXY_B):
_tr(tr),_outDir(outDir) 
{
	_expChanges_PosNodeXY = *expChanges_PosNodeXY;
	
	// Type of correlation - assume _EventTypes =(gain, loss, both) and if less options, the last ones are missing
	
	if(gainLossOptions::_isCorrelateWithPearson)
		_isPearson.push_back(true);
	if(gainLossOptions::_isCorrelateWithSpearman)
		_isPearson.push_back(false);
	if(_isPearson.size()==0){
		_isPearson.push_back(true);
		LOGnOUT(4,<<"Pearson correlation is compted since no option is selected"<<endl);
	}

	if(gainLossOptions::_isOnlyCorrelateWithBoth){
		_EventTypes.push_back("gain");
		_EventTypes.push_back("loss");
		_EventTypes.push_back("both");
	}
	else{
		_EventTypes.push_back("gain");
		if(gainLossOptions::_isAlsoCorrelateWithLoss)
			_EventTypes.push_back("loss");
		if(gainLossOptions::_isAlsoCorrelateWithBoth)
			_EventTypes.push_back("both");
	}	
	for (int i = 0; i <_EventTypes.size(); ++i){
		_EventTypesMap[_EventTypes[i]]=i;
		map<string, int> FromTo;
		if(_EventTypes[i] == "gain"){
			FromTo["from"]=0;
			FromTo["to"]=1; 		}
		else if(_EventTypes[i] == "loss"){
			FromTo["from"]=1;
			FromTo["to"]=0;
		}else if(_EventTypes[i] == "both"){
			LOGnOUT(4,<<"Event _EventTypesFromTo is not applicable for "<<_EventTypes[i]<<" both 0->1 and 1->0 are computed"<<endl);
			break;
		}
		_EventTypesFromTo[_EventTypes[i]] = FromTo;
		LOGnOUT(4,<<"Event Type="<<_EventTypes[i]<<endl);
	}

	if(expChanges_PosNodeXY_B){
		_expChanges_PosNodeXY_B = *expChanges_PosNodeXY_B;
		_isTwoSetsOfInputForCorrelation = true;
	}else{
		_isTwoSetsOfInputForCorrelation = false;
	}
	_numOfSamplesInLowRateFirstBin = (int)min(100.0, (double)(_expChanges_PosNodeXY.size()/10.0)); // thus, the best p-value for low Rate 0.01
	if(_numOfSamplesInLowRateFirstBin<1)
		_numOfSamplesInLowRateFirstBin = 1;
	LOGnOUT(4,<<"Lowest pVal for correlation with rate below simulations is "<<1.0/_numOfSamplesInLowRateFirstBin<<endl);
}

/********************************************************************************************
*********************************************************************************************/
computeCorrelations::~computeCorrelations(){
	//clearVVVV(_jointProb_PosNodeXY);
}

/********************************************************************************************
*********************************************************************************************/
computeCorrelations& computeCorrelations::operator=(const computeCorrelations &other){
	if (this != &other) {              // Check for self-assignment
	}
	return *this;
}

/********************************************************************************************
Compute the Pearson / Spearman correlation among sites.
*********************************************************************************************/
void computeCorrelations::runComputeCorrelations(const Vint& selectedPositions, const Vint& numOfGapsTillSite, const bool isNormalizeForBranch)
{
	LOGnOUT(4,<<endl<<"runComputeCorrelations..."<<endl);
	time_t t1,t2;
	time(&t1);

	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node
	int numOfSitesSelected = selectedPositions.size();
	int numOfpositionsIn_A = _expChanges_PosNodeXY.size();
	int numOfpositionsIn_B;
	if(_isTwoSetsOfInputForCorrelation)
		numOfpositionsIn_B = _expChanges_PosNodeXY_B.size();
	else
		numOfpositionsIn_B = _expChanges_PosNodeXY.size(); // if B is not given, it's copy

	if(_isTwoSetsOfInputForCorrelation)
		LOGnOUT(3, <<"NOTE: Two seperate dataset input.\n Compute correl for selectedSites="<<numOfSitesSelected<<" subset of A="<<numOfpositionsIn_A<<" against B="<<numOfpositionsIn_B<<endl);
		
	//// Mapping vectors
	LOGnOUT(4, <<"Fill events vectors..."<<endl);
	// Expectation, keep the duplicated code. Maybe update later
	_expPerPosPerBranchVec.resize(_EventTypes.size());
	_expPerPosPerBranchVec_B.resize(_EventTypes.size());
	for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){
		if(*evnt == "gain" || *evnt == "loss")
			fillMapValPerPosPerBranch(_expPerPosPerBranchVec[_EventTypesMap[*evnt]],*evnt,_expChanges_PosNodeXY,isNormalizeForBranch);	// fill _expPerPosPerBranchVec
		if(*evnt == "both"){
			if(_EventTypes.size()<3)
				errorMsg::reportError("Error: correlation for _EventTypes=both with less than 3 options assume:(gain, loss, both)");
			_expPerPosPerBranchVec[_EventTypesMap[*evnt]] = _expPerPosPerBranchVec[_EventTypesMap["gain"]]; // gain and loss appended (double size vector)
			appendVectors(_expPerPosPerBranchVec[_EventTypesMap[*evnt]], _expPerPosPerBranchVec[_EventTypesMap["loss"]]);
		}
		if(_isTwoSetsOfInputForCorrelation){
			if(*evnt == "gain" || *evnt == "loss")
				fillMapValPerPosPerBranch(_expPerPosPerBranchVec_B[_EventTypesMap[*evnt]],*evnt,_expChanges_PosNodeXY_B,isNormalizeForBranch);	// 
			if(*evnt == "both"){
				_expPerPosPerBranchVec_B[_EventTypesMap[*evnt]] = _expPerPosPerBranchVec_B[_EventTypesMap["gain"]]; // gain and loss appended (double size vector)
				appendVectors(_expPerPosPerBranchVec_B[_EventTypesMap[*evnt]], _expPerPosPerBranchVec_B[_EventTypesMap["loss"]]);
			}
		}else{
			_expPerPosPerBranchVec_B = _expPerPosPerBranchVec;
		}		
	}
	
	if(gainLossOptions::_isOnlyCorrelateWithBoth){ // if "both", gain and loss were used only for the fill-up.
		while(*_EventTypes.begin() == "gain" || *_EventTypes.begin() == "loss")
		 _EventTypes.erase (_EventTypes.begin());		
	}

	//// correlation vectors, filled below
	LOGnOUT(6, <<"Resize correlation vectors vectors"<<endl);
	int numberOfCorrelations = _isPearson.size()*_EventTypes.size();
	_correlationsPerSitePerPosVec.resize(numberOfCorrelations);
	for (int typeC = 0; typeC <numberOfCorrelations; ++typeC)
		resizeMatrix(_correlationsPerSitePerPosVec[typeC], numOfSitesSelected, numOfpositionsIn_B);

	//for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
	//	for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
	//		LOGnOUT(4, <<vecIndex<<" - Compute correl isSpearman="<<*it<<" with type="<<*evnt<<endl);
	//		vecIndex++;
	//	}
	//}

	int vecIndex=0;
	for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
		//int typeIndex=0;
		for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
			Vdouble correlationVecAve; // per correlation type, each item is the Mean for a selected position again all
			Vdouble correlationVecMedian; // per correlation type, each item is the Median for a selected position again all
			LOGnOUT(4, <<"Compute correlation isPearson="<<*it<<" with type="<<*evnt<<endl);
			for (int selectedSiteIndex = 0; selectedSiteIndex <numOfSitesSelected; ++selectedSiteIndex){
				if(selectedSiteIndex%100==0)
					cout<<"*";
				int selectedSite = selectedPositions[selectedSiteIndex];
				int selectedSiteRemovedGaps = selectedSite- numOfGapsTillSite[selectedSiteIndex];
				fillCorrPerSelectedSites(_correlationsPerSitePerPosVec[vecIndex][selectedSiteIndex],_expPerPosPerBranchVec[_EventTypesMap[*evnt]],_expPerPosPerBranchVec_B[_EventTypesMap[*evnt]],selectedSiteRemovedGaps,(*it)); // expPerPosPerBranchVec still have gain,loss,both
				
				correlationVecAve.push_back(computeAverage((_correlationsPerSitePerPosVec[vecIndex][selectedSiteIndex])));
				correlationVecMedian.push_back(computeMedian((_correlationsPerSitePerPosVec[vecIndex][selectedSiteIndex])));
			}
			cout<<"\n"; // end of "*" for this correlation type
			if(gainLossOptions::_selectedSitesForCorrelation=="")
				LOGnOUT(4, <<"Correlation coefficient (mean of Val=Mean/Median per selected) Mean="<<computeAverage(correlationVecAve)<<" Median="<<computeAverage(correlationVecMedian)<<endl);			
			//typeIndex++;
			vecIndex++;
		}
	}	
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE computeCorrelations::computeNminPerPair(const int site_A, const int site_B, const int typeIndex, const VVVdouble&  exp_PosXY){
	MDOUBLE NminVal = 0;
	MDOUBLE siteA_Rate;
	MDOUBLE siteB_Rate;

	if(typeIndex == 2 || gainLossOptions::_isOnlyCorrelateWithBoth){ // both
		MDOUBLE siteA_Gain = exp_PosXY[site_A][0][1];
		MDOUBLE siteA_Loss =  exp_PosXY[site_A][1][0];
		siteA_Rate = computeNminRforCorrelWithGainAndLoss(siteA_Gain,siteA_Loss);

		MDOUBLE siteB_Gain = exp_PosXY[site_B][0][1];
		MDOUBLE siteB_Loss = exp_PosXY[site_B][1][0];
		siteB_Rate = computeNminRforCorrelWithGainAndLoss(siteB_Gain,siteB_Loss);
	}else{ 
		string type =  _EventTypes[typeIndex];
		int from = _EventTypesFromTo[type]["from"];
		int to = _EventTypesFromTo[type]["to"];
		siteA_Rate = exp_PosXY[site_A][from][to];
		siteB_Rate = exp_PosXY[site_B][from][to];		
	}

	NminVal = min(siteA_Rate, siteB_Rate);
	return NminVal;


}



/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::produceSortedVectorsOfAllCorrelations(Vdouble& rate4siteSim){
	LOGnOUT(4,<<endl<<"produceSortedVectorsOfAllCorrelations for simulated data (sort rates and correlations (paired))..."<<endl);
	time_t t1,t2;
	time(&t1);

	if(rate4siteSim.size()==0) // no Rate4site input is provided
		computeRateValPerPos(_expChanges_PosNodeXY,_exp_PosXY);
	int numberOfcorrelationVec = _correlationsPerSitePerPosVec.size();
	int numOfSites_A = _correlationsPerSitePerPosVec[0].size();
	int numOfSites_B = _correlationsPerSitePerPosVec[0][0].size();
	_pairWiseCorrelationsAndNminSim.resize(numberOfcorrelationVec);
	_NminSortedSim.resize(numberOfcorrelationVec);
	
	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		LOGnOUT(4,<<"  ***   corIndex="<<corIndex<<endl);
		int typeIndex = corIndex % _EventTypes.size();	// in case both Spearman and pearson are used	
		int indexAll = 0; 
		Vdouble correlations;
		Vdouble Nmins;
		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			for (int site_B = site_A; site_B <numOfSites_B; ++site_B){
				if(site_A == site_B)
					continue;
				MDOUBLE correlVal = _correlationsPerSitePerPosVec[corIndex][site_A][site_B];
				correlations.push_back(correlVal);
				MDOUBLE NminVal=0;
				if(rate4siteSim.size()==0)
					NminVal = computeNminPerPair(site_A, site_B, typeIndex, _exp_PosXY);					
				else
					NminVal = min(rate4siteSim[site_A],rate4siteSim[site_B]);

				Nmins.push_back(NminVal);				
				indexAll++;
			}
		}
		vector< vecElem<MDOUBLE> > orderVecNmin;
		orderVec(Nmins, orderVecNmin);

		//resizeMatrix(_pairWiseCorrelationsAndNminSim[corIndex],  2  ,indexAll); // pairWiseCorrelationsAndNmin[corrIndex][pairIndex][0/1][val]
		resizeMatrix(_pairWiseCorrelationsAndNminSim[corIndex],  1  ,indexAll); // pairWiseCorrelationsAndNmin[corrIndex][pairIndex][0/1][val]
		for (int i = 0; i <indexAll; ++i){
			//_pairWiseCorrelationsAndNminSim[corIndex][0][i] = orderVecNmin[i].getValue();
			//_pairWiseCorrelationsAndNminSim[corIndex][1][i] = correlations[orderVecNmin[i].getPlace()];
			_pairWiseCorrelationsAndNminSim[corIndex][0][i] = correlations[orderVecNmin[i].getPlace()];
		}		
		_NminSortedSim[corIndex] = Nmins; // vector copy no resize?
		sort( _NminSortedSim[corIndex].begin(),_NminSortedSim[corIndex].end() );

		LOGnOUT(4,<<"\nSimulated Data correlations frequencies:"<<endl);
		LOGnOUT(4,<<"num of correlations="<<correlations.size()<<endl);
		sort(correlations.begin(),correlations.end());
		printCorrelationsFrequencies(correlations);
		
		LOGnOUT(4,<<"Finish sorting "<<indexAll<<" pairs of correlating sites"<<endl);
		LOGnOUT(4,<<"Minimal rate (Nmin)="<<*(_NminSortedSim[corIndex].begin() )<<" max="<<*(_NminSortedSim[corIndex].end()-1)<<endl);
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes (sort vectors)"<<endl);
}


/********************************************************************************************
Each bin the index within sorted Nmin vector (Rate) is progressed (+++)
if _isSortVectorOfCorrelationsBinsByLowerRateBound
1. => the binLimit = LowLimit = f(index), (computed before index +++)
 UpLimit is maxLimit for all bins
else
2. => the binLimit = UpLimit = f(index),  (computed after index +++)
 LowLimit is the UpLimit of previous bin

Note: Two versions exists for the _isSortVectorOfCorrelationsBinsByLowerRateBound (which is not the default)
The assumption (which appears correct in Corr~1) is that the probability of high correlation by chance is smaller with higher rate.

Thus, in the modified one (16/05/12) higher rates had more simulations to compare with (compare with all sim. with lower rate)
(compared with all those below.)

In previous version, higher rates had less simulations to compare with (compare with all sim. with higher rate), but for high Obs. rate,
comparison with "lower rate bins" was allowed to avoid the paradox of smaller pVal of pairs with low rate (with "while" mechanism)

*********************************************************************************************/
int computeCorrelations::produceSortedVectorsOfCorrelationsBinedByRate(MDOUBLE medianNminOfRealData, ofstream* simCorrelStream){
	LOGnOUT(4,<<endl<<"produceSortedVectorsOfCorrelationsBinedByRate for simulated data..."<<endl);
	time_t t1,t2;
	time(&t1);

	int numberOfHighCorrInSimulationOfMedianNminBin = 0;
	int numOfBins = gainLossOptions::_numOfBinsInParametricBootstrapSimulations;

	pair<vector<double>::iterator,vector<double>::iterator> bounds;	
		
	int numberOfcorrelationVec = _correlationsPerSitePerPosVec.size();
	_correlationSubSetsNminLimitValues.resize(numberOfcorrelationVec);
	_correlationsSubSets.resize(numberOfcorrelationVec);
	_extremeValDistributions.resize(numberOfcorrelationVec);
	
	int numOfSimulatedTotalPairs = _NminSortedSim[0].size(); // same for all CorrTypes
	LOGnOUT(4,<<"Num of pairs in simulations="<<numOfSimulatedTotalPairs<<endl);

	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		LOGnOUT(4,<<"For corIndex="<<corIndex<<endl);
		int typeIndex = corIndex % _EventTypes.size();	// in case both Spearman and pearson are used
		MDOUBLE Nmin_min = 0;
		if(Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair")>0){
			Nmin_min = Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
			LOGnOUT(4,<<"Nmin_min by threshold="<<Nmin_min<<endl);
		}
		else
			Nmin_min = *_NminSortedSim[corIndex].begin();	

		bounds = equal_range (_NminSortedSim[corIndex].begin(), _NminSortedSim[corIndex].end(), Nmin_min);
		int IndexOfValAboveNmin_min = int(bounds.first - _NminSortedSim[corIndex].begin());
		LOGnOUT(4,<<"Nmin_min in dataset="<<*(_NminSortedSim[corIndex].begin()+IndexOfValAboveNmin_min)<<endl);
		
		int numOfSimulationPairs = _NminSortedSim[corIndex].size()-IndexOfValAboveNmin_min;
		if(numOfSimulationPairs==0)
			errorMsg::reportError("_minExpThresholdForPValComputationForCorrelatingPair is too high, no simulations above that value");

		LOGnOUT(4,<<"Num of pairs above threshold (valid)="<<numOfSimulationPairs<<endl<<endl);
		int numOfSamplesInBin;
		bool randomOverlapPerIteration = true;
		MDOUBLE overlap = gainLossOptions::_relativeSizeOfOverLappedBins;
		if(randomOverlapPerIteration)
			overlap = gainLossOptions::_relativeSizeOfOverLappedBins + talRandom::giveRandomNumberBetweenTwoPoints(-0.1, 0.1);
		

		LOGnOUT(4,<<"Size of overlapped bins ="<<overlap<<endl<<endl);
		
		if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound)
			numOfSamplesInBin = (int)(numOfSimulationPairs * overlap);
		else
			numOfSamplesInBin = (int)(numOfSimulationPairs/numOfBins);


		MDOUBLE Nmin_max =  *(_NminSortedSim[corIndex].end()-1);
		MDOUBLE w_range = (Nmin_max-Nmin_min)/numOfBins;

		int numOfSamplesInCurrBin = 0;
		MDOUBLE Nmin_lower = Nmin_min;
		MDOUBLE Nmin_upper = Nmin_min;
		MDOUBLE Nmin_mid = Nmin_min; // use with Mid Boundary with overlap

		MDOUBLE NminPerBin = 0;

		int indexOfSamplesForBin = IndexOfValAboveNmin_min;
		int indexOfSamplesForBinPrev = IndexOfValAboveNmin_min;
		int indexOfSamplesForBinUpper = IndexOfValAboveNmin_min;
		int indexOfSamplesForBinMid = IndexOfValAboveNmin_min;

		_correlationsSubSets[corIndex].resize(numOfBins+1); // the actual size may be smaller, if break
		//_correlationSubSetsNminLimitValues[corIndex].resize(numOfBins+1); // to Zero bin
		
		vector<MDOUBLE>::iterator it = _pairWiseCorrelationsAndNminSim[corIndex][0].begin(); // correlation part of vector		
		

		// elevate Nmin Threshold if: (A) freqOfHighCorr was too high (B) freqOfHighCorr is reduced consistently with higher Nmin (C) new Nmin is lower than medianNminOfRealData
		MDOUBLE minExpTBeforeChange = (double)Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");
		MDOUBLE freqOfHighCorr = 0;
		MDOUBLE freqOfHighCorrPrev = 0;
		MDOUBLE expextedFreq;
		if(gainLossOptions::_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel){
			int isHigherNminReducedFreqOfHighCorr = false;
			int numOfBranches = _expChanges_PosNodeXY[0].size();
			int numOfBranches99 = (int)(numOfBranches*0.99);
			int combo99 = BinomialCoeff(numOfBranches,numOfBranches99);
			expextedFreq = 0.01 /  (double)combo99;			
			LOGnOUT(3,<<"Allowed fraction of high correlation. Computed with number of branch "<<numOfBranches<<" is "<<expextedFreq<<endl);
		}

		for (int binIndex = 0; binIndex < numOfBins; ++binIndex){
			int indexOfIncrementation;
			int numOfSamplesToDivideAmondBins =  numOfSimulationPairs-numOfSamplesInBin;
			
			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound){
				indexOfIncrementation = (int)(numOfSamplesToDivideAmondBins*binIndex/(numOfBins-1) );
				indexOfSamplesForBinPrev = indexOfIncrementation +IndexOfValAboveNmin_min;
			}
			else
				indexOfSamplesForBinPrev = indexOfSamplesForBin;

			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound)
				Nmin_lower = *(_NminSortedSim[corIndex].begin()+indexOfSamplesForBinPrev); // Low is computed before bin-related ++ of index			
			else			
				Nmin_lower = Nmin_upper;

			// +++ index for bin
			if(gainLossOptions::_isDivideBinsByRange){
				NminPerBin =  Nmin_min +(w_range*binIndex);
				bounds = equal_range (_NminSortedSim[corIndex].begin(), _NminSortedSim[corIndex].end(), NminPerBin); // assume sorted, 
				indexOfSamplesForBin = int(bounds.first - _NminSortedSim[corIndex].begin())-1;				// Nmin_endIndex = int(bounds.second - _NminSortedSim[corIndex].begin());				
			}
			else if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound)				
				indexOfSamplesForBin = indexOfIncrementation +numOfSamplesInBin +IndexOfValAboveNmin_min -1; // 			
			else
				indexOfSamplesForBin = (int)(numOfSimulationPairs*(binIndex+1)/numOfBins) +IndexOfValAboveNmin_min -1; // is -1 for bin=0			

			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound)
				Nmin_lower = *(_NminSortedSim[corIndex].begin()+indexOfSamplesForBin);			


			// compute numOfSamples per Bin
			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound){
				numOfSamplesInCurrBin =  indexOfSamplesForBin; //_NminSortedSim[corIndex].size()-indexOfSamplesForBinPrev;
				Nmin_upper = Nmin_max; // UpperIsFixedAtMax
				indexOfSamplesForBinUpper = indexOfSamplesForBin; //_NminSortedSim[corIndex].size()-1;			
			}			
			else{
				numOfSamplesInCurrBin = indexOfSamplesForBin-indexOfSamplesForBinPrev;
				Nmin_upper = *(_NminSortedSim[corIndex].begin()+indexOfSamplesForBin); // Up is computed after bin-related ++ of index;
				indexOfSamplesForBinUpper = indexOfSamplesForBin;
			}			

			if(numOfSamplesInCurrBin<10) //  no samples in this range, for median
				break;
			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound && numOfSamplesInCurrBin<numOfSimulationPairs*0.05 ) //  at least 5% of simulation to start new bin, otherwise previous was last
				break;
			
			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound){
				indexOfSamplesForBinMid = (int)(indexOfSamplesForBinUpper+indexOfSamplesForBinPrev)/2;
				Nmin_mid = *(_NminSortedSim[corIndex].begin()+indexOfSamplesForBinMid); 
			}

			// assign limits utility vector
			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound)
				_correlationSubSetsNminLimitValues[corIndex].push_back(Nmin_lower); // UpperIsFixedAtMax
			else if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound)
				_correlationSubSetsNminLimitValues[corIndex].push_back(Nmin_mid);				
			else
				_correlationSubSetsNminLimitValues[corIndex].push_back(Nmin_upper);				
			
			_correlationsSubSets[corIndex][binIndex].resize(numOfSamplesInCurrBin);
			copy(it+indexOfSamplesForBinPrev, it+indexOfSamplesForBinUpper ,_correlationsSubSets[corIndex][binIndex].begin());
			sort(_correlationsSubSets[corIndex][binIndex].begin(),_correlationsSubSets[corIndex][binIndex].end());


			extremeValDistribution distr;
			MDOUBLE averageCorr = computeAverage(_correlationsSubSets[corIndex][binIndex]);
			MDOUBLE stdCorr = computeStd(_correlationsSubSets[corIndex][binIndex]);
			distr.fitParametersFromMoments(averageCorr, stdCorr);
			_extremeValDistributions[corIndex].push_back(distr);

			pair<vector<double>::iterator,vector<double>::iterator> boundsOne;
			boundsOne = equal_range (_correlationsSubSets[corIndex][binIndex].begin(),_correlationsSubSets[corIndex][binIndex].end(), 0.99999);
			int indexOfpairEq1_first =  int(boundsOne.first - _correlationsSubSets[corIndex][binIndex].begin());
			int numOfpairWithCorrEq1 = numOfSamplesInCurrBin - indexOfpairEq1_first;
			
			boundsOne = equal_range (_correlationsSubSets[corIndex][binIndex].begin(),_correlationsSubSets[corIndex][binIndex].end(), 0.99);
			int indexOfpairEq99_first =  int(boundsOne.first - _correlationsSubSets[corIndex][binIndex].begin());
			int numOfpairWithCorrEq99 = numOfSamplesInCurrBin - indexOfpairEq99_first;

			boundsOne = equal_range (_correlationsSubSets[corIndex][binIndex].begin(),_correlationsSubSets[corIndex][binIndex].end(), 0.9); 
			int indexOfpairEq9_first =  int(boundsOne.first - _correlationsSubSets[corIndex][binIndex].begin());
			int numOfpairWithCorrEq9 = numOfSamplesInCurrBin - indexOfpairEq9_first;

			// elevate Nmin Threshold if: (A) freqOfHighCorr was too high (B) freqOfHighCorr is reduced consistently with higher Nmin (C) new Nmin is lower than medianNminOfRealData
			if(gainLossOptions::_isUpdateMinExpThresholdGivenHighFractionOfHighCorrel){
				freqOfHighCorrPrev = freqOfHighCorr;
				freqOfHighCorr = (double)numOfpairWithCorrEq99/numOfSamplesInCurrBin;				
				if(freqOfHighCorr>expextedFreq  && freqOfHighCorr<freqOfHighCorrPrev && Nmin_lower < medianNminOfRealData){
					LOGnOUT(3,<<"Fraction of high (0.99) correlation prev="<<freqOfHighCorrPrev<<" reduced to "<<freqOfHighCorr<<endl);
					LOGnOUT(3,<<"  Update MinExpThreshold Given highCorrlation in previous Nmin to "<<Nmin_lower<<endl);
					Parameters::updateParameter("_minExpThresholdForPValComputationForCorrelatingPair",double2string(Nmin_lower).c_str());
				}
				if(freqOfHighCorr>freqOfHighCorrPrev){ // revert back
					LOGnOUT(3,<<"Fraction of high (0.99) correlation prev="<<freqOfHighCorrPrev<<" elevated to "<<freqOfHighCorr<<endl);
					LOGnOUT(3,<<"  Revert to "<<minExpTBeforeChange<<endl);
					Parameters::updateParameter("_minExpThresholdForPValComputationForCorrelatingPair",double2string(minExpTBeforeChange).c_str());
				}
			}
			
			if(Nmin_lower>=medianNminOfRealData)
				numberOfHighCorrInSimulationOfMedianNminBin = max((double)numberOfHighCorrInSimulationOfMedianNminBin,(double)numOfpairWithCorrEq1);
			
			*simCorrelStream<<"Bin = "<< binIndex+1 <<"\n";
			printCorrelationsFrequencies(_correlationsSubSets[corIndex][binIndex], simCorrelStream);
			
			LOGnOUT(4,<<binIndex+1<<" Bin.\t#samples=\t"<<numOfSamplesInCurrBin<<".\tFrom rate:\t"<<Nmin_lower<<"\t-\t"<<Nmin_upper
				<<".\tis with corr:\t"<<*_correlationsSubSets[corIndex][binIndex].begin()<<"\t-\t"<<*(_correlationsSubSets[corIndex][binIndex].end()-1)				
				<<"\tAve=\t"<<computeAverage(_correlationsSubSets[corIndex][binIndex])<<"\tMedian=\t"<< computeMedian(_correlationsSubSets[corIndex][binIndex]) 
				<<"\tratioOfpairWithCorrEq1=\t"<<(double)numOfpairWithCorrEq1/numOfSamplesInCurrBin<<" ("<<numOfpairWithCorrEq1<<")"				
				<<"\tratioOfpairWithCorrEq0.99=\t"<<(double)numOfpairWithCorrEq99/numOfSamplesInCurrBin<<" ("<<numOfpairWithCorrEq99<<")"
				<<"\tratioOfpairWithCorrEq0.9=\t"<<(double)numOfpairWithCorrEq9/numOfSamplesInCurrBin<<" ("<<numOfpairWithCorrEq9<<")"<<endl);

			if(gainLossOptions::_isSortVectorOfCorrelationsBinsByMidRateBound)
				LOGnOUT(4,<<"  Mid rate= "<<Nmin_mid<<endl);

			// Util
			//bool isPrintCorrListForEachBin = false;
			//if(isPrintCorrListForEachBin){
			//	string debugS = _outDir + "//"+int2string(corIndex)+int2string(binIndex)+ "Rofbins.txt"; // D
			//	ofstream debugSStream(debugS.c_str()); // D
			//	debugSStream<<" Bin "<<binIndex<<" from rate: "<<Nmin_lower<<" to "<<Nmin_max<<endl;
			//	for(vector<double>::iterator it = _correlationsSubSets[corIndex][binIndex].begin(); it<_correlationsSubSets[corIndex][binIndex].end();++it){
			//		debugSStream<<*it<<"\n";
			//	}
			//}
		}		 
	}
	_pairWiseCorrelationsAndNminSim.clear(); // clear huge vector when not required
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	return numberOfHighCorrInSimulationOfMedianNminBin;
}

/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::printCorrelationsFrequencies(Vdouble& correlationsVecSorted, ofstream* simCorrelStream){
	
	float valsToCheck [] = {0.95,0.99,0.999,0.999999999}; // NOTE - if change size => change in loop!
	
	int numOfCorrelations = correlationsVecSorted.size();	
	pair<vector<double>::iterator,vector<double>::iterator> bounds;
	if(!simCorrelStream==NULL)
		*simCorrelStream<<"Corr eq/above\tratioOfCorAbove\tnumAboveEq\n";
	else
		LOGnOUT(4,<<"Corr eq/above\tratioOfCorAbove\tnumAboveEq"<< endl);
	for (MDOUBLE val=-0.9; val<=0.9; val+=0.1){
		bounds = equal_range (correlationsVecSorted.begin(), correlationsVecSorted.end(), val);			
		int lastIndexWithPValBiggerThanThreshold = int(bounds.first - correlationsVecSorted.begin());
		int numAboveEq = numOfCorrelations-lastIndexWithPValBiggerThanThreshold;
		MDOUBLE ratioOfCorAbove = double(numAboveEq)/numOfCorrelations;
		MDOUBLE rounded = floorf(val * pow(10.0,4) + 0.5) / pow(10.0,4); // if not rounded, perfect correlations may return 1.000002, for example
		if(!simCorrelStream==NULL)
			*simCorrelStream<<rounded<<"\t"<<ratioOfCorAbove<<"\t("<<numAboveEq<<")\n";
		else
			LOGnOUT(4,<<rounded<<"\t"<<ratioOfCorAbove<<"\t("<<numAboveEq<<")"<< endl);
	}	
	for (int i=0; i<4; ++i){
		bounds = equal_range (correlationsVecSorted.begin(), correlationsVecSorted.end(), valsToCheck[i]);			
		int lastIndexWithPValBiggerThanThreshold = int(bounds.first - correlationsVecSorted.begin());
		int numAboveEq = numOfCorrelations-lastIndexWithPValBiggerThanThreshold;
		MDOUBLE ratioOfCorAbove = double(numAboveEq)/numOfCorrelations;
		if(!simCorrelStream==NULL)
			*simCorrelStream<<valsToCheck[i]<<"\t"<<ratioOfCorAbove<<"\t("<<numAboveEq<<")\n";
		else
			LOGnOUT(4,<<valsToCheck[i]<<"\t"<<ratioOfCorAbove<<"\t("<<numAboveEq<<")"<< endl);
	}
	if(!simCorrelStream==NULL)
		*simCorrelStream<<"\n";
	else
		LOGnOUT(4,<< endl);
}




/********************************************************************************************
*********************************************************************************************/
int computeCorrelations::computedCorrelationsPValBasedOnSimulatedDataCoMapBins(VVVdouble& correlationPerSitePerPosReal,vector<vector<bool> >& isComputePairWithRateAboveNim,VVVVdouble& expChanges_PosXYReal, VVVdouble& correlationPerSitePerPos_Pval
		,map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData, Vdouble& rate4siteReal, Vint& selectedSites, Vint& numOfGapsTillSite, Vint& evolvingSites, bool isLastIteration){
	LOGnOUT(4,<<endl<<"computedCorrelationsPValBasedOnSimulatedDataCoMapBins..."<<endl);
	time_t t1,t2;
	time(&t1);				
	
	int numOfpairsWithRateAboveMinRequiredExp = 0;
	string pairWiseCorrelationsAndNmin = _outDir + "//" + "pairWiseCorrelationsAndNmin.txt";
	ofstream corrSigStream(pairWiseCorrelationsAndNmin.c_str());
	corrSigStream<<"site_A"<<"\t"<<"site_B"<<"\t"<<"Nmin_obs"<<"\t"<<"Corr_obs"<<"\n";	

	int numberOfcorrelationVec = correlationPerSitePerPosReal.size();
	int numOfSites_A = correlationPerSitePerPosReal[0].size();
	int numOfSites_B = correlationPerSitePerPosReal[0][0].size();
	_corrVector.resize(numberOfcorrelationVec);

	VVVdouble map_PosXY;
	if(rate4siteReal.size()==0)
		computeRateValPerPos(expChanges_PosXYReal,map_PosXY);

	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		LOGnOUT(4,<<"  ***   corIndex="<<corIndex<<endl);
		int typeIndex = corIndex % _EventTypes.size();	// in case both Spearman and pearson are used				
		int numOfpairsWithRateBelowSimulation = 0;
		int numOfpairsWithRateAboveSimulation = 0;
		int numOfpairsWithRateBelowMinRequiredExp = 0;
		int pairNum = 0;
		bool computePValForPairWithNminAboveMin = true;
		MDOUBLE minExpThresholdForPValComputationForCorrelatingPair;
		minExpThresholdForPValComputationForCorrelatingPair = Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair");

		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			int site_A_original = selectedSites[site_A];
			int site_A_RemovedGaps = site_A_original- numOfGapsTillSite[site_A];
			if((site_A)%100==0)
				cout<<"*";
			for (int site_B = site_A; site_B <numOfSites_B; ++site_B){
				int site_B_original = evolvingSites[site_B];
				computePValForPairWithNminAboveMin = true; // reset as new for each pair
				if(site_A_original == site_B_original){
					correlationPerSitePerPos_Pval[corIndex][site_A][site_B] = 0;
					continue;
				}				
				pairNum++;
				MDOUBLE Corr_obs = correlationPerSitePerPosReal[corIndex][site_A][site_B]; // Real correlation from Input variable
				_corrVector[corIndex].push_back(Corr_obs);
				MDOUBLE Nmin_obs = 0;

				if(rate4siteReal.size()==0)
					Nmin_obs = computeNminPerPair(site_A_RemovedGaps, site_B, typeIndex, map_PosXY);					
				else
					Nmin_obs = min(rate4siteReal[site_A_original] , rate4siteReal[site_B_original]);

				if(gainLossOptions::_isPrintpairWiseCorrelationsAndNmin)					
					corrSigStream<<site_A_original+1<<"\t"<<site_B_original+1<<"\t"<<Nmin_obs<<"\t"<<Corr_obs<<"\n";

				// find the bin with highest simulated Rate suitable of obsNmin
				int binForNmin_obs = 0;
				while(Nmin_obs>=_correlationSubSetsNminLimitValues[corIndex][binForNmin_obs] && binForNmin_obs<_correlationSubSetsNminLimitValues[corIndex].size()-1)
					binForNmin_obs++;

				if(Nmin_obs < minExpThresholdForPValComputationForCorrelatingPair){
					computePValForPairWithNminAboveMin = false;
					numOfpairsWithRateBelowMinRequiredExp++;
				}				
				if(Nmin_obs<*(_NminSortedSim[corIndex].begin())){
					LOGnOUT(7,<<"WARN: low Nmin_obs="<<Nmin_obs<<" Since no simulation support this rate,  pVal computed as "<<1.0/_numOfSamplesInLowRateFirstBin<<" for site_A="<<site_A_original<<" and site_B="<<site_B_original<<" with corr="<<Corr_obs<<endl);					
					//computePValForPairWithNminAboveMin = false;
					if(corIndex == 0) // done for only one type of correlation
						numOfpairsWithRateBelowSimulation++;
				}
				if(Nmin_obs> *(_NminSortedSim[corIndex].end()-1)){
					LOGnOUT(7,<<"WARN: high Nmin_obs="<<Nmin_obs<<" pVal is computed with lower Nmin simulations as referece for site_A="<<site_A_original<<" and site_B="<<site_B_original<<" with corr="<<Corr_obs<<endl);
					if(corIndex == 0) // done for only one type of correlation
						numOfpairsWithRateAboveSimulation++;
				}				
				
				MDOUBLE pVal = 0.99;
				MDOUBLE pValEVD = 0.99;
				int NumberOfSimulationsInRange = 0;
				int NumberOfSimulationPointsMoreExtremeOrEqToCorr = 0;
				int prevNumberOfSimulationsInRange = 1;
				int prevNumberOfSimulationPointsMoreExtremeOrEqToCorr = 1;

				int numOfBinsWithLowerSignificance = 0; // allow 2 "bin-iteration" even with lower significance to mitigate chance of "missing" higher significance in lower bin
				bool isNextLowerBinAllowed = true;
				isComputePairWithRateAboveNim[site_A][site_B] = computePValForPairWithNminAboveMin;
				if(computePValForPairWithNminAboveMin){					
					//MDOUBLE pVal_prev = 1;
					//while(binForNmin_obs>=0 && isNextLowerBinAllowed ){
						//pVal_prev = pVal;
						prevNumberOfSimulationsInRange = NumberOfSimulationsInRange;
						prevNumberOfSimulationPointsMoreExtremeOrEqToCorr = NumberOfSimulationPointsMoreExtremeOrEqToCorr;
						NumberOfSimulationsInRange =  _correlationsSubSets[corIndex][binForNmin_obs].size();
						NumberOfSimulationPointsMoreExtremeOrEqToCorr = 0;
						pair<vector<double>::iterator,vector<double>::iterator> bounds;
						vector<double>::iterator startCorV =  _correlationsSubSets[corIndex][binForNmin_obs].begin();
						vector<double>::iterator endCorV =  _correlationsSubSets[corIndex][binForNmin_obs].end();
						bounds = equal_range (startCorV, endCorV, Corr_obs);				
						//cout << "bounds at positions " << int(bounds.first - startCorV) << " and " << int(bounds.second - startCorV) << endl;				
						NumberOfSimulationPointsMoreExtremeOrEqToCorr = NumberOfSimulationsInRange-int(bounds.first - startCorV);
						
						if(gainLossOptions::_isConsiderNegativeCorrelations){
							int NumberOfSimulationPointsMoreExtremeOrEqToCorrNegative = int(bounds.second - startCorV);
							NumberOfSimulationPointsMoreExtremeOrEqToCorr = min(NumberOfSimulationPointsMoreExtremeOrEqToCorr, NumberOfSimulationPointsMoreExtremeOrEqToCorrNegative);
							pVal = (double(NumberOfSimulationPointsMoreExtremeOrEqToCorr+1)/(NumberOfSimulationsInRange+1)) *2; // multiplied by 2, since it's two-sided
						}else
							pVal = double(NumberOfSimulationPointsMoreExtremeOrEqToCorr+1)/(NumberOfSimulationsInRange+1);

						if(gainLossOptions::_isCompExtremeValDistribution)
							pValEVD = 1- _extremeValDistributions[corIndex][binForNmin_obs].getCDF(Corr_obs);

						//if(pVal_prev<pVal){
						//	pVal= pVal_prev;
						//	NumberOfSimulationsInRange = prevNumberOfSimulationsInRange;
						//	NumberOfSimulationPointsMoreExtremeOrEqToCorr = prevNumberOfSimulationPointsMoreExtremeOrEqToCorr;
						//	++numOfBinsWithLowerSignificance; // the upper bin had more significant pVal. Count few such "steps down" and quite
						//}												
						//binForNmin_obs--;
						//if(!gainLossOptions::_isSortVectorOfCorrelationsBinsByLowerRateBound || (pVal<pVal_prev &&  numOfBinsWithLowerSignificance<3))
						//	isNextLowerBinAllowed = false;
					//}					
				}
				else					
					pVal = 1; // value that is not possible with computation				

				// Only pairs with pVal < cuttoff are re-computed in iterations 	
				if(pVal<=gainLossOptions::_pValueCutOffForBootStrap || gainLossOptions::_selectedSitesForCorrelation!="" ){ //TEMP for selected sites, fill map for all correlations
					//cout<<site_A<<" "<<site_B<<" "<<pVal<<" "<<Nmin_obs<<" "<<Corr_obs<<endl;
					correlationsData[site_A_original][site_B_original][int2string(corIndex)]["R"] = Corr_obs;
					correlationsData[site_A_original][site_B_original][int2string(corIndex)]["Rate"] = Nmin_obs;
					
					correlationsData[site_A_original][site_B_original][int2string(corIndex)]["SimTotal"] += NumberOfSimulationsInRange;
					correlationsData[site_A_original][site_B_original][int2string(corIndex)]["SimExtreme"] += NumberOfSimulationPointsMoreExtremeOrEqToCorr;
					pVal = (correlationsData[site_A_original][site_B_original][int2string(corIndex)]["SimExtreme"]+1)/(correlationsData[site_A_original][site_B_original][int2string(corIndex)]["SimTotal"]+1);
					correlationsData[site_A_original][site_B_original][int2string(corIndex)]["pVal"] = pVal;
					
					// take the higher pVal from all iterations

					bool isFirstEstimation = false;
					if(correlationsData[site_A_original][site_B_original][int2string(corIndex)]["SimTotal"] == NumberOfSimulationsInRange)
						isFirstEstimation = true;
					if(gainLossOptions::_isCompExtremeValDistribution){
						if(!isFirstEstimation)
							pValEVD = max(pValEVD, correlationsData[site_A_original][site_B_original][int2string(corIndex)]["pValEVD"]);					
						correlationsData[site_A_original][site_B_original][int2string(corIndex)]["pValEVD"] = pValEVD;
					}

				}
				if(gainLossOptions::_selectedSitesForCorrelation==""){
					correlationPerSitePerPos_Pval[corIndex][site_A][site_B] = pVal;
					correlationPerSitePerPos_Pval[corIndex][site_B][site_A] = pVal;
				}
			}
		}
		cout<<"\n";
		numOfpairsWithRateAboveMinRequiredExp = pairNum-numOfpairsWithRateBelowMinRequiredExp;
		LOGnOUT(4,<<"numOfpairs With Rate - Below Simulation="<<numOfpairsWithRateBelowSimulation<<" - Above Simulation="<<numOfpairsWithRateAboveSimulation<<endl);

		if(isLastIteration){
			LOGnOUT(4,<<"numOfpairs="<<pairNum<<endl);
			if(Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair") > 0){
				LOGnOUT(4,<<"numOfpairs With Rate below minimal Threshold="<<Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair")<<" are "<<numOfpairsWithRateBelowMinRequiredExp<<endl);
			}
			LOGnOUT(4,<<"\nReal Data correlations frequencies:"<<endl);
			LOGnOUT(4,<<"num of correlations="<<_corrVector[corIndex].size()<<endl);
			sort(_corrVector[corIndex].begin(),_corrVector[corIndex].end());
			printCorrelationsFrequencies(_corrVector[corIndex]);
		}
	}	
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return numOfpairsWithRateAboveMinRequiredExp;
}





/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::computedCorrelationsRankBasedOnSimulatedData(const Vint& selectedPositions, VVVdouble& correlationPerSitePerPos, VVVdouble& correlationPerSitePerPos_Simulations, VVVdouble& correlationPerSitePerPos_Pval){
	LOGnOUT(4,<<endl<<"computedCorrelationsRankBasedOnSimulatedData..."<<endl);
	time_t t1,t2;
	time(&t1);
	
	int numberOfcorrelationVec = correlationPerSitePerPos.size();
	int numOfSites_A = correlationPerSitePerPos[0].size();
	int numOfSites_B = correlationPerSitePerPos[0][0].size();
	int numberOfSimulation = correlationPerSitePerPos_Simulations[0][0].size();

	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			int selectedSite = site_A;	//??
			for (int site_B = 0; site_B <numOfSites_B; ++site_B){
				MDOUBLE rank = 0; //numberOfSimulation
				MDOUBLE correlVal = correlationPerSitePerPos[corIndex][selectedSite][site_B];
				for (int pos = 0; pos<numberOfSimulation; ++pos){
					MDOUBLE correlSim = correlationPerSitePerPos_Simulations[corIndex][selectedSite][pos];
					if((correlVal>correlSim && correlVal>=0 ) || (correlVal<correlSim && correlVal<0 ))
						rank++; // --
				}
				//MDOUBLE pVal = double((rank+1)/numberOfSimulation);
				correlationPerSitePerPos_Pval[corIndex][selectedSite][site_B] = rank; // pVal
			}
		}
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}

/********************************************************************************************
fill correlationPerSitePerPos_Pval, 
Compute pVal for correlationPerSitePerPos, taking into account expChanges_PosXY
*********************************************************************************************/
void computeCorrelations::computedCorrelationsPValBasedOnSimulatedDataCoMap(VVVdouble& correlationPerSitePerPosReal,VVVVdouble& expChanges_PosXYReal, VVVdouble& correlationPerSitePerPos_Pval){
	LOGnOUT(4,<<endl<<"computedCorrelationsRankBasedOnSimulatedDataCoMap..."<<endl);
	time_t t1,t2;
	time(&t1);

	MDOUBLE theWfactor = 5;
	int numberOfcorrelationVec = correlationPerSitePerPosReal.size();
	int numOfSites_A = correlationPerSitePerPosReal[0].size();
	int numOfSites_B = correlationPerSitePerPosReal[0][0].size();

	VVVdouble map_PosXY;
	computeRateValPerPos(expChanges_PosXYReal,map_PosXY);

	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		LOGnOUT(4,<<"  ***   corIndex="<<corIndex<<endl);
		MDOUBLE Nmin_min = *_NminSortedSim[corIndex].begin();
		MDOUBLE Nmin_max =  *(_NminSortedSim[corIndex].end()-1);
		MDOUBLE w_range = (Nmin_max-Nmin_min)/theWfactor;

		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			if((site_A)%100==0)
				cout<<"*";
			for (int site_B = site_A; site_B <numOfSites_B; ++site_B){
				if(site_A == site_B){
					correlationPerSitePerPos_Pval[corIndex][site_A][site_B] = 0;
					continue;
				}
				MDOUBLE Corr_obs = correlationPerSitePerPosReal[corIndex][site_A][site_B]; // Real correlation from Input variable
				int typeIndex = corIndex%_EventTypes.size();
				string type =  _EventTypes[typeIndex];
				int from = _EventTypesFromTo[type]["from"];
				int to = _EventTypesFromTo[type]["to"];
				MDOUBLE Nmin_obs = min(map_PosXY[site_A][from][to],map_PosXY[site_B][from][to]); // Real Nmin from Input variable

				MDOUBLE Nmin_lower = Nmin_obs-w_range/2;
				MDOUBLE Nmin_upper = Nmin_obs+w_range/2;

				pair<vector<double>::iterator,vector<double>::iterator> bounds;
				bounds = equal_range (_NminSortedSim[corIndex].begin(), _NminSortedSim[corIndex].end(), Nmin_lower);
				int Nmin_startIndex = int(bounds.first - _NminSortedSim[corIndex].begin());
				bounds = equal_range (_NminSortedSim[corIndex].begin(), _NminSortedSim[corIndex].end(), Nmin_upper);
				int Nmin_endIndex = int(bounds.second - _NminSortedSim[corIndex].begin());
				//cout <<Nmin_obs<< " is rang is " << Nmin_startIndex<< " and " << Nmin_endIndex << endl;

				int NumberOfSimulationsInRange = Nmin_endIndex-Nmin_startIndex+1;				
				
				//for (int i = Nmin_startIndex; i <    Nmin_endIndex; ++i){
				//	CorrelationsSubSet.push_back(_pairWiseCorrelationsAndNminSim[corIndex][1][i]); // simulations based data
				//}
				//sort(CorrelationsSubSet.begin(),CorrelationsSubSet.end());
				//bounds = equal_range (CorrelationsSubSet.begin(), CorrelationsSubSet.end(), Corr_obs);
				//int NumberOfSimulationPointsGreaterOrEqToCorr = NumberOfSimulationsInRange-int(bounds.first - CorrelationsSubSet.begin());
				//cout <<"corr= "<<Corr_obs<<" is ranked " << NumberOfSimulationPointsGreaterOrEqToCorr<< " out of " << NumberOfSimulationsInRange << endl;
				
				int NumberOfSimulationPointsGreaterOrEqToCorr =0;
				for (int i = 0; i<NumberOfSimulationsInRange; ++i){
					MDOUBLE correlSim = _pairWiseCorrelationsAndNminSim[corIndex][0][i+Nmin_startIndex];					
					if(Corr_obs<=correlSim && Corr_obs>=0 )
						NumberOfSimulationPointsGreaterOrEqToCorr++; // --
					if(gainLossOptions::_isConsiderNegativeCorrelations && Corr_obs<0 && Corr_obs>=correlSim )
						NumberOfSimulationPointsGreaterOrEqToCorr++;
				}
				
				MDOUBLE pVal = double(NumberOfSimulationPointsGreaterOrEqToCorr+1)/(NumberOfSimulationsInRange+1);
				//cout << "pVal="<<pVal<<endl;
				correlationPerSitePerPos_Pval[corIndex][site_A][site_B] = pVal;
				correlationPerSitePerPos_Pval[corIndex][site_B][site_A] = pVal;
			}
		}
		cout<<"\n";
	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
}



/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::produceSymeticMatrix(VVVdouble& correlationPerSitePerPos, bool isMin){
	LOGnOUT(4,<<endl<<"produceSymeticMatrix..."<<endl);
	int numberOfcorrelationVec = correlationPerSitePerPos.size();
	int numOfSites_A = correlationPerSitePerPos[0].size();
	int numOfSites_B = correlationPerSitePerPos[0][0].size();
	if(!numOfSites_A == numOfSites_B){
		LOGnOUT(6, <<"WARN dim not equal in produceSymeticMatrix "<<numOfSites_A<<" vs "<<numOfSites_B<<endl);
		return;
	}
	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			int selectedSite = site_A;	//??
			for (int site_B = 0; site_B <numOfSites_B; ++site_B){
				MDOUBLE minVal = correlationPerSitePerPos[corIndex][selectedSite][site_B];
				if(correlationPerSitePerPos[corIndex][site_B][selectedSite]<minVal)
					correlationPerSitePerPos[corIndex][selectedSite][site_B] = correlationPerSitePerPos[corIndex][site_B][selectedSite];
			}
		}
	}
}

/********************************************************************************************
PrintExpPerPosPerBranchMatrix (CoMap input)
NOTE!!! this version only consist of gain or loss values
Alternatively, (1) abs(gain+loss) (2) gain-loss (3) separate gain and loss matrices
*********************************************************************************************/
void computeCorrelations::printComputedCorrelationsData(const bool isNormalizeForBranch, const bool correlationForZscore
													,map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData, Vdouble& T_BH, bool isPairsAboveBH)
{
		LOGnOUT(4,<<endl<<"print Correlation data All significant sites..."<<endl);		
		int precisionCorr = 8;
		string pairsAboveBH = "";
		if(isPairsAboveBH)
			pairsAboveBH = ".pairsAboveBH";

		string corrSigSites = _outDir + "//" + "significantCorrelations.isNormForBr."+int2string(isNormalizeForBranch)+pairsAboveBH + ".txt";
		ofstream corrSigStream(corrSigSites.c_str());
		corrSigStream.precision(precisionCorr);

		//  _correlationsData["i"]["j"]["type"]["R" / "pVal" / "qVal" / "Nmin"]
		typedef map<int,map<int, map<string, map<string, MDOUBLE> > > >::iterator it_A;
		typedef map<int, map<string, map<string, MDOUBLE> > >::iterator it_B;
		typedef map<string, map<string, MDOUBLE> >::iterator it_CorrT;
		typedef map<string, MDOUBLE>::iterator it_valT;

		it_A  it1 = correlationsData.begin();	// COG A
		it_B it2 = it1->second.begin();			// COG B
		it_CorrT it3 = it2->second.begin();		// corrType
		it_valT it4 = it3->second.begin();		// valType, val (["R" / "pVal" / "qVal" / "Nmin"])		

		map<int, map<int,bool> > isPairWithSignificantPValAfterBH;		
		
		//if(!isPairsAboveBH){			
			for(it1 = correlationsData.begin(); it1 != correlationsData.end(); it1++) {				
				for(it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
					if( gainLossOptions::_isAllCorrTypeReqruiedToBeSignificant)
						isPairWithSignificantPValAfterBH[it1->first][it2->first] = true;
					else
						isPairWithSignificantPValAfterBH[it1->first][it2->first] = false;
					
					for(it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {						
						for(it4 = it3->second.begin(); it4 != it3->second.end(); it4++) {
							if( gainLossOptions::_isAllCorrTypeReqruiedToBeSignificant && it4->first == "pVal" && it4->second > T_BH[ string2double(it3->first)])
								isPairWithSignificantPValAfterBH[it1->first][it2->first] = false; // sufficient that one corType results with pVal>BH[corType] not to print
							else if (! gainLossOptions::_isAllCorrTypeReqruiedToBeSignificant && it4->first == "pVal" && it4->second<= T_BH[string2double(it3->first)])
								isPairWithSignificantPValAfterBH[it1->first][it2->first] = true; // sufficient that one corType results with pVal<=BH[corType] to print								
						}
						
					}					
				}				
			}
		//}



		// Reset, before printing Header
		it1 = correlationsData.begin();		
		it2 = it1->second.begin();
		
		// print Header
		corrSigStream<<"posA"<<"\t"<<"posB"<<"\t";
		for(it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
			for(it4 = it3->second.begin(); it4 != it3->second.end(); it4++) { // iterate over all valTypes (["R" / "pVal" / "qVal" / "Nmin"])
				corrSigStream<<it3->first<<"_"<<it4->first<<"\t";			// the combination results with e.g., 0_R	0_pVal	1_R	1_pVal
			}
		}		
		corrSigStream<<"\n";
		
		// print pair-specific computations
		for(it1 = correlationsData.begin(); it1 != correlationsData.end(); it1++) {			
			for(it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
				if(/*isPairsAboveBH ||*/ isPairWithSignificantPValAfterBH[it1->first][it2->first])
					corrSigStream<<it1->first+1<<"\t"<<it2->first+1<<"\t";
				for(it3 = it2->second.begin(); it3 != it2->second.end(); it3++) {
					for(it4 = it3->second.begin(); it4 != it3->second.end(); it4++) {
						if(/*isPairsAboveBH || */isPairWithSignificantPValAfterBH[it1->first][it2->first])
							corrSigStream<<it4->second<<"\t";
					}					
				}
				if(/*isPairsAboveBH ||*/ isPairWithSignificantPValAfterBH[it1->first][it2->first])
					corrSigStream<<"\n";
			}			
		}
		corrSigStream.close();
}



/********************************************************************************************
PrintExpPerPosPerBranchMatrix (CoMap input)
NOTE!!! this version only consist of gain or loss values
Alternatively, (1) abs(gain+loss) (2) gain-loss (3) separate gain and loss matrices
*********************************************************************************************/
void computeCorrelations::printComputedCorrelations(const Vint& selectedPositions,const Vint& evolvingSites, const bool isNormalizeForBranch, const bool correlationForZscore, VVVdouble* correlationsVec, string* valType)
{	
	// OLD version
	
	bool isOldAllAgainstAllVersion = false;
	bool isTransform = false;
	bool isMinForPrint = true;
	bool isPearson = false;
	int precisionCorr = 8;
	MDOUBLE minForPrint = 0.1; // max =1

	string pVal = "";
	if(valType)
	 pVal = *valType;
	VVVdouble correlationsVec2print;
	if(correlationsVec){
		correlationsVec2print = *correlationsVec;
		LOGnOUT(4, <<"Print correlation for external data"<<endl);
	}
	else
		correlationsVec2print = _correlationsPerSitePerPosVec;

	int numOfpositions = correlationsVec2print[0][0].size(); // assume all correlation vectors the same size
	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node

	//// Mapping vectors
	LOGnOUT(6, <<"Copy events vectors"<<endl);
	//////////////////////////////////////////////////////////////////////////	 
	if(!gainLossOptions::_printComputedCorrelationsAllSites){
		LOGnOUT(4,<<"print Correlations selected sites..."<<endl);
		for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
			int selectedSite = selectedPositions[selectedSiteIndex];
			Vdouble MeansVal(_isPearson.size()*_EventTypes.size());
			Vdouble SdVal(_isPearson.size()*_EventTypes.size());

			// for each selectedSite a new file is created
			LOGnOUT(4, <<"Correlations with site="<<selectedSite<<" With NormalizeForBranch "<<isNormalizeForBranch<<" With correlationForZscore "<<correlationForZscore<<endl);
			string corrPerSite = _outDir + "//" + "selectedCorr.Site"+ int2string(selectedSite+1)+ ".isNormForBr."+int2string(isNormalizeForBranch)+pVal+/*+ ".isCorrForZ."+int2string(correlationForZscore)*/+ ".txt";

			ofstream corrPerSiteStream(corrPerSite.c_str());
			corrPerSiteStream.precision(precisionCorr);
			corrPerSiteStream<<"# "<<selectedSite+1<<"\n";
			int vecIndex=0;
			for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
				int typeIndex=0;
				for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
					LOGnOUT(6, <<"Compute correl isPearson="<<*it<<" with type="<<*evnt<<endl);
					MeansVal[vecIndex] = computeAverage(correlationsVec2print[vecIndex][selectedSiteIndex]);
					SdVal[vecIndex] = computeStd(correlationsVec2print[vecIndex][selectedSiteIndex]);
					corrPerSiteStream<<"# Correlation isSpearman="<<*it<<" with type="<<*evnt<<" Mean="<<MeansVal[vecIndex]<<" Sd="<<SdVal[vecIndex]<<"\n";
					typeIndex++;
					vecIndex++;
				}
			}
			corrPerSiteStream<<"pos";
			vecIndex=0;
			for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
				for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
					corrPerSiteStream<<"\t"<<*evnt+int2string(*it);						
					vecIndex++;
				}
			}
			corrPerSiteStream<<"\n";
			for (int posIndex = 0; posIndex<numOfpositions; ++posIndex){
				int evolvingSite = evolvingSites[posIndex];
				if(selectedSite == evolvingSite)	// since selectedSite starts from 1
					continue;
				bool isPosOneOfSelectedSites = false;
				if(gainLossOptions::_isIgnoreCorrelationAmongSelectedSites){
					for (int selectedSiteI = 0; selectedSiteI <selectedPositions.size(); ++selectedSiteI){
						int selectedS = selectedPositions[selectedSiteI];
						if(selectedS == evolvingSite){
							isPosOneOfSelectedSites = true;
							continue;
						}				
					}
					if(isPosOneOfSelectedSites)
						continue;
				}
				corrPerSiteStream<<evolvingSite+1;
				int vecIndex=0;
				for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
					for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
						corrPerSiteStream<<"\t"<<correlationsVec2print[vecIndex][selectedSiteIndex][posIndex];						
						vecIndex++;
					}
				}
				corrPerSiteStream<<"\n";
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////	All-against-all different format
	else if(isOldAllAgainstAllVersion){
		LOGnOUT(4,<<endl<<"print Correlations All sites (old version)..."<<endl);
		string corrAllSites = _outDir + "//" + "allCorrelations.isNormForBr."+int2string(isNormalizeForBranch)+pVal+/* ".isCorrForZ."+int2string(correlationForZscore)+*/ ".txt";
		ofstream corrAllStream(corrAllSites.c_str());
		corrAllStream.precision(precisionCorr);
		corrAllStream<<"#COGA"<<"\t"<<"COGB"<<"\t"<<"posGainGain"<<"\t"<<"posLossLoss"<<"\t"<<"negGainGain"<<"\t"<<"negLossLoss"<<"\n";
		for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
			int selectedSite = selectedPositions[selectedSiteIndex];

			MDOUBLE meanCorrGainGain = computeAverage(_correlationsPerSitePerPosVec[0][selectedSiteIndex]);
			MDOUBLE stdCorrGainGain = computeStd(_correlationsPerSitePerPosVec[0][selectedSiteIndex]);
			MDOUBLE meanCorrLossLoss = computeAverage(_correlationsPerSitePerPosVec[1][selectedSiteIndex]);
			MDOUBLE stdCorrLossLoss = computeStd(_correlationsPerSitePerPosVec[1][selectedSiteIndex]);

			for (int posIndex = 0; posIndex<numOfpositions; ++posIndex){
				int evolvingSite = evolvingSites[posIndex];
				if(selectedSite == evolvingSite)
					continue;
				MDOUBLE correlationGainGain = _correlationsPerSitePerPosVec[0][selectedSiteIndex][posIndex];
				MDOUBLE correlationLossLoss = _correlationsPerSitePerPosVec[1][selectedSiteIndex][posIndex];

				if(correlationForZscore){
					correlationGainGain = (correlationGainGain - meanCorrGainGain)/stdCorrGainGain;
					correlationLossLoss = (correlationLossLoss - meanCorrLossLoss)/stdCorrLossLoss;
				}
				if(isMinForPrint && max(abs(correlationGainGain),abs(correlationLossLoss))<minForPrint)
					continue;
				MDOUBLE posCorrelationGainGain = (correlationGainGain >=0) ? correlationGainGain*1000-1 : 0;
				MDOUBLE negCorrelationGainGain = (correlationGainGain < 0) ? correlationGainGain*1000-1 : 0;
				MDOUBLE posCorrelationLossLoss = (correlationLossLoss >=0) ? correlationLossLoss*1000-1 : 0;
				MDOUBLE negCorrelationLossLoss = (correlationLossLoss < 0) ? correlationLossLoss*1000-1 : 0;
				if(isTransform){
					posCorrelationGainGain = pow(posCorrelationGainGain/10,2)/10;
					negCorrelationGainGain = pow(negCorrelationGainGain/10,2)/10;
					posCorrelationLossLoss = pow(posCorrelationLossLoss/10,2)/10;
					negCorrelationLossLoss = pow(negCorrelationLossLoss/10,2)/10;
				}
				corrAllStream<<selectedSiteIndex+1<<"\t"<<evolvingSite+1<<"\t"<<(int)posCorrelationGainGain<<"\t"<<(int)posCorrelationLossLoss<<"\t"<<(int)negCorrelationGainGain<<"\t"<<(int)negCorrelationLossLoss<<"\n";
			}
		}
	}
	else{
		LOGnOUT(4,<<"print Correlations All sites ..."<<endl);
		string corrAllSites = _outDir + "//" + "allCorrelations.isNormForBr."+int2string(isNormalizeForBranch)+pVal+ /*".isCorrForZ."+int2string(correlationForZscore)+*/ ".txt";
		ofstream corrAllStream(corrAllSites.c_str());
		corrAllStream.precision(precisionCorr);
		corrAllStream<<"siteA"<<"\t"<<"siteB";
		int vecIndex=0;
		for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
			for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
				corrAllStream<<"\t"<<*evnt+int2string(*it);						
				vecIndex++;
			}
		}
		corrAllStream<<"\n";

		for (int selectedSiteIndex = 0; selectedSiteIndex <selectedPositions.size(); ++selectedSiteIndex){
			int selectedSite = selectedPositions[selectedSiteIndex];
			for (int posIndex = 0; posIndex<numOfpositions; ++posIndex){
				int evolvingSite = evolvingSites[posIndex];
				if(selectedSite == evolvingSite)
					continue;
				corrAllStream<<selectedSite+1<<"\t"<<evolvingSite+1;
				int vecIndex=0;
				for (vector<bool>::iterator it=_isPearson.begin() ; it < _isPearson.end(); it++ ){
					for (vector<string>::iterator evnt=_EventTypes.begin() ; evnt < _EventTypes.end(); evnt++ ){ // could be done with int
						corrAllStream<<"\t"<<correlationsVec2print[vecIndex][selectedSiteIndex][posIndex];						
						vecIndex++;
					}
				}
				corrAllStream<<"\n";
			}
		}
	}
}


/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::fillCorrPerSelectedSites(Vdouble& correlationPerPos,VVdouble& expEventsPerPosPerBranch,VVdouble& expEventsPerPosPerBranch_B,const int selectedSite, const bool isPearson){
	int numOfpositions = expEventsPerPosPerBranch_B.size();
	//correlationPerPos.resize(numOfpositions);
	
	
	for (int pos = 0; pos <numOfpositions; ++pos){
		MDOUBLE correlation = 0;
		if(isMinEQMaxInVector(expEventsPerPosPerBranch[selectedSite]) || isMinEQMaxInVector(expEventsPerPosPerBranch_B[pos]))
			correlationPerPos[pos]=-99; // can't compute correlation
		else{
			if(isPearson)
				correlation = calcPearsonCorrelation(expEventsPerPosPerBranch[selectedSite], expEventsPerPosPerBranch_B[pos]);
			else{
				//correlation = calcRankCorrelation(expEventsPerPosPerBranch[selectedSite], expEventsPerPosPerBranch_B[pos]); // seems to be problematic, diffrent results from R, Matlab
				correlation = calcRankCorrelation2(expEventsPerPosPerBranch[selectedSite], expEventsPerPosPerBranch_B[pos]);				
			}
			correlationPerPos[pos]=correlation;
		}
	}		
}



/********************************************************************************************
fill expEventsPerPosPerBranch
*********************************************************************************************/
void computeCorrelations::fillMapValPerPosPerBranch(VVdouble& expEventsPerPosPerBranch,const string type, VVVVdouble& expChanges_PosNodeXY
												   ,const bool isNormalizeForBranch, MDOUBLE* cutOff_p){

	
	int numOfpositions = expChanges_PosNodeXY.size();
	int numOfbranches = _tr.getNodesNum()-1; // was -1, minus the root node

	int from = _EventTypesFromTo[type]["from"];
	int to = _EventTypesFromTo[type]["to"];
	expEventsPerPosPerBranch.resize(numOfpositions);
	treeIterTopDownConst tIt(_tr);
	for (int pos = 0; pos <numOfpositions; ++pos){
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) 
		{
			if(mynode->isRoot())
				continue;
			MDOUBLE val = 0;
			MDOUBLE normalizationFactor = 1.0;
			if(isNormalizeForBranch){				
				if(gainLossOptions::_isNormalizeByExpectationPerBranch){
					if(_expChanges_NodeXY.size()==0)
						sumExpectationPerBranch(expChanges_PosNodeXY, _expChanges_NodeXY); // filled once for both 0->1 and 1->0					
					normalizationFactor =  _expChanges_NodeXY[mynode->id()][from][to]/numOfbranches; // mynode->dis2father()
				}else
					normalizationFactor = mynode->dis2father();				
			}
			val = (expChanges_PosNodeXY[pos][mynode->id()][from][to] ) / normalizationFactor;
			if(cutOff_p){
				if(val>= *cutOff_p)
					expEventsPerPosPerBranch[pos].push_back(1);
				else
					expEventsPerPosPerBranch[pos].push_back(0);
			}
			else
				expEventsPerPosPerBranch[pos].push_back(val);			
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
void computeCorrelations::sumExpectationPerBranch(VVVVdouble& expChanges_PosNodeXY, VVVdouble& map_NodeXY){
	int numOfPositions = expChanges_PosNodeXY.size();
	int numOfBranches = expChanges_PosNodeXY[0].size();
	int AlphSize = expChanges_PosNodeXY[0][0].size(); // =2

	treeIterTopDownConst tIt(_tr);
	resizeVVV(numOfBranches,AlphSize,AlphSize,map_NodeXY);
	for (int pos = 0; pos <numOfPositions; ++pos){
		//int i=0;
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()){			
		//for(int i=0;i<numOfBranches;++i){
			for(int j=0;j<AlphSize;++j){
				for(int k=0;k<AlphSize;++k){
					map_NodeXY[mynode->id()][j][k] += expChanges_PosNodeXY[pos][mynode->id()][j][k];
				}
			}			
			//cout<<i<<" "<<mynode->id()<<endl; // DEBUG
			//++i;
		}
	}
}

/********************************************************************************************
//Compute p-values of each statistic: P1
//, P2
//, P3
//, ??? , PN
//• Sort these: P(1)
//? P(2)
//? P(3)
//? ??? ? P(N)
//{subscript
//()
//? sorted}
//• For k =1..N, q(k)
//= minm ? k
//[ N?P(m)
///m]
//? Easily computed from sorted p-values by looping
//downwards from k= N to k =1
*********************************************************************************************/
VVVdouble computeCorrelations::pVals2qVals(VVVdouble& pValVec, map<int, map<int, map<string,  map<string, MDOUBLE > > > >& correlationsData
	, vector<vector<bool> >& isComputePairWithRateAboveNim, Vdouble& T_BH, Vint& selectedSites, Vint& evolvingSites)
{
	LOGnOUT(4,<<endl<<"pVals2qVals..."<<endl);
	time_t t1,t2;
	time(&t1);
	
	VVVdouble qValsVec;
	if(gainLossOptions::_isComputeQVals)		
		qValsVec = pValVec; // instead of re-size, not required if not computed

	int numberOfcorrelationVec = pValVec.size();
	int numOfSites_A = pValVec[0].size();
	int numOfSites_B = pValVec[0][0].size();

	typedef map<int,map<int, map<string, map<string, MDOUBLE> > > >::iterator it_A;
	typedef map<int, map<string, map<string, MDOUBLE> > >::iterator it_B;

	it_A it_siteA = correlationsData.begin();
	it_B it_siteB = it_siteA->second.begin();

	for (int corIndex = 0; corIndex <numberOfcorrelationVec; ++corIndex){
		LOGnOUT(4,<<"  ***   corIndex="<<corIndex<<endl);

		Vdouble pVals;
		Vdouble qVals;

		// get pVals
		LOGnOUT(6,<<"get pVals..."<<endl);
		
		for (int site_A = 0; site_A <numOfSites_A; ++site_A){
			int site_A_original = selectedSites[site_A];
			for (int site_B = site_A; site_B <numOfSites_B; ++site_B){
				int site_B_original = evolvingSites[site_B];
				if(site_A_original == site_B_original){
					continue;
				}
				MDOUBLE pVal = 1;
				if(isComputePairWithRateAboveNim[site_A][site_B]){
					if(gainLossOptions::_selectedSitesForCorrelation!="" ){
					// consider only pairs with min Rate
					//if(correlationsData[site_A_original][site_B][int2string(corIndex)]["Rate"]>Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair") ){ //TEMP
						pVal = correlationsData[site_A_original][site_B_original][int2string(corIndex)]["pVal"]; // 
						//if(pValVec[corIndex][site_A][site_B] != pVal)
						//	cout<<"ERRRRRRR diff pval\n";
						//if(correlationsData[site_A_original][site_B][int2string(corIndex)]["Rate"]>Parameters::getFloat("_minExpThresholdForPValComputationForCorrelatingPair") && pVal > 1.99)
						//	cout<<"ERRRRRRR diff pval\n";
					}
					else
						pVal = pValVec[corIndex][site_A][site_B]; // Real correlation from Input variable

					pVals.push_back(pVal);
				}
				//if(!(pVal > 1)) // pair is with Nmin below T, and ignored, since it's removed from both simulations and real data no need to correct for this hypothesis				
				//	pVals.push_back(pVal);					
									
			}
		}
		// sort pVal
		vector< vecElem<MDOUBLE> > orderVecPVal;
		orderVec(pVals, orderVecPVal);
		qVals.resize(pVals.size(),1);
		
		sort(pVals.begin(),pVals.end()); // faster than using the "getValue"
		pair<vector<double>::iterator,vector<double>::iterator> bounds;
		
		float pVals2checkBeforeFDR [] = {gainLossOptions::_pValueCutOffForBootStrap, 0.05, 0.01, 0.005, 0.001, 0.0001};
		
		int lastIndexWithPVal2check;
		for (int i=0; i<6; ++i){
			bounds = equal_range (pVals.begin(), pVals.end(), pVals2checkBeforeFDR[i]);
			if(i==0)
				lastIndexWithPVal2check = int(bounds.second - pVals.begin());
			int lastIndexWithPValBiggerThanThreshold = int(bounds.second - pVals.begin());
			LOGnOUT(4,<<"Before FDR correction there are "<<lastIndexWithPValBiggerThanThreshold<<" pairs with significant pVal="<< pVals2checkBeforeFDR[i]<<endl);
		}		
		
		LOGnOUT(4,<<"Compute BH threshold for number of multiple tests="<<pVals.size()<<  " ..."<<endl);		
		T_BH[corIndex] = 0;
		T_BH[corIndex] = computeFDRthreshold(pVals, gainLossOptions::_pValueCutOffForBootStrap, true);


		//for (int i=0; i<pVals.size(); ++i){
		//	MDOUBLE correctedVal = (double)(i+1)/(double)indexAll *pValcutOff;
		//	if( pVals[i] <= correctedVal){
		//		T_BH[corIndex] = pVals[i];
		//	}			
		//}
		bounds = equal_range (pVals.begin(), pVals.end(),T_BH[corIndex]);		
		if(T_BH[corIndex] > 0.0){
			LOGnOUT(4,<<"For FDR level of "<<gainLossOptions::_pValueCutOffForBootStrap<<" BH threshold="<<T_BH[corIndex]<<" with "<<int(bounds.first - pVals.begin())<<" "<<int(bounds.second - pVals.begin())<<" significant values"<<endl);}
		else{
			LOGnOUT(4,<<"For FDR level of "<<gainLossOptions::_pValueCutOffForBootStrap<<" BH threshold="<<T_BH[corIndex]<<" with no significant values"<<endl<<endl);}

		// additional BH thresholds
		float additionalFDRlevels [] = {0.1, 0.05, 0.01, 0.001};
		for (int i=0; i< 4 ; ++i){ // must be length of additionalFDRlevels
			if(gainLossOptions::_pValueCutOffForBootStrap == additionalFDRlevels[i])
				continue;
			MDOUBLE BH = computeFDRthreshold(pVals,additionalFDRlevels[i], true);
			LOGnOUT(4,<<"For FDR level of "<<additionalFDRlevels[i]<<" BH threshold is "<<BH<<endl);
		}

		// compute q-vals
		VVVdouble qValsVec;
		if(gainLossOptions::_isComputeQVals){
			// produce qVals by FDR, assume the pVal vector is sorted P_1<=P_2<=...<=P_N			 
			LOGnOUT(4,<<"Compute FDR correction - get qVals..."<<endl);
			for (int k=1; k<=lastIndexWithPVal2check; ++k){
				if(k%1000==0)
					cout<<"*";
				int m = k;
				MDOUBLE pVal = pVals[k-1];
				MDOUBLE qVal = pVal;
				//cout<<"pVal "<<k<<" "<<orderVecPVal[m-1].getValue()<<endl; //DEB
				//MDOUBLE qVal = (double)indexAll*pVal/(double)m; // only init
				if(pVal < gainLossOptions::_pValueCutOffForBootStrap && qVal < gainLossOptions::_pValueCutOffForBootStrap){ // since pVals are sorted, if last qVal computation yielded >0.05, no need to compute
					qVal = 1; // init, not corrected
					for (m=k; m<= lastIndexWithPVal2check; ++m){
						MDOUBLE pValtemp = pVals[m-1];
						MDOUBLE qValtemp = (double)pVals.size()*pValtemp/(double)m;
						if(qValtemp < qVal)
							qVal = qValtemp;
					}
				}else{
					break;
				}
				//cout<<"pVal="<<pVal<<" and qVal="<<qVal<<"\n"; //DEB
				qVals[orderVecPVal[k-1].getPlace()] = qVal; // fill values, related to original order of pVals
			}
			cout<<"\n";

			// assign qVals
			LOGnOUT(4,<<"Assign qVals..."<<endl);
			int ind = 0; 
			for (int site_A = 0; site_A <numOfSites_A; ++site_A){
				for (int site_B = site_A; site_B <numOfSites_B; ++site_B){
					if(site_A == site_B){
						continue;
					}
					MDOUBLE qVal = qVals[ind]; // this works only if qVals order is the same as the original pVal
					qValsVec[corIndex][site_A][site_B] = qVal;
					qValsVec[corIndex][site_B][site_A] = qVal;
					//map<string, Vdouble>::iterator iterTerm = _totalTerminals.find(nodeName);

					it_A iterA = correlationsData.find(site_A);
					it_B iterB = correlationsData[site_A].find(site_B);

					if (!(iterA==correlationsData.end()) && !(iterB==correlationsData[site_A].end())){				
						//cout<<site_A<<" "<<site_B<<"\n";
						correlationsData[site_A][site_B][int2string(corIndex)]["qVal"] = qVal;
					}
					ind++;
				}
			}
		}

	}
	time(&t2);
	LOGnOUT(4,<<"TIME = "<<(t2-t1)/60.0<<" minutes"<<endl<<endl);
	return qValsVec;
}
