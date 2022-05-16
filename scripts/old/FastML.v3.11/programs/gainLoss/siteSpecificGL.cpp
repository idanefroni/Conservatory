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

#include "siteSpecificGL.h"
#include "definitions.h"
#include "numRec.h"
#include "matrixUtils.h"
#include "seqContainerTreeMap.h"
#include "gainLossUtils.h"
#include "gainLossModel.h"
#include "gainLossOptions.h"



// THE BAYESIAN EB_EXP PART OF gain and loss ESTIMATION. //

/*************************************
This function computes the expectation of 
the posterior gain and loss distribution for a specific site
as well as the confidence interval 
*************************************/
// per all sites computation
void computeEB_EXP_siteSpecificGL(Vdouble & GainLossV,
									Vdouble & stdV,
									Vdouble & lowerBoundV,
									Vdouble & upperBoundV,
									VVdouble & posteriorsV,
									const sequenceContainer& sc,
									const vector<vector<stochasticProcess*> >& spVVec,
									const tree& tr,
									const distribution * gainDist,
									const distribution * lossDist,
									const distribution * distPrim,
									const MDOUBLE alphaConf,
									VVVdouble & postProbPerSpPerCatPerPos,	//2 fill (*postProbPerSpPerCatPerPos)[sp][pos]
									unObservableData* unObservableData_p)
{
	LOG(5,<<"Calculating posterior and expectation of posterior values for all sites"<<endl);
	int seqLen = sc.seqLen();
	GainLossV.resize(seqLen);
	stdV.resize(seqLen);
	lowerBoundV.resize(seqLen);
	upperBoundV.resize(seqLen);	
	int numOfSPs = gainDist->categories()*lossDist->categories();
	resizeMatrix(posteriorsV,seqLen,numOfSPs);
	//computePijGam cpg;
	//cpg._V.resize(numOfSPs);
	//for (int i=0; i < numOfSPs; ++i) {
	//	int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
	//	int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
	//	cpg._V[i].fillPij(tr,*spVVec[gainIndex][lossIndex]);
	//}
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		computeEB_EXP_siteSpecificGL(pos, sc, spVVec, tr, gainDist,lossDist,distPrim,posteriorsV[pos],	//cpg
			GainLossV[pos], stdV[pos], lowerBoundV[pos], upperBoundV[pos], alphaConf, postProbPerSpPerCatPerPos,unObservableData_p);
	}
}


/********************************************************************************************
*********************************************************************************************/
void computeEB_EXP_siteSpecificGL(int pos,
									const sequenceContainer& sc,
									const vector<vector<stochasticProcess*> >& spVVec,
									//const computePijGam& cpg,
									const tree &tr,
									const distribution * gainDist,
									const distribution * lossDist,
									const distribution * distPrim,
									Vdouble & posteriorV,
									MDOUBLE& GainLossExpectation,
									MDOUBLE & stdGainLoss,
									MDOUBLE & lowerConf,
									MDOUBLE & upperConf,
									const MDOUBLE alphaConf,
									VVVdouble & postProbPerSpPerCatPerPos,	//2 fill (*postProbPerSpPerCatPerPos)[sp][pos]
									unObservableData* unObservableData_p) // alpha of 0.05 is considered 0.025 for each side.
{
	bool isLpostPerSpPerCatComputed =false;
	if(postProbPerSpPerCatPerPos[0][0][pos]>0)
		isLpostPerSpPerCatComputed =true;


	// here we compute the posterior P(r|data)
	int numOfRateCat = (*spVVec[0][0]).categories();	// ver2
	int numOfSPs = gainDist->categories()*lossDist->categories();

	posteriorV.resize(distPrim->categories(),0.0);
	// ver2
	VVdoubleRep PosteriorVVRateCat;
	resizeMatrix(PosteriorVVRateCat,numOfSPs,numOfRateCat);

	doubleRep dRepTotalLikelihood(0.0);// temporary dblRep for total likelihood

	for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
		int gainIndex =fromIndex2gainIndex(spIndex,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(spIndex,gainDist->categories(),lossDist->categories());
		
		//int primIndex;
		//if(distPrim == gainDist)
		//	primIndex = gainIndex;
		//else
		//	primIndex = lossIndex;
		
		computePijGam pi;
		pi.fillPij(tr,*spVVec[gainIndex][lossIndex]);
		
		// ver1 - no rate dist in rate computation
		//dblRepPosteriorV[primIndex] += likelihoodComputation::getLofPos(pos,tr,sc,pi,*spVVec[gainIndex][lossIndex])* gainDist->ratesProb(gainIndex)*lossDist->ratesProb(lossIndex);

		// ver2 - with rate dist
		for (int rateInd=0; rateInd < numOfRateCat; ++rateInd) {
			PosteriorVVRateCat[spIndex][rateInd] += likelihoodComputation::getLofPos(pos,tr,sc,pi[rateInd],*spVVec[gainIndex][lossIndex],unObservableData_p)
					* gainDist->ratesProb(gainIndex) * lossDist->ratesProb(lossIndex) * spVVec[gainIndex][lossIndex]->ratesProb(rateInd);
		}		
	}

	// here we compute sigma r * P(r | data)
	GainLossExpectation = 0.0;
	MDOUBLE sumOfSquares = 0.0; // this is the sum of squares. this will be used to compute the variance
	
	// ver1 - no rate dist in rate computation
	//for (int i=0; i < distPrim->categories(); ++i) {
	//	dblRepTotalLikelihood+=dblRepPosteriorV[i];
	//}
	//for (int j=0; j < distPrim->categories(); ++j) {
	//	dblRepPosteriorV[j]/=dblRepTotalLikelihood; // so that posteriorV is probability.
	//	if(unObservableData_p){
	//		dblRepPosteriorV[j] = dblRepPosteriorV[j]/(1- exp(unObservableData_p->getlogLforMissingData()));	// Note: each postProbCat corrected by unObs of all cat
	//	}
	//	posteriorV[j] = convert(dblRepPosteriorV[j]); // revert back to DOUBLE
	//	MDOUBLE tmp = posteriorV[j]*distPrim->rates(j);
	//	GainLossExpectation += tmp;
	//	sumOfSquares += (tmp*distPrim->rates(j));
	//}

	// ver2
	for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
		for (int i=0; i < numOfRateCat; ++i) {
			dRepTotalLikelihood+=PosteriorVVRateCat[spIndex][i];
		}
	}
	
	
	for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
		int gainIndex =fromIndex2gainIndex(spIndex,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(spIndex,gainDist->categories(),lossDist->categories());

		int primIndex;
		if(distPrim == gainDist)
			primIndex = gainIndex;
		else
			primIndex = lossIndex;

		for (int i=0; i < numOfRateCat; ++i) {
			PosteriorVVRateCat[spIndex][i]/=convert(dRepTotalLikelihood); // so that posteriorV is probability.
			posteriorV[primIndex] += convert(PosteriorVVRateCat[spIndex][i]);
			MDOUBLE tmp = convert(PosteriorVVRateCat[spIndex][i]) * distPrim->rates(primIndex) * spVVec[0][0]->rates(i); // the rateVal
			GainLossExpectation += tmp;
			sumOfSquares += (tmp * distPrim->rates(primIndex) * spVVec[0][0]->rates(i));	// ???
		}
	}
////////////////////////////////////////////////////////////////////////// ?
	if(!isLpostPerSpPerCatComputed){
		for (int spIndex=0; spIndex < numOfSPs; ++spIndex) {
			for (int rateInd=0; rateInd < numOfRateCat; ++rateInd) {
				postProbPerSpPerCatPerPos[spIndex][rateInd][pos] = convert(PosteriorVVRateCat[spIndex][rateInd]);
			}
		}
	}
	MDOUBLE variance = sumOfSquares - GainLossExpectation*GainLossExpectation; // variance
	//if (!(variance!=0))
	//	errorMsg::reportError("Error in computeEB_EXP_siteSpecificGainLoss, variance = 0");
	stdGainLoss = sqrt(variance); // standard deviation of inferred Ka/Ks

	// detecting the confidence intervals.
	MDOUBLE oneSideConfAlpha = alphaConf/2.0; // because we are computing the two tail.
	MDOUBLE cdf = 0.0; // cumulative density function.
	int k=0;
	while (k < distPrim->categories()){
		cdf += posteriorV[k];
		if (cdf >oneSideConfAlpha) {
			lowerConf = distPrim->rates(k);
			break;
		} 
		k++;
	}
	while (k < distPrim->categories()) {
		if (cdf >(1.0-oneSideConfAlpha)) {
			upperConf = distPrim->rates(k);
			break;
		}
		++k;
		cdf += posteriorV[k];
	}
	if (k==distPrim->categories()) 
		upperConf = distPrim->rates(k-1);
}



