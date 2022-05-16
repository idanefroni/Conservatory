// $Id: siteSpecificRate.cpp 11008 2012-10-16 21:54:04Z rubi $

#include "siteSpecificRate.h"
#include "numRec.h"
#include "checkcovFanctors.h"
#include "definitions.h"


/********************************************************************************************
ML - full data (1)
*********************************************************************************************/
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & likelihoodsV,
								   const sequenceContainer& sc,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE maxRate,//20.0f
								   const MDOUBLE tol){//=0.0001f;

	ratesV.resize(sc.seqLen());
	likelihoodsV.resize(sc.seqLen());
	MDOUBLE Lsum = 0.0;

	for (int pos=0; pos < sc.seqLen(); ++pos) {
		computeML_siteSpecificRate(pos,sc,sp,et,ratesV[pos],likelihoodsV[pos],maxRate,tol);
		assert(log(likelihoodsV[pos])>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(6,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}
/********************************************************************************************
ML - per Pos (1.1)
*********************************************************************************************/
// note that this places the likelihood, rather then the *log*likelihood into posL
void computeML_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE& posL,
								 const MDOUBLE maxRate,
								 const MDOUBLE tol) {
	LOG(6,<<".");
	MDOUBLE ax=0.00001f,bx=maxRate*0.25,cx=maxRate;	// MN
	posL=-brent(ax,bx,cx,Cevaluate_L_given_r(sc,et,sp,pos),tol,&bestRate);
}


/********************************************************************************************
ML - full data AttributesVecs (1)
*********************************************************************************************/
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const vector<const stochasticProcess *> & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	MDOUBLE Lsum = 0.0;
	ratesV.resize(sc.seqLen()); // the rates themselves
	likelihoodsV.resize(sc.seqLen()); // the log likelihood of each position
	
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		LOG(6,<<".");
		MDOUBLE bestR=-1.0; // tree1
		//		MDOUBLE LmaxR1=0;
		
		// getting the right tree for the specific position:
		const tree*  treeForThisPosition=NULL;
		if ((etVec.size() >0 ) && (treeAttributesVec[pos]>0)) {
			treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("tree vector is empty, or treeAttribute is empty, or treeAttribute[pos] is zero (it should be one)");
		}

		// getting the right stochastic process for the specific position:

		const stochasticProcess* spForThisPosition=NULL;

		if ((spVec.size() >0 ) && (spAttributesVec[pos]>0)) {
			spForThisPosition = spVec[ spAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("stochastic process vector is empty, or spAttributesVec is empty, or spAttribute[pos] is zero (it should be one)");
		}

		computeML_siteSpecificRate(pos,sc,*spForThisPosition,*treeForThisPosition,bestR,likelihoodsV[pos],maxRate,tol);
		ratesV[pos] = bestR;
		assert(log(likelihoodsV[pos])>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(6,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}
/********************************************************************************************
ML - AttributesVecs (1.1)
*********************************************************************************************/
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& treeAttributesVec,	//treeAttributesVec
						const vector<tree> & etVec,
						const stochasticProcess& sp,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol) {
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	return computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}
/********************************************************************************************
ML - AttributesVecs (1.1)
*********************************************************************************************/
MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,	// spAttributesVec
						const tree & et,
						const vector<const stochasticProcess* > & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	Vint treeAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	return computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}



// THE BAYESIAN EB_EXP PART OF RATE ESTIMATION. //
/********************************************************************************************
EB_EXP - full data (1)
*********************************************************************************************/
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
										Vdouble & stdV,
										Vdouble & lowerBoundV,
										Vdouble & upperBoundV,
										const sequenceContainer& sc,
										const stochasticProcess& sp,
										const tree& et,
										const MDOUBLE alphaConf,
										VVdouble* LpostPerCat,	//2 fill (*LpostPerCat)[cat][pos]
										unObservableData* unObservableData_p)
{
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());

	computePijGam cpg;
	cpg.fillPij(et,sp);
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		computeEB_EXP_siteSpecificRate(pos,sc,sp,cpg, et,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf,LpostPerCat,unObservableData_p);
		LOG(6,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(6,<<" number of sites: "<<sc.seqLen()<<endl);
}


/********************************************************************************************
EB_EXP - per Pos (1.1)
*********************************************************************************************/
void computeEB_EXP_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const computePijGam& cpg,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE & stdRate,
								 MDOUBLE & lowerConf,
								 MDOUBLE & upperConf,
								 const MDOUBLE alphaConf, // alpha of 0.05 is considered 0.025 for each side.
								 VVdouble* LpostPerCat,		//2 fill (*LpostPerCat)[cat][pos]
								 unObservableData* unObservableData_p)
{
	// here we compute P(r | data)
	VdoubleRep pGivenR(sp.categories(),0.0);
	doubleRep sum=0;
	doubleRep LofPos_givenRateCat;
	LOG(8,<<pos+1<<"\t"); //DEBUG
	for (int cat=0; cat < sp.categories(); ++cat) {
		LofPos_givenRateCat = likelihoodComputation::getLofPos(pos,et,sc,cpg[cat],sp);

// ver1 - fix likelihoodForEachCat by LforMissingDataPerCat
		//if(unObservableData_p){
		//	LofPos_givenRateCat = LofPos_givenRateCat/(1- unObservableData_p->getLforMissingDataPerCat()[cat]);
		//}
// ver2 - fix likelihoodForEachCat by LforMissingDataAll
		if(unObservableData_p){
			LofPos_givenRateCat = LofPos_givenRateCat/(1- exp(unObservableData_p->getlogLforMissingData()));
		}
		pGivenR[cat] = LofPos_givenRateCat * sp.ratesProb(cat); 
		sum+=pGivenR[cat]; 
	}
	LOG(8,<<"\n"); //DEBUG
	assert(sum!=0);
	
	// here we compute sigma r * P(r | data)
	doubleRep sumOfSquares(0.0);
	doubleRep bestRate_dblRep(0.0);
	
	LOG(6,<<"Pos "<<pos<<" content = "<<sc[0][pos]<<" ,total likelihood = "<<sum<<endl); //DEBUG
	
	for (int j=0; j < sp.categories(); ++j) {
		pGivenR[j]/=sum; // So that pGivenR is probability.
		                 // From here on we can convert it back
		                 // to MDOUBLE because it's not a very
		                 // small likelihood any more

// ver3 - fix likelihoodForEachCat after multiplied by Prob - Error??
		//if(unObservableData_p){
		//	pGivenR[j] = pGivenR[j]/(1- (unObservableData_p->getLforMissingDataPerCat())[j]) ;	// Note: each postProbCat corrected by unObs of a
		//}
		
		if (LpostPerCat){
			(*LpostPerCat)[j][pos]= convert(pGivenR[j]);
		}
		doubleRep tmp = pGivenR[j]*sp.rates(j);
		LOG(8,<<j<<"\t"<<sp.rates(j)<<"\t"<<convert(pGivenR[j])<<"\t"); //DEBUG
		bestRate_dblRep += tmp;
		sumOfSquares += (tmp*sp.rates(j));
	}

	bestRate = convert(bestRate_dblRep);
	MDOUBLE varRate = convert(sumOfSquares) - convert(bestRate*bestRate);
	MDOUBLE tolerance = 0.0001; // tolerance for variance is not very exact, and also exact computation not very important
	if (varRate<-tolerance)
		LOGnOUT(3,<<"Error in computeEB_EXP_siteSpecificRate pos="<<pos<<", varRate="<<varRate<<" (< 0) \n");
	if ((varRate<0) && (varRate>=-tolerance))
	    varRate = 0;
	stdRate = sqrt(varRate);

	// detecting the confidence intervals.
	MDOUBLE oneSideConfAlpha = alphaConf/2.0; // because we are computing the two tail.
	MDOUBLE cdf = 0.0; // cumulative density function.
	MDOUBLE lower_interval = 0;
	MDOUBLE total_interval = 0;
	int k=0;
	while (k < sp.categories()){
		cdf += convert(pGivenR[k]);
		if (cdf >oneSideConfAlpha) {
			if(k>0) {
				lowerConf = sp.rates(k-1);
				lower_interval = convert(pGivenR[k-1]);
			}
			else {
				lowerConf = 0;
				lower_interval = 0;
			}
			break;
		} 
		k++;
	}
	while (k < sp.categories()) {
		if (cdf >(1.0-oneSideConfAlpha)) {
			upperConf = sp.rates(k);
			total_interval = cdf - lower_interval;
			break;
		}
		++k;
		cdf += convert(pGivenR[k]);
	}
	if (k==sp.categories()) {
		upperConf = sp.rates(k-1);
		total_interval = 1.0 - lower_interval;
	}
	LOG(7,<<"Pos: "<<pos<<", conf_interval= "<<total_interval<<endl); 
}

// THE PROPORTIONAL BAYESIAN EB_EXP PART OF RATE ESTIMATION. //
/********************************************************************************************
EB_EXP - full data (1)
*********************************************************************************************/
void computeEB_EXP_siteSpecificRateProportional(Vdouble & ratesV,
										Vdouble & stdV,
										Vdouble & lowerBoundV,
										Vdouble & upperBoundV,
										const vector<sequenceContainer>& scVec,
										multipleStochasticProcess& msp,
										const gammaDistribution* pProportionDist,
										const tree& et,
										const MDOUBLE alphaConf,
										VVdouble* LpostPerCat)	//2 fill (*LpostPerCat)[globalRate][localRate][pos]
{
	ratesV.resize(scVec.size());
	stdV.resize(scVec.size());
	lowerBoundV.resize(scVec.size());
	upperBoundV.resize(scVec.size());

	for (int gene=0; gene < scVec.size(); ++gene) {
		computeEB_EXP_siteSpecificRateProportional(gene,scVec[gene],msp,pProportionDist,et,ratesV[gene],stdV[gene],lowerBoundV[gene],upperBoundV[gene],alphaConf,LpostPerCat);
		LOG(6,<<" rate of gene "<<gene<<" = "<<ratesV[gene]<<endl);
	}
	LOG(6,<<" number of genes: "<<scVec.size()<<endl);
}


/********************************************************************************************
EB_EXP - per Pos (1.1)
*********************************************************************************************/
void computeEB_EXP_siteSpecificRateProportional(int gene,
												const sequenceContainer& sc,
												multipleStochasticProcess& msp,
												const gammaDistribution* pProportionDist,
												const tree &et,
												MDOUBLE& bestRate,
												MDOUBLE & stdRate,
												MDOUBLE & lowerConf,
												MDOUBLE & upperConf,
												const MDOUBLE alphaConf, 
												VVdouble* LpostPerCat)	
{
	
	// here we compute P(r | data)
	VdoubleRep pGivenR(pProportionDist->categories(),0.0);
	doubleRep sum=0;
	doubleRep LofGene_givenRateCat = 0.0;
	LOG(8,<<gene+1<<"\t"); //DEBUG
	for (int cat=0; cat < pProportionDist->categories(); ++cat) {
		msp.getSp(gene)->setGlobalRate(pProportionDist->rates(cat));
		computePijGam cpg;
		cpg.fillPij(et,*msp.getSp(gene));
		for (int k=0; k < sc.seqLen(); ++k) {
			 LofGene_givenRateCat += likelihoodComputation::getLofPosProportional(k,//pos,
					et,		//const tree& 
					sc,		// sequenceContainer& sc,
					cpg,		//const computePijGam& ,
					*msp.getSp(gene)); //removed the prior of the globar rate categ cause it is multiplied below
		}
		pGivenR[cat] = LofGene_givenRateCat*pProportionDist->ratesProb(cat);
		sum+=pGivenR[cat];
	}
	LOG(8,<<"\n"); //DEBUG
	assert(sum!=0);
	
	// here we compute sigma r * P(r | data)
    doubleRep sumOfSquares(0.0);
	doubleRep bestRate_dblRep(0.0);
	
	for (int j=0; j < pProportionDist->categories(); ++j) {
		pGivenR[j]/=sum; // So that pGivenR is probability.
		                 // From here on we can convert it back
		                 // to MDOUBLE because it's not a very
		                 // small likelihood any more

		if (LpostPerCat){
			(*LpostPerCat)[j][gene]= convert(pGivenR[j]);
		}
		doubleRep tmp = pGivenR[j]*pProportionDist->rates(j);
		LOG(8,<<j<<"\t"<<pProportionDist->rates(j)<<"\t"<<convert(pGivenR[j])<<"\t"); //DEBUG
		bestRate_dblRep += tmp;
		sumOfSquares += (tmp*pProportionDist->rates(j));
	}
	
	bestRate = convert(bestRate_dblRep);
	MDOUBLE varRate = convert(sumOfSquares) - convert(bestRate*bestRate);
	MDOUBLE tolerance = 0.0001; // tolerance for variance is not very exact, and also exact computation not very important
	if (varRate<-tolerance)
		LOGnOUT(3,<<"Error in computeEB_EXP_siteSpecificRateProportional gene="<<gene<<", varRate="<<varRate<<" (< 0) \n");
	if ((varRate<0) && (varRate>=-tolerance))
	    varRate = 0;
	stdRate = sqrt(varRate);


	// detecting the confidence intervals.
	MDOUBLE oneSideConfAlpha = alphaConf/2.0; // because we are computing the two tail.
	MDOUBLE cdf = 0.0; // cumulative density function.
	MDOUBLE lower_interval = 0;
	MDOUBLE total_interval = 0;
	int k=0;
	while (k < pProportionDist->categories()){
		cdf += convert(pGivenR[k]);
		if (cdf >oneSideConfAlpha) {
			if(k>0) {
				lowerConf = pProportionDist->rates(k-1);
				lower_interval = convert(pGivenR[k-1]);
			}
			else {
				lowerConf = 0;
				lower_interval = 0;
			}
			break;
		} 
		k++;
	}
	while (k < pProportionDist->categories()) {
		if (cdf >(1.0-oneSideConfAlpha)) {
			upperConf = pProportionDist->rates(k);
			total_interval = cdf - lower_interval;
			break;
		}
		++k;
		cdf += convert(pGivenR[k]);
	}
	if (k==pProportionDist->categories()) {
		upperConf = pProportionDist->rates(k-1);
		total_interval = 1.0 - lower_interval;
	}
	LOG(7,<<"Gene: "<<gene<<", conf_interval= "<<total_interval<<endl); 
}

/********************************************************************************************
EB_EXP - full data AttributesVecs (1)
*********************************************************************************************/
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());
	for (int treeNum=0; treeNum<etVec.size(); ++treeNum) {
		for (int spNum = 0; spNum<spVec.size(); ++spNum) {
            computePijGam cpg;
		    cpg.fillPij(etVec[treeNum],*(spVec[spNum]));
			for (int pos=0; pos < sc.seqLen(); ++pos) {
				if (((spAttributesVec[pos]-1)!=spNum ) || ((treeAttributesVec[pos]-1)!=treeNum )) continue;
				const tree*  treeForThisPosition=NULL;
				assert ((etVec.size() >0 ) && (treeAttributesVec[pos]>0));
				treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
				const stochasticProcess* spForThisPosition=NULL;
				assert ((spVec.size() >0 ) && (spAttributesVec[pos]>0));
				spForThisPosition = spVec[ spAttributesVec[pos] -1];
				computeEB_EXP_siteSpecificRate(pos,sc,*spForThisPosition,cpg, *treeForThisPosition,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf);
				LOG(6,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
			}
		}
	}
	LOG(6,<<" number of sites: "<<sc.seqLen()<<endl);
}

/********************************************************************************************
EB_EXP -  AttributesVecs  - one tree many sps
*********************************************************************************************/
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
							       const sequenceContainer& sc,
								   const tree & et,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	Vint etAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,etAttributesVec,sc,etVec,spVec,alphaConf);
}

/********************************************************************************************
EB_EXP -  AttributesVecs  - one sp many trees
*********************************************************************************************/
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const stochasticProcess & sp,
								   const MDOUBLE alphaConf){
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,treeAttributesVec,sc,etVec,spVec,alphaConf);
}

