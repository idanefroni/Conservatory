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
#include "rate4siteGL.h"
#include "gainLossUtils.h"
#include "gainLossAlphabet.h"
#include <cstring>
/********************************************************************************************
gainLoss4site
*********************************************************************************************/
rate4siteGL::rate4siteGL(sequenceContainer& sc, tree& tr, stochasticProcess* sp,  string& outDir, unObservableData* unObservableData_p):
_tr(tr),_sp(sp),_sc(sc),_outDir(outDir),_unObservableData_p(unObservableData_p)

//init:
{
	fillReferenceSequence();
	_alphaConf = 0.05;
}


rate4siteGL& rate4siteGL::operator=(const rate4siteGL &other){
	if (this != &other) {              // Check for self-assignment
	}
	return *this;
}


/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::run()
{
	LOGnOUT(4,<<"Running rate4site..."<<endl);
	computeRate4site();
	computeAveAndStd(); // put them in ave, and std
	normalizeRates(); // change also the ave, the std the quantiles, etc.
}

/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::printRates()
{
	string r4sNonNorm = _outDir + "//" + "rate4siteOrig.txt";
	ofstream nonNormalizedOutStream(r4sNonNorm.c_str());
	nonNormalizedOutStream.precision(PRECISION);//PRECISION
	printRates(nonNormalizedOutStream,_rates);
	nonNormalizedOutStream.close();
}
/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::printRatesNormalized()
{
	string r4sNorm = _outDir + "//" + "rate4site.txt";
	ofstream normalizedOutStream(r4sNorm.c_str());
	normalizedOutStream.precision(PRECISION);
	normalizedOutStream<<"# Rate values were normalized to Z score (mean rate=0, +/-1 rat= +/-standard error)"<<endl;
	printRates(normalizedOutStream,_normalizedRates);
	normalizedOutStream.close();	
}


/********************************************************************************************
computeRate4site
*********************************************************************************************/
Vdouble  rate4siteGL::computeRate4site()
{
	time_t t1;
	time(&t1);
	time_t t2;

	if (gainLossOptions::_rateEstimationMethod == gainLossOptions::ebExp) {
		LOGnOUT (4,<<"perform computeEB_EXP_siteSpecificRate... while computing posteriorProb PerCategory PerPosition"<<endl);
		_postProbPerCatPerPos.resize(_sp->categories());
		for (int rateIndex=0 ; rateIndex<_sp->categories(); ++rateIndex){
			_postProbPerCatPerPos[rateIndex].resize(_sc.seqLen());
		}
		computeEB_EXP_siteSpecificRate(_rates,_BayesianSTD,_BayesianLowerBound,_BayesianUpperBound,_sc,*_sp,_tr,_alphaConf,&_postProbPerCatPerPos,_unObservableData_p);
	}
	else if (gainLossOptions::_rateEstimationMethod == gainLossOptions::mlRate) {
		LOGnOUT (4,<<"perform computeML_siteSpecificRate with maxRate= "<<gainLossOptions::_maxRateForML<<endl);
		computeML_siteSpecificRate(_rates,_Lrate,_sc, *_sp,_tr, gainLossOptions::_maxRateForML);		
	}
	else 
		errorMsg::reportError("non such method for rate inference, in function void rate4site::computeRate4site()");

	time(&t2);
	LOGnOUT(4,<<"computeRate4site RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	return _rates;
}

/********************************************************************************************
printRates
*********************************************************************************************/
void rate4siteGL::printRates(ostream & out, const Vdouble & rate2print) {

	if (gainLossOptions::_rateDistributionType == gainLossOptions::GAMMA_MIXTURE){
		mixtureDistribution* pMixture = static_cast<mixtureDistribution*>(_sp->distr());
		pMixture->printParams(out);
	}
	switch (gainLossOptions::_rateEstimationMethod){
		case (gainLossOptions::ebExp):  
			printRatesBayes(out,rate2print);
			break;
		case (gainLossOptions::mlRate):
			printRatesML(out,rate2print);
			break;
	}
	printAveAndStd(out);
}

/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::printRatesML(ostream& out, const Vdouble & rate2print) {
	out<<"#Rates were calculated using Maximum Likelihood"<<endl;
	out<<"#SEQ: The presence(1) or Absence(0) in the reference sequence."<<"Displayed on sequence "<<_refSeq->name()<<endl;
	out<<"#SCORE: The conservation scores. lower value = higher conservation."<<endl;
	out<<"#MSA DATA: The number of aligned sequences having a character from the overall number of sequences at each position."<<endl;
	out<<endl;
	out<<"========================================================================================================================================================="<<endl;
	out<<"#POS"<<"\t"<<"SEQ"<<"\t"<<"SCORE"<<"\t"<<"MSA DATA"<<endl; // note position start from 1.
	out<<"========================================================================================================================================================="<<endl;

#ifdef unix
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		out<<pos+1<<"\t"<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<"\t"<<setprecision(7)<<rate2print[pos]<<"\t";
		out<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1.
	}
#else
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		out<<left<<pos+1<<left<<"\t"<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<"\t";
		out<<left<<setprecision(7)<<fixed<<rate2print[pos]<<"\t";
		out<<right<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1. 
	}
#endif
}
/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::printRatesBayes(ostream& out, const Vdouble & rate2print) {
	int precisionHigh = 5;
	int precisionLow = 3;

	out<<"# Rates were calculated using the expectation of the posterior rate distribution"<<endl;
	out<<"# Prior distribution is Gamma with "<<gainLossOptions::_numberOfRateCategories<<" discrete categories"<<endl;
	//out<<"# SEQ: The presence(1) or Absence(0) in the reference sequence."<<"Displayed on sequence "<<_refSeq->name()<<endl;
	//out<<"# SCORE: The conservation scores. lower value = higher conservation."<<endl;
	//out<<"# QQ-INTERVAL: the confidence interval for the rate estimates. The default interval is 25-75 percentiles"<<endl;
	//out<<"# STD: the standard deviation of the posterior rate distribution."<<endl;
	//out<<"# MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position."<<endl;
	MDOUBLE AlphaRate = getRateAlpha(_sp->distr());
	out<<"# The alpha parameter "<<AlphaRate<<endl;
	int k=0;
	while (k < _sp->categories()){
		out<<"# sp.rates(j)  j= " <<k<<"\t"<<_sp->rates(k)<<"\t"<<_sp->ratesProb(k)<<endl;
		k++;
	}		

	out<<endl;
	out<<"========================================================================================================================================================="<<endl;
	//out<<"# POS"<<"\t"<<"SEQ"<<"\t"<<"SCORE"<<"\t"<<"QQ-INTERVAL"<<"\t"<<"STD"<<"\t"<<"MSA DATA"<<endl; // note position start from 1.
	out<<"# POS"<<"\t"<<"RATE"<<endl; // note position start from 1.
	//out<<"========================================================================================================================================================="<<endl;


	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		out<<pos+1<<"\t"<<rate2print[pos]<<endl;
	}

//#ifdef unix	
//	for (int pos=0; pos < _sc.seqLen(); ++pos) {
//		out<<pos+1<<"\t"
//			//<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<"\t"
//		out<<setprecision(precisionHigh)<<rate2print[pos]<<"\t";
//		//out<<"["<<setprecision(precisionLow)<<_BayesianLowerBound[pos]<<","<<setprecision(precisionLow)<<_BayesianUpperBound[pos]<<"]"<<"\t";
//		//out<<setprecision(precisionLow)<<_BayesianSTD[pos]<<"\t";
//		//out<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()
//		out<<endl; // note position start from 1.
//	}
//#else
//	for (int pos=0; pos < _sc.seqLen(); ++pos) {
//		out<<left<<pos+1;
//		//out<<left<<"\t"<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<"\t";
//		out<<left<<setprecision(precisionHigh)<<fixed<<rate2print[pos]<<"\t";
//		//out<<right<<"["<<setprecision(precisionLow)<<left<<_BayesianLowerBound[pos]<<","<<setprecision(precisionLow)<<right<<_BayesianUpperBound[pos]<<"]"<<"\t";
//		//out<<right<<setprecision(precisionLow)<<_BayesianSTD[pos];
//		//out<<right<<"\t"<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()
//		out<<endl; // note position start from 1.
//	}
//#endif
}
/********************************************************************************************
*********************************************************************************************/
void rate4siteGL::printAveAndStd(ostream& out) {
	out<<"# Average = "<<_ave<<endl;
	out<<"# Standard Deviation = "<<_std<<endl;
}
/********************************************************************************************
computeAveAndStd
*********************************************************************************************/
void rate4siteGL::computeAveAndStd(){
	MDOUBLE sum = 0;
	MDOUBLE sumSqr=0.0;
	for (int i=0; i < _sc.seqLen(); ++i) {
		sum+=_rates[i];
		sumSqr+=(_rates[i]*_rates[i]);
	}
	_ave = sum/_sc.seqLen();
	_std= sumSqr-(sum*sum/_sc.seqLen());
	_std /= (_sc.seqLen()-1.0);
	_std = sqrt(_std);
	if (((_ave<1e-9)) && (_ave>(-(1e-9)))) _ave=0;
	if ((_std>(1-(1e-9))) && (_std< (1.0+(1e-9)))) _std=1.0;
}
/********************************************************************************************
normalizeRates
*********************************************************************************************/
void rate4siteGL::normalizeRates() {
	int i=0;
	if (_std==0){
		LOGnOUT(4,<<"ERROR:\n std = 0 in function normalizeRates\n");
	}
	_normalizedRates.resize(_sc.seqLen(),0.0);
	for (i=0;i<_normalizedRates.size();++i) {
		_normalizedRates[i]=(_rates[i]-_ave)/_std;
	}

	if (gainLossOptions::_rateEstimationMethod == gainLossOptions::ebExp) {
		for (int k=0; k < _sc.seqLen(); ++k) {
			_BayesianUpperBound[k] = (_BayesianUpperBound[k] - _ave)/_std;
			_BayesianLowerBound[k] = (_BayesianLowerBound[k] - _ave)/_std;
			_BayesianSTD[k] = (_BayesianSTD[k])/_std;
		}
	}
	//_ave = 0.0;
	//_std = 1.0;
}

/********************************************************************************************
normalizeRates
*********************************************************************************************/
void rate4siteGL::fillReferenceSequence(){
	if (strcmp(gainLossOptions::_referenceSeq.c_str(),"non")==0) {
		_refSeq = &(_sc[0]);
	}
	else {
		int id1 = _sc.getId(gainLossOptions::_referenceSeq,true);
		_refSeq = (&_sc[id1]);
	}
}
