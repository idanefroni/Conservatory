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
#include "gainLoss4site.h"

/********************************************************************************************
gainLoss4site
*********************************************************************************************/
gainLoss4site::gainLoss4site(sequenceContainer& sc, tree& tr, 
							 vector<vector<stochasticProcess*> > spVVec,distribution* gainDist,distribution* lossDist,   
							 string& outDir, unObservableData* unObservableData_p, MDOUBLE alphaConf):
_tr(tr),_spVVec(spVVec),_gainDist(gainDist),_lossDist(lossDist),_sc(sc),_outDir(outDir),_unObservableData_p(unObservableData_p),_alphaConf(alphaConf)
{
	//init:
	_refSeq = &(_sc[0]);
}

gainLoss4site& gainLoss4site::operator=(const gainLoss4site &other){
	if (this != &other) {              // Check for self-assignment
	}
	return *this;
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::computeGain4Site()
{
	LOGnOUT (4,<<"perform computeGain4Site... while computing posteriorProb PerCategory PerPosition"<<endl);
	
	initializeLpostPerSpPerCat();
	computeEB_EXP_siteSpecificGL(_gainV, _stdGainV, _lowerBoundGainV, _upperBoundGainV, _posteriorsGainV, _sc, _spVVec, _tr, _gainDist,_lossDist,_gainDist,
		_alphaConf,_postProbPerSpPerCatPerPos,_unObservableData_p);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::computeLoss4Site()
{
	LOGnOUT (4,<<"perform computeLoss4Site... while computing posteriorProb PerCategory PerPosition"<<endl);
	initializeLpostPerSpPerCat();
	computeEB_EXP_siteSpecificGL(_lossV, _stdLossV, _lowerBoundLossV, _upperBoundLossV, _posteriorsLossV, _sc, _spVVec, _tr, _gainDist,_lossDist,_lossDist,
		_alphaConf,_postProbPerSpPerCatPerPos,_unObservableData_p);
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::printGain4Site()
{
	//ofstream outGain(gainLossOptions::_outFileGain4Site.c_str());
	string g4s = _outDir + "//" + "gain4site.txt";
	ofstream outGain(g4s.c_str());
	outGain.precision(PRECISION);
	printGainLossBayes(outGain,_gainV,_lowerBoundGainV,_upperBoundGainV, _posteriorsGainV, _gainDist,_spVVec[0][0]);
	outGain.close();

}
/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::printLoss4Site()
{
	//ofstream outLoss(gainLossOptions::_outFileLoss4Site.c_str());
	string l4s = _outDir + "//" + "loss4site.txt";
	ofstream outLoss(l4s.c_str());
	outLoss.precision(PRECISION);
	printGainLossBayes(outLoss,_lossV,_lowerBoundLossV,_upperBoundLossV, _posteriorsLossV, _lossDist,_spVVec[0][0]);
	outLoss.close();
}

/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::printGainLossBayes(ostream& out, const Vdouble& rate2printV, const Vdouble& lowerBoundV, const Vdouble& upperBoundV,const VVdouble& posteriorV, const distribution* dist, const stochasticProcess* sp) 
{	
	out.precision(7);
	out<<"# Empirical Bayesian Rates"<<endl;
	//out<<"#Displayed on sequence "<<_refSeq->name()<<endl;
	if(sp->categories() > 1){
		out<<"# each sp with overall rate distribution cat: ";
		for (int cat = 0; cat < sp->categories(); ++cat)
			out<<sp->rates(cat)<<" ";
		out<<endl;
	}
	out<<"# Rate in each gamma category: ";
	for (int cat = 0; cat <dist->categories(); ++cat)
		out<<"category "<<cat+1<<", Rate= "<<dist->rates(cat)<<" ";
	out<<endl;
	out<<"# Posterior probability for each category, and each position is given.\n";

	//out<<"========================================================================================================================================================="<<endl;
	out<<"POS\t"<<"Rate\t";//<<"[Confidence Interval]\t";
	for (int cat = 0; cat <dist->categories(); ++cat)
		out<<"Categ "<<cat+1<<"\t";
	out<<endl;
	int numOfCategories = dist->categories();
	for (int i=0;i<_sc.seqLen();i++){	 
		//string aaStr = _refSeq->getAlphabet()->fromInt((*_refSeq)[i]);
		out<<i+1 /*<<"\t"<<aaStr*/<<"\t"<< rate2printV[i]<<"\t";
		//<<"["<<lowerBoundV[i]<<","<<upperBoundV[i]<<"]\t";
		//if (lowerBoundV[i]>1) out <<"*"; //significance indicator: if entire confidence interval >1
		for (int cat = 0; cat < numOfCategories; ++cat)
			out<<posteriorV[i][cat]<<"\t";
		out<<endl;
	}
}	

/********************************************************************************************
*********************************************************************************************/
void gainLoss4site::initializeLpostPerSpPerCat()
{
	int numOfSPs = _gainDist->categories()*_lossDist->categories();
	int rateCategories = _spVVec[0][0]->categories();
	if(_postProbPerSpPerCatPerPos.size()==0){
		resizeVVV(numOfSPs,rateCategories,_sc.seqLen(),_postProbPerSpPerCatPerPos);
	}
}


