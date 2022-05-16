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


#ifndef ___RATE_4_site___GL
#define ___RATE_4_site___GL

#include "definitions.h"
#include "replacementModel.h"
#include "gainLoss.h"
#include "unObservableData.h"

/********************************************************************************************
rate4siteGL
*********************************************************************************************/
class rate4siteGL{
public:
	explicit rate4siteGL(sequenceContainer& sc, tree& tr, stochasticProcess* sp,  string& outDir, unObservableData* unObservableData_p);

	rate4siteGL(const rate4siteGL& other) {*this = other;}	
	rate4siteGL& operator=(const rate4siteGL &other);
	virtual ~rate4siteGL() {;}
	void run();
	VVdouble getLpostPerCat() {return _postProbPerCatPerPos;}
	Vdouble getRates() {return _rates;}
	Vdouble getNormalizedRates() {return _normalizedRates;}

	void printRatesNormalized();
	void printRates();

protected:
//func
	Vdouble computeRate4site();
	void computeAveAndStd();
	void normalizeRates();
	void printRates(ostream & out, const Vdouble & rate2print);
	void printRatesML(ostream& out, const Vdouble & rate2print);
	void printRatesBayes(ostream& out, const Vdouble & rate2print);
	void printAveAndStd(ostream& out= cout);
	void fillReferenceSequence();

protected:
//members
	stochasticProcess *_sp; 
	tree _tr;
	sequenceContainer _sc;
	sequence* _refSeq; // the reference sequence
	string _outDir;

	Vdouble _rates;// the rates themselves
	Vdouble _Lrate;// the log likelihood of each position
	VVdouble _postProbPerCatPerPos; // the posterior probability for each category for each site
	Vdouble _normalizedRates; // the rates when their ave = 0 and std = 1.
	MDOUBLE _ave; // the average over all rates.
	MDOUBLE _std; // the std over all rates.
	Vdouble _BayesianSTD;// the std of the Bayesian rates
	Vdouble _BayesianLowerBound;// lower bound of rate in Bayesian inference
	Vdouble _BayesianUpperBound;// upper bound of rate in Bayesian inference
	MDOUBLE _alphaConf; // the alpha confidence interval of Bayesian rates (set to 0.5). interval - rates that are in the 95% area under the curve.
	unObservableData* _unObservableData_p;		// 
};


#endif
