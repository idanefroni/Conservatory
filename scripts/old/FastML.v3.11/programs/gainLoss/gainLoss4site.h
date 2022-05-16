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


#ifndef ___GAIN_LOSS_4site
#define ___GAIN_LOSS_4site

#include "definitions.h"
#include "replacementModel.h"
#include "gainLoss.h"
#include "siteSpecificRate.h"
#include "siteSpecificGL.h"

/********************************************************************************************
gainLoss4site
*********************************************************************************************/
class gainLoss4site{
public:
	explicit gainLoss4site(sequenceContainer& sc, tree& tr, vector<vector<stochasticProcess*> > spVVec,distribution* gainDist,distribution* lossDist,
		string& outDir, unObservableData* unObservableData_p, MDOUBLE alphaConf= 0.05); 
	gainLoss4site(const gainLoss4site& other) {*this = other;}	
	gainLoss4site& operator=(const gainLoss4site &other);
	virtual ~gainLoss4site() {;}
	
	void computeGain4Site();
	void computeLoss4Site();
	void printGain4Site();
	void printLoss4Site();

	Vdouble  get_gainV(){return _gainV;};
	Vdouble  get_lossV(){return _lossV;};

	Vdouble  get_stdGainV(){return _stdGainV;};
	Vdouble  get_stdLossV(){return _stdLossV;};
	
	VVdouble  get_posteriorsGainV(){return _posteriorsGainV;};
	VVdouble  get_posteriorsLossV(){return _posteriorsLossV;};
	VVVdouble getLpostPerSpPerCat() {return _postProbPerSpPerCatPerPos;}
	void initializeLpostPerSpPerCat();


protected:
//func
	void printGainLossBayes(ostream& out, const Vdouble& rate2printV, const Vdouble& lowerBoundV, const Vdouble& upperBoundV,const VVdouble& posteriorV, const distribution* dist,const stochasticProcess* sp);

protected:
	vector<vector<stochasticProcess*> > _spVVec; //save stochasticProcess for each category
	distribution* _gainDist;
	distribution* _lossDist;

	VVVdouble _postProbPerSpPerCatPerPos; // the posterior probability for each stochastic process for each Cat for each site


	tree _tr;
	sequenceContainer _sc;
	sequence* _refSeq; // the reference sequence
	string _outDir;

	Vdouble  _gainV,_stdGainV,_lowerBoundGainV,_upperBoundGainV;							  
	VVdouble _posteriorsGainV;	
	
	Vdouble  _lossV,_stdLossV,_lowerBoundLossV,_upperBoundLossV;							  
	VVdouble _posteriorsLossV;
	MDOUBLE _alphaConf;
	unObservableData* _unObservableData_p;		// 

};


#endif
