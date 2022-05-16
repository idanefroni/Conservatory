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


#ifndef ___GAIN_LOSS_MODEL
#define ___GAIN_LOSS_MODEL

#include "definitions.h"
#include "replacementModel.h"
#include "fromQtoPt.h"
#include "errorMsg.h"
#include "matrixUtils.h"
#include "gainLossUtils.h"
#include "gainLossOptions.h"

/********************************************************************************************
Q is a matrix of the following form:
(where 0 and 1 stand for absence or presence)
for a reversible case,
0			1
0	-m1			m1
1	m1*pi0/pi1	-m1*pi0/pi1

and without assuming reversibility,
0			1
0	-m1			m1(gain)
1	m2(loss)	-m2


1. The gainLossModel class is derived from the general replacementModel class - it models the stochastic process with one param gain=loss
2. Additionally we use the gainLossModelNonReversible class which is derived from gainLossModel class - we get the second param - gain!=loss
*********************************************************************************************/

/********************************************************************************************
gainLossModel
*********************************************************************************************/
class gainLossModel : public replacementModel {
public:
	explicit gainLossModel(const MDOUBLE m1, const Vdouble freq, bool isRootFreqEQstationary, bool isReversible, bool isHGT_normal_Pij,  bool _isHGT_with_Q); 
	virtual replacementModel* clone() const { 
		return new gainLossModel(*this); 
	}
	gainLossModel(const gainLossModel& other): _q2pt(NULL) {*this = other;}
	virtual gainLossModel& operator=(const gainLossModel &other);

	virtual ~gainLossModel() {if (_q2pt) delete _q2pt; }
	const int alphabetSize() const {return 2;} // assumes only absence or presence
	const MDOUBLE err_allow_for_pijt_function() const {return 1e-4;} // same as q2p definitions
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const int i) const;
	void setMu1(const MDOUBLE val, bool isReversible);
	MDOUBLE setTheta(const MDOUBLE val);
	MDOUBLE getTheta() const {return _freq[1];}

	bool isReversible() const {return _isReversible;}
	bool isRootFreqEQstationary() const {return _isRootFreqEQstationary;}	
	bool isHGT_normal_Pij() const {return _isHGT_normal_Pij;}
	bool isHGT_with_Q() const {return _isHGT_with_Q;}

	const VVdouble getQ() const {return _Q;}
	const MDOUBLE getMu1() const {return _Q[0][1];}
	const MDOUBLE getMu2() const {return _Q[1][0];}

	void norm(const MDOUBLE scale);
	MDOUBLE sumPijQij();

	//const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const{
	//	return _q2pt->Pij_t(i,j,d);
	//}
	//const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
	//	return _q2pt->dPij_dt(i,j,d);
	//}
	//const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
	//	return _q2pt->d2Pij_dt2(i,j,d);
	//}


protected:
	virtual MDOUBLE updateQ(bool isReversible);
	virtual MDOUBLE normalizeQ();

	bool pijt_is_prob_value(MDOUBLE val) const {
		if ((abs(val)+err_allow_for_pijt_function()<0) || (val>1+err_allow_for_pijt_function()))
			return false;
		else
			return true;
	}
protected:
	Vdouble _freq;
	Vdouble _freqQ;
	MDOUBLE _rQ;
	MDOUBLE _gain;	// _Q[0][1]
	VVdouble _Q;	
	q2pt *_q2pt; // dont use q2p
	bool _isReversible;
	bool _isRootFreqEQstationary;
	bool _isHGT_normal_Pij;
	bool _isHGT_with_Q;

};



/********************************************************************************************
gainLossModelNonReversible
All the methods of this class are implemented in the header
*********************************************************************************************/
class gainLossModelNonReversible : public gainLossModel {
public:
	//////////////////////////////////////////////////////////////////////////
	explicit gainLossModelNonReversible(const MDOUBLE m1, const MDOUBLE m2, const Vdouble freq,bool isRootFreqEQstationary, bool isHGT_normal_Pij,  bool _isHGT_with_Q)
		:_loss(m2),gainLossModel(m1,freq,isRootFreqEQstationary,false,isHGT_normal_Pij,_isHGT_with_Q) 	
	{
		updateQ();
		if(_isRootFreqEQstationary) {	
			setTheta(getMu1()/(getMu1()+getMu2()));
		}
	}
	//////////////////////////////////////////////////////////////////////////
	virtual replacementModel* clone() const { 
		return new gainLossModelNonReversible(*this); 
	}
	gainLossModelNonReversible(const gainLossModelNonReversible& other)  : gainLossModel(other)
	{
		_loss = other._loss;
	}	
	virtual ~gainLossModelNonReversible(){
		//cout<<"gainLossModelNonReversible Deleted\n";
	}	
	//gainLossModelNonReversible& operator=(const gainLossModelNonReversible &other) 
	//{
	//	_loss = other._loss;
	//	return *this;
	//}


	//////////////////////////////////////////////////////////////////////////
	void setMu2(const MDOUBLE val) {
		_loss = val; 
		updateQ();
		if(_isRootFreqEQstationary) {	
			setTheta(getMu1()/(getMu1()+getMu2()));
		}
		//if(gainLossOptions::_isNormalizeQ)  // part of update Q
		//	normalizeQ();

	}
	//const MDOUBLE getMu2() const {return _loss;}	// moved to gainLossModel
	//const VVdouble getQ() const {return _Q;} // moved to gainLossModel


protected:
	//virtual void updateQ();
	//////////////////////////////////////////////////////////////////////////
	void updateQ(){
		//gainLossModel::updateQ(false);
		_Q[1][1] = -_loss;
		_Q[1][0] =  _loss;
		if(gainLossOptions::_isNormalizeQ && !gainLossOptions::_gainLossDist && (_Q[1][0]>0))//? 
			normalizeQ(); 
	}

	//bool pijt_is_prob_value(MDOUBLE val) const { // moved to gainLossModel
	//	if ((abs(val)+err_allow_for_pijt_function()<0) || (val>1+err_allow_for_pijt_function()))
	//		return false;
	//	else
	//		return true;
	//}
private:
	MDOUBLE _loss; // _Q[1][0]
};

#endif
