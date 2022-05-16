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
#include "gainLossModel.h"

/********************************************************************************************
gainLossModel
Note: All gainLossOptions parameter are sent
to the c'tor as a preperation for the model to be part of the Lib. 
*********************************************************************************************/
gainLossModel::gainLossModel(const MDOUBLE m1, const Vdouble freq, bool isRootFreqEQstationary, bool isReversible, bool isHGT_normal_Pij, bool isHGT_with_Q):
_gain(m1),_freq(freq),_isRootFreqEQstationary(isRootFreqEQstationary),_isReversible(isReversible),_isHGT_normal_Pij(isHGT_normal_Pij),_isHGT_with_Q(isHGT_with_Q),_q2pt(NULL){
	if (freq.size() != alphabetSize())
		errorMsg::reportError("Error in gainLossModel, size of frequency vector must be as in alphabet");
	for(int i=0; i<freq.size(); ++i)
		if(freq[i]<0 || freq[i]>1)
			errorMsg::reportError("Freq not within [0,1]\n");
	if(!_isHGT_with_Q){_gain = 0;}
	resizeMatrix(_Q,alphabetSize(),alphabetSize());
	updateQ(_isReversible);
	//setTheta(_freq[1]);	// no Need
	if(_isRootFreqEQstationary) {	
		setTheta(getMu1()/(getMu1()+getMu2()));
	} 
}
/********************************************************************************************
*********************************************************************************************/
gainLossModel& gainLossModel::operator=(const gainLossModel &other){
	if (this != &other) {              // Check for self-assignment
		if (_q2pt) delete _q2pt;
		if (other._q2pt != NULL) 
			_q2pt = (q2pt*)(other._q2pt->clone());
	}
	_isReversible = other.isReversible();
	_isRootFreqEQstationary = other.isRootFreqEQstationary();
	_isHGT_normal_Pij = other.isHGT_normal_Pij();
	_isHGT_with_Q = other.isHGT_with_Q();
	_gain = other._gain;
	_freq = other._freq;
	_Q = other._Q;
	return *this;
}
/********************************************************************************************
*********************************************************************************************/
void gainLossModel::setMu1(const MDOUBLE val, bool isReversible) { 
	if(_isHGT_with_Q) {_gain = val;}
	updateQ(isReversible);
	if(_isRootFreqEQstationary) {	
		setTheta(getMu1()/(getMu1()+getMu2()));
	}
	//if(gainLossOptions::_isNormalizeQ) // part of update Q
	//	normalizeQ();
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLossModel::setTheta(const MDOUBLE val) {
	if(val<0 || val>1)
		errorMsg::reportError("Freq not within [0,1]\n");
	_freq[1]=val; 
	_freq[0]= 1-val; 
	MDOUBLE normFactor =  updateQ(_isReversible);
	return normFactor;
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLossModel::updateQ(bool isReversible){
	MDOUBLE normFactor=1;
	_Q[0][1] = _gain;
	_Q[0][0] = -_Q[0][1];

	if (isReversible) {
		_Q[1][0] = _Q[0][1] * _freq[0] / _freq[1]; // m1*pi0/pi1
		_Q[1][1] = -_Q[1][0]; 
	}
	//else{
	//	_Q[1][0] = 1;	//To be overwritten by gainLossModelNonReversible
	//	_Q[1][1] = -1;	//To be overwritten by gainLossModelNonReversible
	//}
	if (gainLossOptions::_gainEQloss) {
		_Q[1][0] = _gain;
		_Q[1][1] = -_Q[1][0]; 
	}
	if (gainLossOptions::_gainLossRateAreFreq) {
		_Q[1][0] = 1 - _gain;
		_Q[1][1] = -_Q[1][0]; 
	}

	for (int i=0; i<_Q.size();i++) {
		MDOUBLE sum = _Q[i][0]+_Q[i][1];
		if ((abs(sum)>err_allow_for_pijt_function()))
			errorMsg::reportError("Error in gainLossModel::updateQ, sum of row is not 0");
	}
	//if (isReversible){
	//	if (!_q2pt) 
	//		_q2pt = new q2pt();
	//	_q2pt->fillFromRateMatrix(_freq,_Q); 
	//}	
	if(gainLossOptions::_isNormalizeQ && !gainLossOptions::_gainLossDist  && (_Q[1][0]>0)) // 
		normFactor= normalizeQ();
	return normFactor;
}
/********************************************************************************************
*********************************************************************************************/
const MDOUBLE gainLossModel::freq(const int i) const {
	if (i >= _freq.size()) 
		errorMsg::reportError("Error in gainLossModel::freq, i > size of frequency vector");
	return _freq[i];
}
/********************************************************************************************
// normalize Q so that sum of changes = 1
*********************************************************************************************/
MDOUBLE gainLossModel::normalizeQ(){ 
	MDOUBLE norm_factor=0.0;
	for (int i=0;i<_Q.size();i++)
		norm_factor+=(_freq[i]*_Q[i][i]);
	MDOUBLE fac = -1.0/norm_factor;
	_Q = multiplyMatrixByScalar(_Q,fac);
	return fac;
}
/********************************************************************************************
*********************************************************************************************/
void gainLossModel::norm(const MDOUBLE scale)
{
	for (int i=0; i < _Q.size(); ++i) {
		for (int j=0; j < _Q.size(); ++j) {
			_Q[i][j] *= scale; 		
		}
	}
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE gainLossModel::sumPijQij(){
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= (_Q[i][i])*_freq[i];
	}
	return sum;
}

/********************************************************************************************
Pij_t - Based on Analytic solution
*********************************************************************************************/
const MDOUBLE gainLossModel::Pij_t(const int i,const int j, const MDOUBLE d)  const {
	MDOUBLE gain = getMu1();
	MDOUBLE loss = getMu2();
	MDOUBLE eigenvalue =  -(gain + loss);
	bool withHGT = isHGT_normal_Pij();
	
	MDOUBLE noHGTfactor = 0.0001;

	VVdouble Pt;
	resizeMatrix(Pt,_Q.size(),_Q.size());
	int caseNum = i + j*2;
	switch (caseNum) {
		case 0 : Pt[0][0] =  loss/(-eigenvalue) + exp(eigenvalue*d)*(1 - loss/(-eigenvalue)); break;
		case 1 : Pt[1][0] =  loss/(-eigenvalue) - exp(eigenvalue*d)*(1 - gain/(-eigenvalue)); break;
		case 2 : if(withHGT)
				 { Pt[0][1] =  gain/(-eigenvalue) - exp(eigenvalue*d)*(1 - loss/(-eigenvalue));}
				 else 
				 { Pt[0][1] =  (gain/(-eigenvalue) - exp(eigenvalue*d)*(1 - loss/(-eigenvalue)))*noHGTfactor;} break;
		case 3 : Pt[1][1] =  gain/(-eigenvalue) + exp(eigenvalue*d)*(1 - gain/(-eigenvalue));  break;
	}
	MDOUBLE val = (Pt[i][j]);
	if (!pijt_is_prob_value(val)){
		string err = "Error in gainLossModelNonReversible::Pij_t, pijt <0 or >1. val=";
		err+=double2string(val);
		err+=" d=";
		err+=double2string(d);		
		LOG(4,<<err<<endl);		//errorMsg::reportError(err);
	}
	if(!(val>VERYSMALL))
		val = VERYSMALL;
	LOG(10,<<"for gain "<<gain<<" loss "<<loss<<"  P"<<i<<j<<"("<<d<<")  "<<val<<endl;)
	return val; 
}

/********************************************************************************************
dPij_t - Based on Analytic solution
*********************************************************************************************/
const MDOUBLE gainLossModel::dPij_dt(const int i,const int j, const MDOUBLE d)  const {
	MDOUBLE gain = getMu1();;
	MDOUBLE loss = getMu2();;
	MDOUBLE eigenvalue =  -(gain + loss);

	VVdouble Pt;
	resizeMatrix(Pt,_Q.size(),_Q.size());
	int caseNum = i + j*2;
	switch (caseNum) {
		case 0 : Pt[0][0] =    exp(eigenvalue*d)*(eigenvalue + loss); break;
		case 1 : Pt[1][0] =  -(exp(eigenvalue*d)*(eigenvalue + gain)); break;
		case 2 : Pt[0][1] =  -(exp(eigenvalue*d)*(eigenvalue + loss));  break;
		case 3 : Pt[1][1] =    exp(eigenvalue*d)*(eigenvalue + gain);  break;
	}
	MDOUBLE val = (Pt[i][j]);
	//if (!pijt_is_prob_value(val)){
	//	string err = "Error in gainLossModelNonReversible::dPij_t_dt, pijt <0 or >1. val=";
	//	err+=double2string(val);
	//	err+=" d=";
	//	err+=double2string(d);		
	//	LOG(6,<<err<<endl);	//errorMsg::reportError(err);
	//}
	return val; 
}
/********************************************************************************************
d2Pij_dt2 - Based on Analytic solution
*********************************************************************************************/
const MDOUBLE gainLossModel::d2Pij_dt2(const int i,const int j, const MDOUBLE d)  const {
	MDOUBLE gain = getMu1();;
	MDOUBLE loss = getMu2();;
	MDOUBLE eigenvalue =  -(gain + loss);

	VVdouble Pt;
	resizeMatrix(Pt,_Q.size(),_Q.size());
	int caseNum = i + j*2;
	switch (caseNum) {
		case 0 : Pt[0][0] =    exp(eigenvalue*d)*(eigenvalue + loss)*eigenvalue; break;
		case 1 : Pt[1][0] =  -(exp(eigenvalue*d)*(eigenvalue + gain))*eigenvalue; break;
		case 2 : Pt[0][1] =  -(exp(eigenvalue*d)*(eigenvalue + loss))*eigenvalue;  break;
		case 3 : Pt[1][1] =    exp(eigenvalue*d)*(eigenvalue + gain)*eigenvalue;  break;
	}
	MDOUBLE val = (Pt[i][j]);
	//if (!pijt_is_prob_value(val)){
	//	string err = "Error in gainLossModelNonReversible::d2Pij_t_dt2, pijt <0 or >1. val=";
	//	err+=double2string(val);
	//	LOG(6,<<err<<endl);			//errorMsg::reportError(err);
	//}
	return val; 
}



/********************************************************************************************
non reversible model
updateQ
*********************************************************************************************/
//void gainLossModelNonReversible::updateQ(){
//	//gainLossModel::updateQ(false);
//	_Q[1][1] = -_loss;
//	_Q[1][0] =  _loss;
//	//normalizeQ();
//}





/********************************************************************************************
Pij_t - converging series
 IMPORTANT NOTE: this function is VERY inefficient. It calculates all of Pt for every call of Pijt
 this is unimportant for a small dataset (one position) but pre-processing should be done for larger datasets:

 SOLUTION: save the computed Pijt matrix each time it is called. In every call of Pij_t, check if a saved value exists
*********************************************************************************************/
//const MDOUBLE gainLossModelNonReversible::Pij_t(const int i,const int j, const MDOUBLE d)  const {
//
//	VVdoubleRep QdblRep; 
//	resizeMatrix(QdblRep,_Q.size(),_Q.size());
//	for (int row=0;row<_Q.size();row++){
//			for (int col=0;col<_Q[row].size();col++)
//				QdblRep[row][col]=convert(_Q[row][col]);
//	}
//	VVdoubleRep Qt = multiplyMatrixByScalar(QdblRep,d);
//	VVdoubleRep unit;
//	unitMatrix(unit,_Q.size());
//	VVdoubleRep Pt = add(unit,Qt) ; // I + Qt
//	VVdoubleRep Qt_power = Qt;
//	doubleRep old_val = Pt[i][j];
//	doubleRep diff(1.0);
//	int n=2;
//	while ((diff>err_allow_for_pijt_function()) || (!pijt_is_prob_value(convert(Pt[i][j])))){//(abs(old_val-new_val) > err_allow_for_pijt_function()){
//		old_val = Pt[i][j];
//		Qt_power = multiplyMatrixes(Qt_power,multiplyMatrixByScalar(Qt,1.0/n));
//		Pt= add(Pt,Qt_power); // I + Qt + Qt^2/2! + ....  + Qt^n/n!
//
//		diff = Pt[i][j]-old_val; // difference is measured by diff between P[0][0] vals (a little primitive...)
//		if (diff<0) diff=-diff;
//		n++;
//		if (n>200) { 
//			string err = "Error in gainLossModelNonReversible::Pij_t, too many (>n=200) iterations for t = " + double2string(d);
//			cerr<<diff<<endl;
//			errorMsg::reportError(err);
//		}
//	}
//	MDOUBLE val = convert(Pt[i][j]);
//	if (!pijt_is_prob_value(val))
//		errorMsg::reportError("Error in gainLossModelNonReversible::Pij_t, pijt <0 or >1");
//	LOG(10,<<"for gain "<<getMu1()<<" loss "<<getMu2()<<"  P"<<i<<j<<"("<<d<<")  "<<val<<endl;)
//
//	return val; 
//}
//
