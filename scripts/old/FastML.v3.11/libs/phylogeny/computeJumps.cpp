#include "computeJumps.h"
#include "talRandom.h"
#include "someUtil.h"
#include "matrixUtils.h"
#include <algorithm>


computeJumps::computeJumps(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r, const int maxNumOfChangesPerBranchSum)
: _Lambda1(Lambda1), _Lambda2(Lambda2),_maxNumOfChangesPerBranchSum(maxNumOfChangesPerBranchSum)
{
	if(_Lambda1==_Lambda2)
		_Lambda1+=EPSILON; // Patch: fix a BUG, if gain==loss the probability of transition from 0 to 1 given states start==End==1, is NA, thus add epsilon

	_gFuncStart0		= gFunc(_Lambda1, _Lambda2, r);
	_gFuncStart0MinusR	= gFunc(_Lambda1, _Lambda2, -r);
	_gFuncStart1		= gFunc(_Lambda2, _Lambda1, r);
	_gFuncStart1MinusR	= gFunc(_Lambda2, _Lambda1, -r);
}
computeJumps::~computeJumps()
{
}


/********************************************************************************************
getExpectation
*********************************************************************************************/
MDOUBLE computeJumps::getExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId)
{
	if(BranchLength>=0){
		if(fromId==0 && toId==1){ // Gain
			if(terminalStart==0 && terminalEnd==1)
				return gainExpGiven01(BranchLength);
			if(terminalStart==0 && terminalEnd==0)
				return gainExpGiven00(BranchLength);
			if(terminalStart==1 && terminalEnd==1)
				return gainExpGiven11(BranchLength);
			else //(terminalStart==1 && terminalEnd==0)
				return gainExpGiven10(BranchLength);	
		}
		if(fromId==1 && toId==0){ // Loss
			if(terminalStart==0 && terminalEnd==1)
				return lossExpGiven01(BranchLength);	
			if(terminalStart==0 && terminalEnd==0)
				return lossExpGiven00(BranchLength);
			if(terminalStart==1 && terminalEnd==1)
				return lossExpGiven11(BranchLength);
			else //(terminalStart==1 && terminalEnd==0)
				return lossExpGiven10(BranchLength);
		}
		else
			return 0;
	}
	else
		return 0;

}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE computeJumps::getTotalExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd)
{
	if(BranchLength>=0){
			if(terminalStart==0 && terminalEnd==1)
				return m01(BranchLength);
			if(terminalStart==0 && terminalEnd==0)
				return m00(BranchLength);
			if(terminalStart==1 && terminalEnd==1)
				return m11(BranchLength);
			else //(terminalStart==1 && terminalEnd==0)
				return m10(BranchLength);
	}
	else
		return 0;

}


/********************************************************************************************
gainExpGivenXY lossExpGivenXY
// Note: divide by Pij, since the computation is gainExp and End=0 given start=0
*********************************************************************************************/
MDOUBLE computeJumps::gainExpGiven01(MDOUBLE BranchLength){
	return 0.5*(m01(BranchLength) +Pij_t(0,1,BranchLength))/Pij_t(0,1,BranchLength);
}
MDOUBLE computeJumps::gainExpGiven00(MDOUBLE BranchLength){
	return 0.5*(m00(BranchLength)/Pij_t(0,0,BranchLength));
}
MDOUBLE computeJumps::gainExpGiven11(MDOUBLE BranchLength){
	return 0.5*(m11(BranchLength)/Pij_t(1,1,BranchLength) ); //???
}
MDOUBLE computeJumps::gainExpGiven10(MDOUBLE BranchLength){
	return m10(BranchLength)/Pij_t(1,0,BranchLength) - lossExpGiven10(BranchLength); //???
}
//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::lossExpGiven01(MDOUBLE BranchLength){
	return m01(BranchLength)/Pij_t(0,1,BranchLength) - gainExpGiven01(BranchLength); //???
}
MDOUBLE computeJumps::lossExpGiven00(MDOUBLE BranchLength){
	return m00(BranchLength)/Pij_t(0,0,BranchLength)  - gainExpGiven00(BranchLength); //???
}
MDOUBLE computeJumps::lossExpGiven11(MDOUBLE BranchLength){
	return m11(BranchLength)/Pij_t(1,1,BranchLength) - gainExpGiven11(BranchLength); //???
}
MDOUBLE computeJumps::lossExpGiven10(MDOUBLE BranchLength){
	return 0.5*(m10(BranchLength) + Pij_t(1,0,BranchLength) )/Pij_t(1,0,BranchLength); //???
	//return m10(BranchLength)/Pij_t(1,0,BranchLength) - gainExpGiven10(BranchLength); //???
}



/********************************************************************************************
getProbability
*********************************************************************************************/
MDOUBLE computeJumps::getProb(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId)
{
	if(BranchLength>=0){
		if(fromId==0 && toId==1){ // Gain
			if(terminalStart==0 && terminalEnd==1)
				return gainProbGiven01(BranchLength);
			if(terminalStart==0 && terminalEnd==0)
				return gainProbGiven00(BranchLength);
			if(terminalStart==1 && terminalEnd==1)
				return gainProbGiven11(BranchLength);
			else //(terminalStart==1 && terminalEnd==0)
				return gainProbGiven10(BranchLength); // if g=l, return -NaN
		}
		if(fromId==1 && toId==0){ // Loss
			if(terminalStart==0 && terminalEnd==1)
				return lossProbGiven01(BranchLength); // if g=l, return -NaN
			if(terminalStart==0 && terminalEnd==0)
				return lossProbGiven00(BranchLength);
			if(terminalStart==1 && terminalEnd==1)
				return lossProbGiven11(BranchLength);
			else //(terminalStart==1 && terminalEnd==0)
				return lossProbGiven10(BranchLength);
		}
		else
			return 0;
	}
	else
		return 0;

}
//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gainProbGiven01(MDOUBLE BranchLength){
	MDOUBLE probSum = 1.0;
	return probSum;
}
MDOUBLE computeJumps::gainProbGiven00(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(2,4,6,...) changes
	//for(int k = 1; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart0.qFunc_2k(BranchLength,k);
	//}
	//B 1 - Sum(uneven changes) - zeroEvenChanges
	probSum = 1 - 0.5*(_gFuncStart0.gFunc_(BranchLength) - _gFuncStart0MinusR.gFunc_(BranchLength)) - _gFuncStart0.qFunc_2k(BranchLength,0);
	return probSum/Pij_t(0,0,BranchLength);
}
MDOUBLE computeJumps::gainProbGiven11(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(2,4,6,...) changes
	//for(int k = 1; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart1.qFunc_2k(BranchLength,k);	//? _gFuncStart1 or _gFuncStart0
	//}
	//B 1 - Sum(uneven changes) - zeroEvenChanges
	probSum = 1 - 0.5*(_gFuncStart1.gFunc_(BranchLength) - _gFuncStart1MinusR.gFunc_(BranchLength)) - _gFuncStart1.qFunc_2k(BranchLength,0);
	return probSum/Pij_t(1,1,BranchLength);
}
MDOUBLE computeJumps::gainProbGiven10(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(3,5,7,...) changes
	//for(int k = 2; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart1.qFunc_2k_1(BranchLength,k);
	//}
	//B 1 - Sum(even changes) - oneUnEvenChanges
	probSum = 1 - 0.5*(_gFuncStart1.gFunc_(BranchLength) + _gFuncStart1MinusR.gFunc_(BranchLength)) - _gFuncStart1.qFunc_2k_1(BranchLength,1);
	return probSum/Pij_t(1,0,BranchLength);
}

//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::lossProbGiven01(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(3,5,7,...) changes
	//for(int k = 2; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart0.qFunc_2k_1(BranchLength,k);
	//}
	//B 1 - Sum(even changes) - oneUnEvenChanges
	probSum = 1 - 0.5*(_gFuncStart0.gFunc_(BranchLength) + _gFuncStart0MinusR.gFunc_(BranchLength)) - _gFuncStart0.qFunc_2k_1(BranchLength,1);
	return probSum/Pij_t(0,1,BranchLength);
}
MDOUBLE computeJumps::lossProbGiven00(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(2,4,6,...) changes
	//for(int k = 1; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart0.qFunc_2k(BranchLength,k);
	//}
	//B 1 - Sum(uneven changes) - zeroEvenChanges
	probSum = 1 - 0.5*(_gFuncStart0.gFunc_(BranchLength) - _gFuncStart0MinusR.gFunc_(BranchLength)) - _gFuncStart0.qFunc_2k(BranchLength,0);
	return probSum/Pij_t(0,0,BranchLength);
}

MDOUBLE computeJumps::lossProbGiven11(MDOUBLE BranchLength){
	MDOUBLE probSum = 0.0;
	//A Sum(2,4,6,...) changes
	//for(int k = 1; k<=_maxNumOfChangesPerBranchSum; ++k){
	//	probSum += _gFuncStart1.qFunc_2k(BranchLength,k);	//? _gFuncStart1 or _gFuncStart0
	//}
	//B 1 - Sum(uneven changes) - zeroEvenChanges
	probSum = 1 - 0.5*(_gFuncStart1.gFunc_(BranchLength) - _gFuncStart1MinusR.gFunc_(BranchLength)) - _gFuncStart1.qFunc_2k(BranchLength,0);
	return probSum/Pij_t(1,1,BranchLength);
}
MDOUBLE computeJumps::lossProbGiven10(MDOUBLE BranchLength){
	MDOUBLE probSum = 1.0;
	return probSum;
}


/********************************************************************************************
// mij(t) = E(N, end=j | start=i)
*********************************************************************************************/
MDOUBLE computeJumps::m01(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart0.gFunc_dr(BranchLength) - _gFuncStart0MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m00(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart0.gFunc_dr(BranchLength) + _gFuncStart0MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m11(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart1.gFunc_dr(BranchLength) + _gFuncStart1MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m10(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart1.gFunc_dr(BranchLength) - _gFuncStart1MinusR.gFunc_dr(BranchLength));
}

/********************************************************************************************
gFunc_dr
*********************************************************************************************/
MDOUBLE computeJumps::gFunc_dr(MDOUBLE BranchLength, int startState){
	// test:
	if(startState == 0){
		return _gFuncStart0.g1Func_dr(BranchLength) + _gFuncStart0.g2Func_dr(BranchLength);
	}
	if(startState == 1)
		return _gFuncStart1.g1Func_dr(BranchLength) + _gFuncStart1.g2Func_dr(BranchLength);
	else
		return 0;
}







/********************************************************************************************
gFunc
*********************************************************************************************/
computeJumps::gFunc::gFunc(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r)
: _Lambda1(Lambda1), _Lambda2(Lambda2), _r(r)
{
	_delta = sqrt((_Lambda1+_Lambda2)*(_Lambda1+_Lambda2) + 4*(_r*_r - 1)*_Lambda1*_Lambda2);
	_delta_dr = (4*_r*_Lambda1*_Lambda2)/_delta;

	_Alpha1 = 0.5*(-_Lambda1-_Lambda2 +_delta);
	_Alpha2 = 0.5*(-_Lambda1-_Lambda2 -_delta);

	_Alpha1_dr =  0.5*_delta_dr;
	_Alpha2_dr = -0.5*_delta_dr;

	_Alpha1_2    = _delta;		//= _Alpha1 - _Alpha2;
	_Alpha1_2_dr = _delta_dr;	//= _Alpha1_dr - _Alpha2_dr;

	_g1Part = ( (_r-1)*_Lambda1 - _Alpha2)/_Alpha1_2;
	_g2Part = (-(_r-1)*_Lambda1 + _Alpha1)/_Alpha1_2;

	_g1Part_dr = ( _Alpha1_2*( _Lambda1-_Alpha2_dr) - ( (_r-1)*_Lambda1 - _Alpha2)*_Alpha1_2_dr )/(_Alpha1_2*_Alpha1_2);
	_g2Part_dr = ( _Alpha1_2*(-_Lambda1+_Alpha1_dr) - (-(_r-1)*_Lambda1 + _Alpha1)*_Alpha1_2_dr )/(_Alpha1_2*_Alpha1_2);

}
//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gFunc::gFunc_dr(MDOUBLE BranchLength){
	return sign(_r)*(g1Func_dr(BranchLength) + g2Func_dr(BranchLength));
}
MDOUBLE computeJumps::gFunc::g1Func_dr(MDOUBLE BranchLength){
	return _g1Part_dr*g1Exp(BranchLength) + _g1Part*g1Exp(BranchLength)*BranchLength*_Alpha1_dr;
}
MDOUBLE computeJumps::gFunc::g2Func_dr(MDOUBLE BranchLength){
	return _g2Part_dr*g2Exp(BranchLength) +  _g2Part*g2Exp(BranchLength)*BranchLength*_Alpha2_dr;
}

//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gFunc::g1Exp(MDOUBLE BranchLength){
	return exp(_Alpha1*BranchLength);
}
MDOUBLE computeJumps::gFunc::g2Exp(MDOUBLE BranchLength){
	return exp(_Alpha2*BranchLength);
}

MDOUBLE computeJumps::gFunc::gFunc_(MDOUBLE BranchLength){
	return _g1Part*g1Exp(BranchLength) + _g2Part*g2Exp(BranchLength);
};

MDOUBLE computeJumps::gFunc::_A_(int k, int i){return BinomialCoeff((k+i-1),i) * pow(-1.0,i)*pow(_Lambda1,k)*pow(_Lambda2,(k-1)) / pow((_Lambda2-_Lambda1),(k+i))  ;	}
MDOUBLE computeJumps::gFunc::_B_(int k, int i){return BinomialCoeff((k+i-1),i) * pow(-1.0,i)*pow(_Lambda1,k)*pow(_Lambda2,(k-1)) / pow((_Lambda1-_Lambda2),(k+i))  ;	}
MDOUBLE computeJumps::gFunc::_C_(int k, int i){return BinomialCoeff((k+i-1),i) * pow(-1.0,i)*pow(_Lambda1,k)*pow(_Lambda2,(k))   / pow((_Lambda2-_Lambda1),(k+i))  ;	}
MDOUBLE computeJumps::gFunc::_D_(int k, int i){return BinomialCoeff((k+i),i)   * pow(-1.0,i)*pow(_Lambda1,k)*pow(_Lambda2,(k))   / pow((_Lambda1-_Lambda2),(k+i+1));	}

// prob for (2k-1) transitions (gains and losses), given start=0
MDOUBLE computeJumps::gFunc::qFunc_2k_1  (MDOUBLE BranchLength, int k){
	MDOUBLE qSUM = 0.0;
	for(int i=1; i<=k; ++i){
		qSUM += _A_(k,(k-i))* pow(BranchLength,(i-1))/factorial(i-1) * exp(-_Lambda1*BranchLength)
			+ _B_(k,(k-i))* pow(BranchLength,(i-1))/factorial(i-1) * exp(-_Lambda2*BranchLength);
	}
	return qSUM;
}
// prob for (2k) transitions (gains and losses), given start=0
MDOUBLE computeJumps::gFunc::qFunc_2k  (MDOUBLE BranchLength, int k){
	MDOUBLE qSUM = 0.0;
	for(int i=1; i<=(k+1); ++i){
		qSUM += _C_(k,(k-i+1))* pow(BranchLength,(i-1))/factorial(i-1)*exp(-_Lambda1*BranchLength);
	}
	for(int i=1; i<=k; ++i){
		qSUM +=  _D_(k,(k-i))* pow(BranchLength,(i-1))/factorial(i-1)*exp(-_Lambda2*BranchLength);
	}
	return qSUM;
}






/********************************************************************************************
Pij_t - Based on Analytic solution
*********************************************************************************************/
MDOUBLE computeJumps::Pij_t(const int i,const int j, const MDOUBLE d)  {
	MDOUBLE gain = _Lambda1;
	MDOUBLE loss = _Lambda2;
	MDOUBLE eigenvalue =  -(gain + loss);


	VVdouble Pt;
	int AlphaSize = 2;
	resizeMatrix(Pt,AlphaSize,AlphaSize);
	int caseNum = i + j*2;
	switch (caseNum) {
		case 0 : Pt[0][0] =  loss/(-eigenvalue) + exp(eigenvalue*d)*(1 - loss/(-eigenvalue)); break;
		case 1 : Pt[1][0] =  loss/(-eigenvalue) - exp(eigenvalue*d)*(1 - gain/(-eigenvalue)); break;
		case 2 : Pt[0][1] =  gain/(-eigenvalue) - exp(eigenvalue*d)*(1 - loss/(-eigenvalue)); break;
		case 3 : Pt[1][1] =  gain/(-eigenvalue) + exp(eigenvalue*d)*(1 - gain/(-eigenvalue));  break;
	}
	MDOUBLE val = (Pt[i][j]);
	return val; 
}
