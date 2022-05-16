#include "oneTwoMoreModel.h"
#include "matrixUtils.h"
#include "someUtil.h"

///////////////////////////////////////////////////////////
//non reversible model
///////////////////////////////////////////////////////////

const MDOUBLE EPSILON_3STATEMODEL = 1e-04;


oneTwoMoreModel::oneTwoMoreModel(const MDOUBLE m1, const MDOUBLE m2,const MDOUBLE m3, const MDOUBLE m4
								,const Vdouble &freq, bool useMarkovLimiting)
	:_gain(m1),_more(m2), _less(m3),_loss(m4),_freq(freq),_useMarkovLimiting(useMarkovLimiting){
		resizeMatrix(_Q,alphabetSize(),alphabetSize());
		resizeMatrix(_lastPtCalculated, alphabetSize(), alphabetSize());
		updateQ();
	}

oneTwoMoreModel& oneTwoMoreModel::operator=(const oneTwoMoreModel &other){
	_gain = other._gain;
	_more = other._more;
	_less = other._less;
	_loss = other._loss;
	_freq = other._freq;
	_useMarkovLimiting = other._useMarkovLimiting;
	_Q = other._Q;
	_bQchanged = other._bQchanged;
	_lastPtCalculated = other._lastPtCalculated;
	_lastTcalculated = other._lastTcalculated;

	return *this;
}

void oneTwoMoreModel::updateQ(){
	setEpsilonForZeroParams();
	_Q[0][0] = -_gain;
	_Q[0][1] = _gain;
	_Q[0][2] = 0;
	_Q[1][0] = _loss;
	_Q[1][1] = -_more-_loss;
	_Q[1][2] = _more;
	_Q[2][0] = 0;
	_Q[2][1] = _less;
	_Q[2][2] = -_less;
	for (int i=0; i<_Q.size();i++) {
		MDOUBLE sum = _Q[i][0]+_Q[i][1]+_Q[i][2];
		if ((abs(sum)>err_allow_for_pijt_function()))
			errorMsg::reportError("Error in oneTwoMoreModel::updateQ, sum of row is not 0");
	}
	if ((!checkIsNullModel()) && (_useMarkovLimiting))
		computeMarkovLimitingDistribution();
	_bQchanged = true;
}

// when Q matrix parameters are zero the lib code underflows and the likelihood is set to EPSILON
void oneTwoMoreModel::setEpsilonForZeroParams(){
	if (DEQUAL(_more,0.0,EPSILON_3STATEMODEL))
		_more = EPSILON_3STATEMODEL;
	if (DEQUAL(_gain,0.0,EPSILON_3STATEMODEL))
		_gain = EPSILON_3STATEMODEL;
	if (DEQUAL(_loss,0.0,EPSILON_3STATEMODEL))
		_loss = EPSILON_3STATEMODEL;
	if (DEQUAL(_less,0.0,EPSILON_3STATEMODEL))
		_less = EPSILON_3STATEMODEL;
}

void oneTwoMoreModel::setMu1(const MDOUBLE val) { 
	_gain = val; 
	updateQ();
}

void oneTwoMoreModel::setMu2(const MDOUBLE val) { 
	_more = val; 
	updateQ();
}

void oneTwoMoreModel::setMu3(const MDOUBLE val) { 
	_less = val; 
	updateQ();
}

void oneTwoMoreModel::setMu4(const MDOUBLE val) { 
	_loss = val; 
	updateQ();
}



bool oneTwoMoreModel::pijt_is_prob_value(MDOUBLE val) const {
	if ((abs(val)+err_allow_for_pijt_function()<0) || (val>1+err_allow_for_pijt_function()))
		return false;
	else
		return true;
}

bool oneTwoMoreModel::areFreqsValid(Vdouble freq) const{
	MDOUBLE sum=0.0;
	for (int i=0; i<freq.size(); ++i){
		if (freq[i]<0.0)
			return false;
		sum+=freq[i];
	}
	if (!DEQUAL(sum,1.0)) {
		return false;
	}
	return true;
}

bool oneTwoMoreModel::checkIsNullModel(){
	if (_more!=EPSILON_3STATEMODEL)
		return false;
	if (_more!=EPSILON_3STATEMODEL)
		return false;
	if (!(DEQUAL(_freq[2],1.0,EPSILON_3STATEMODEL)))
		return false;
	return true;
}

void oneTwoMoreModel::setFreq(const Vdouble &freq){
	if (freq.size()!=_freq.size()) {
		errorMsg::reportError("Error in oneTwoMoreModel::setFreq, size of freq is different than member");
	}
	
	if (!areFreqsValid(freq)) {
		string strErr = "Error in oneTwoMoreModel::setFreq, sum of freq is different than 1 or negative freq value";
		errorMsg::reportError(strErr);
	}
	for (int i=0; i<freq.size(); ++i){
		_freq[i] = freq[i];
	}
}






void oneTwoMoreModel::computeMarkovLimitingDistribution(){

	VVdouble P;
	int as = alphabetSize();
	resizeMatrix(P,as, as);
	// initializing P with P at time 1
	for (int i=0; i< as; ++i) {
		for (int j=0; j< as; ++j) {
			P[i][j]=Pij_t(i,j,1.0);
		}
	}
	VVdouble previous_P = P;
	int numIterations = 0;
	Vdouble freqs(3,-1.0);
	bool converged = false;
	MDOUBLE epsilon=0.000001;
	int row, col;
	
	while ( converged==false ) {
		previous_P = P;
		P = multiplyMatrixes(P,P);
		// due to rounding errors, we set the diagonal to be 1-(the rest)	
		P[0][0]=1.0-P[0][1]-P[0][2];
		P[1][1]=1.0-P[1][0]-P[1][2];
		P[2][2]=1.0-P[2][0]-P[2][1];
		for (int d=0; d<as;++d){
			freqs[d] = P[0][d];// ** taking the freqs as the first row; this is not necessarily correct if 3 rows are different
		}
		converged = true;
		for (row = 0; row < P.size(); ++row) {
			for (col = 0; col < P.size(); ++col)
			{
				MDOUBLE diff = abs(convert(previous_P[row][col] - P[row][col]));
				if ( ( ( ( !DEQUAL(diff,0.0,epsilon) ) || (!areFreqsValid(freqs) ) ) )){
					converged = false;
				}
			}
		}
		numIterations++;
		if (numIterations>100) { 
			string err = "Error in oneTwoMoreModel::computeMarkovLimitingDistribution, too many iterations =" + double2string(numIterations);
			errorMsg::reportError(err);
		}

	}
//making sure that the three rows are the same
	for (row =1; row < P.size(); ++row) {
	    for (col = 0; col < P.size(); ++col)
	    {
		if (!(DEQUAL(P[row][col],P[row-1][col],epsilon))) {
		    errorMsg::reportError("Error in oneTwoMoreModel::computeMarkovLimitingDistribution, rows are not equal"    );
		    
		}
		
	    }
	    
	}
		
	setFreq(freqs);
}

// new implementation copied from Itay Mayrose which saves the last values of t computed
const MDOUBLE oneTwoMoreModel::Pij_t(const int i,const int j, const MDOUBLE d) const
{
	if (!_bQchanged && DEQUAL(d, _lastTcalculated))
		return convert(_lastPtCalculated[i][j]);
	// converting Q into doubleRep format
	VVdouble QdblRep; 
	resizeMatrix(QdblRep,_Q.size(),_Q.size());
	for (int row=0;row<_Q.size();row++){
			for (int col=0;col<_Q[row].size();col++)
				QdblRep[row][col]=convert(_Q[row][col]);
	}

	VVdouble Qt = multiplyMatrixByScalar(QdblRep, d);
	VVdouble unit;
	unitMatrix(unit,_Q.size());
	_lastPtCalculated = add(unit,Qt) ; // I + Qt
	VVdouble Qt_power = Qt;
	VVdouble prevIter_matrix = _lastPtCalculated;
	VVdouble diffM = _lastPtCalculated; //init to whatever
	int n=2;
	bool bConverged = false;
	while (bConverged == false) 
	{
		prevIter_matrix = _lastPtCalculated;
		VVdouble tempQ = multiplyMatrixByScalar(Qt,1.0/n);
		Qt_power = multiplyMatrixes(Qt_power,tempQ);
		_lastPtCalculated = add(_lastPtCalculated,Qt_power); // I + Qt + Qt^2/2! + ....  + Qt^n/n!
		//check if the difference between the cur and prev iteration is smaller than the allowed error of all matrix entries
		bConverged = true;
		for (int row = 0; row < _lastPtCalculated.size(); ++row) {
			for (int col = 0; col < _lastPtCalculated.size(); ++col)
			{
				MDOUBLE diff = abs(convert(_lastPtCalculated[row][col] - prevIter_matrix[row][col]));
				if ((diff > err_allow_for_pijt_function()) || (!pijt_is_prob_value(convert(_lastPtCalculated[i][j]))))
					bConverged = false;
			}
		}
		n++;
		if (n>150) { 
			string err = "Error in oneTwoMoreModel::Pij_t, too many iterations for t = " + double2string(d);
			//cerr<<diff<<endl;
			errorMsg::reportError(err);
		}
	}
	MDOUBLE val = convert(_lastPtCalculated[i][j]);
	if (!pijt_is_prob_value(val))
		errorMsg::reportError("Error in oneTwoMoreModel::Pij_t, pijt <0 or >1");
	if (val<0.0)
		val = EPSILON; // absolute zero creates a problem later on in computations
	if (val>1.0)
		val = 1.0;
	_bQchanged = false;
	return val; 
}
//////////////////////////////////////////////////////////////////////////
MDOUBLE oneTwoMoreModel::sumPijQij(){
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= (_Q[i][i])*_freq[i];
	}
	return sum;
}
//////////////////////////////////////////////////////////////////////////
void oneTwoMoreModel::norm(const MDOUBLE scale)
{
	for (int i=0; i < _Q.size(); ++i) {
		for (int j=0; j < _Q.size(); ++j) {
			_Q[i][j] *= scale; 		
		}
	}
}
