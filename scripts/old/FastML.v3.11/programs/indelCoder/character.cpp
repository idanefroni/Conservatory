#include "character.h"
#include "gaps.h"

void character::checkForTriangleInequality(int st1, int st2){
	int longestGapStIndex = getLongestGapStIndex();
	if(_stepmatrix[st1][st2] > _stepmatrix[st1][longestGapStIndex]+_stepmatrix[st2][longestGapStIndex]){
		_isTriangleInequalityCorrectionNeeded = true;
		_stepmatrixTriagleInCorrected = _stepmatrix;
		++_stepmatrixTriagleInCorrected[st1][ longestGapStIndex];
		++_stepmatrixTriagleInCorrected[longestGapStIndex][st1];
		++_stepmatrixTriagleInCorrected[st2][ longestGapStIndex];
		++_stepmatrixTriagleInCorrected[longestGapStIndex][st2];
	}
};

int character::getLongestGapStIndex(){
	int longestGapStIndex;
	int longestGapNumOfZeros = 0;
	int characterLength = _states[0].size();
	for(int st = 0; st<_states.size(); ++st){
		int gapNumOfZeros = 0;
		for(int ind = 0; ind<characterLength; ++ind){
			gapNumOfZeros +=  _states[st][ind];
		}
		if(gapNumOfZeros > longestGapNumOfZeros)
			longestGapStIndex = st;
	}
	return longestGapStIndex;
};

	//********************************************************************************************
	//computeNumOfSteps
	// Foreach c in character_1:character_M
	//  Foreach st_x and st_y in state_0:state_c_ST (There are ST states in character c) (go over all state combinations)
	//   Do A to E steps for the pair st_x and st_y:
	//   A) translate into 01 to X set of 5'-3' coordinats of the X gaps within st_x and st_y
	//   B) ignore 0-0 (cost_c_x_y =- #0-0 colomns)
	//   C) merge adjacent 0-1 and 1-0 ((cost_c_x_y =- #adjacent 0-1 and 1-0 colomns)
	//   D) ignore 1-1 (cost_c_x_y =- #1-1 colomns)
	//********************************************************************************************
int character::computeNumOfSteps(int st1, int st2){
	int numOfSteps =_states[st1].size();
	vector<int> state1(_states[st1].size());
	state1 = _states[st1];
	vector<int> state2(_states[st2].size());
	state2 = _states[st2];

	vector<int>::iterator iter1 = state1.begin();
	vector<int>::iterator iter2 = state2.begin();
	vector<int>::iterator iterLastCounted1 = iter1;
	vector<int>::iterator iterLastCounted2 = iter2;
	
	LOGnOUT(6,<<" step "<<st1<<" "<<st2<<endl);	// DEBUG
	int i = 0;
	while( iter1!=state1.end() ){ // both same length
		if(*iter1 == *iter2  
			|| (iter1 != state1.begin() && *iter1 == *(iter1-1) && *iter2 == *(iter2-1)) 
			|| (i>0 && *iter1 == *iterLastCounted1 && *iter2 == *iterLastCounted2 && *(iter1-1)==0 )
			)
		{
			LOGnOUT(6,<<i<<" "<<*iter1<<" "<<*iter2<<endl);	// DEBUG
			//state1.erase(iter1);
			//state2.erase(iter2);
			--numOfSteps;
		}
		else{
			iterLastCounted1 = iter1;
			iterLastCounted2 = iter2;
			LOGnOUT(6,<<"Count step "<<i<<" "<<*iter1<<" "<<*iter2<<endl);	// DEBUG
		}
		++iter1;
		++iter2;
		++i;
	}

	if(state1.size() != state2.size())
		cout<<"error"<<endl;

	//numOfSteps = state1.size();
	return numOfSteps;
};

