#ifndef ___CHARACTER__
#define ___CHARACTER__

#include "definitions.h"
#include "gaps.h"
#include "matrixUtils.h"
#include "indelCoderOptions.h"


using namespace std;


class character {
public:

	explicit character(int coord_5p, int coord_3p, int numOfSquencs, int numOfStates=0):_coord_5p(coord_5p), _coord_3p(coord_3p), _numOfSequences(numOfSquencs)
	{
		_numOfStates = 1;
		_isTriangleInequalityCorrectionNeeded = false;
		//_sc_states.resize(numOfSquencs);	// No need - done later
	};

	~character() {};
	
	int getCoord5(){return _coord_5p;}
	int getCoord3(){return _coord_3p;}	
	void setCoord3(int coord_3p){ _coord_3p = coord_3p;}
	int getNumOfGaps(){return _gaps.numOfGaps();}	
	int getNumOfStates(){return _numOfStates;}
	vector<int> getGapsIndices() const {return _gapsIndices;}


	void addGap(gaps::gap* gap_p, int gapIndex){
		_gaps.insertNewGap(gap_p);
		_gapsIndices.push_back(gapIndex);
	}
	
	void addZeroState(){
		vector<int> zeroState((int)(_coord_3p-_coord_5p+1),1);	// vector of ones, length of character
		_states.push_back(zeroState);
	};
	
	void addState(vector<int> state){
		_states.push_back(state);
		++_numOfStates;
	};

	void resizeSc_states(){resizeMatrix(_sc_states,_numOfSequences,(int)(_coord_3p-_coord_5p+1)); oneMatrix(_sc_states); };
	void resizeStepMatrix(){resizeMatrix(_stepmatrix,_numOfStates,_numOfStates); };
	void setGapsInSc_states(int seqId, int coord5, int coord3){
		for (int i = coord5; i<=coord3; ++i){
			_sc_states[seqId][i-_coord_5p] = 0;
		}
	};
	vector< vector<int> > getScStates(){return _sc_states;};
	
	vector< vector<int> > getStates(){return _states;};

	
	//********************************************************************************************
	//isTriangleInequalityCorrectionNeeded
	//********************************************************************************************
	bool isTriangleInequalityCorrectionNeeded(){return _isTriangleInequalityCorrectionNeeded;};
	void checkForTriangleInequality(int st1, int st2);
	int getLongestGapStIndex();


	int computeNumOfSteps(int st1, int st2);
	
	
	//*******************************************************************************************
	//printScStates
	//*******************************************************************************************
	void printScStates(){
		cout<<"ScStates:"<<endl;
		for(int s=0; s<_sc_states.size(); ++s){
			for(int alp=0; alp<_sc_states[0].size(); ++alp){
				cout<<_sc_states[s][alp];
			}
			cout<<endl;
		}			
	}
	
	//********************************************************************************************
	//printStates
	//*********************************************************************************************
	void printStates(){
		cout<<"States:"<<endl;
		for(int st=0; st<_numOfStates; ++st){
			for(int alp=0; alp<_states[0].size(); ++alp){
				cout<<_states[st][alp];
			}
			cout<<endl;
		}			
	}


	//********************************************************************************************
	//determinationStepsMatrix
	//*********************************************************************************************
	void determinationStepsMatrix(){		
		resizeStepMatrix();
		for(int st1 = 0; st1< _numOfStates; ++st1){
			for(int st2= st1; st2< _numOfStates; ++st2){
				if(st1==st2)
					_stepmatrix[st1][st2] = 0;
				else{
					_stepmatrix[st1][st2] = computeNumOfSteps(st1,st2);
					_stepmatrix[st2][st1] = _stepmatrix[st1][st2];
				}
				if(indelCoderOptions::_isCheckForTriangleInequality)
					checkForTriangleInequality(st1,st2);
			}
		}		
	}
	
	//********************************************************************************************
	//printStepsMatrix
	//*********************************************************************************************
	void printStepsMatrix(ostream& out = cout, bool isCorrectdForTriangleInEq = false){
		out<<"    ";
		for(int st1 = 0; st1< _numOfStates; ++st1){
			out<<st1<<" ";
		}
		out<<endl;

		for(int st1 = 0; st1< _numOfStates; ++st1){
			out<<"["<<st1<<"] ";
			for(int st2= 0; st2< _numOfStates; ++st2){
				if(isCorrectdForTriangleInEq)
					out<<_stepmatrixTriagleInCorrected[st1][st2]<<" ";
				else
                    out<<_stepmatrix[st1][st2]<<" ";
			}
			out<<endl;
		}		
	}




private:
	int _coord_5p; 
	int _coord_3p;	
	int _numOfStates;
	int _numOfSequences;
	
	gaps _gaps;							// gaps included in this character	
	vector<int> _gapsIndices;			// since all gaps are indexed, here you find the indices of the gaps included in this character
	vector< vector<int> > _stepmatrix;
	vector< vector<int> > _sc_states;	// matrix - species X lengthOfCharacter
	vector< vector<int> > _states;

	bool _isTriangleInequalityCorrectionNeeded;
	vector< vector<int> > _stepmatrixTriagleInCorrected;

};

#endif
