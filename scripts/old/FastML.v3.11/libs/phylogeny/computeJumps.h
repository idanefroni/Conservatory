#ifndef ___COMPUTE_JUMPS__
#define ___COMPUTE_JUMPS__

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "alphabet.h"
#include "someUtil.h"
#include <math.h>

#include <map>
#include <vector>
using namespace std;

/******************************************************************
This class compute jumps (events) by Suchard equations along differing branch lengths (according to a 
given tree), with the aim of giving the expectation of the number of jumps
from state a to state b given that the terminal states at the end of the branch are
x and y.
*******************************************************************/

class computeJumps  {
public:
	computeJumps(const MDOUBLE Lambda1, const MDOUBLE Lambda2, const MDOUBLE r=1, const int maxNumOfChangesPerBranchSum=5);
	virtual ~computeJumps();

	/******************************************************************
	Foreach computeJumps, for gFunc objects are needed:
	inner class gFunc, if startState=0, Lambda1=gain, Lambda2= loss 
	if startState=1, Lambda1=loss, Lambda2= gain.
	For both with use +r and -r versions
	*******************************************************************/
	class gFunc  {
	public:
		gFunc(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r);
		gFunc(){};
		~gFunc(){};

		MDOUBLE gFunc_dr(MDOUBLE BranchLength);
		MDOUBLE g1Func_dr(MDOUBLE BranchLength);
		MDOUBLE g2Func_dr(MDOUBLE BranchLength);
		MDOUBLE g1Exp(MDOUBLE BranchLength);
		MDOUBLE g2Exp(MDOUBLE BranchLength);
		MDOUBLE gFunc_(MDOUBLE BranchLength);

		//////////////////////////////////////////////////////////////////////////
		MDOUBLE _A_(int k, int i);
		MDOUBLE _B_(int k, int i);
		MDOUBLE _C_(int k, int i);
		MDOUBLE _D_(int k, int i);
		// prob for (2k-1) transitions (gains and losses), given start=0
		MDOUBLE qFunc_2k_1  (MDOUBLE BranchLength, int k=1);
		// prob for (2k) transitions (gains and losses), given start=0
		MDOUBLE qFunc_2k  (MDOUBLE BranchLength, int k=0);

	private:
		MDOUBLE _r;
		MDOUBLE _Lambda1;
		MDOUBLE _Lambda2;

		MDOUBLE _Alpha1;
		MDOUBLE _Alpha2;
		MDOUBLE _Alpha1_dr;
		MDOUBLE _Alpha2_dr;

		MDOUBLE _Alpha1_2;
		MDOUBLE _Alpha1_2_dr;

		MDOUBLE _delta;
		MDOUBLE _delta_dr;

		MDOUBLE _g1Part;
		MDOUBLE _g2Part;
		MDOUBLE _g1Part_dr;
		MDOUBLE _g2Part_dr;

	};
	//////////////////////////////////////////////////////////////////////////
	
	MDOUBLE getExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId);
	MDOUBLE getTotalExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd);

	MDOUBLE gainExp(MDOUBLE BranchLength,MDOUBLE prob01,MDOUBLE prob11);
	
	MDOUBLE gainExpGiven01(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven00(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven11(MDOUBLE BranchLength);
	MDOUBLE gainExpGiven10(MDOUBLE BranchLength);

	MDOUBLE lossExpGiven01(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven00(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven11(MDOUBLE BranchLength);
	MDOUBLE lossExpGiven10(MDOUBLE BranchLength);

	MDOUBLE getProb(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId);
	MDOUBLE gainProbGiven01(MDOUBLE BranchLength);
	MDOUBLE gainProbGiven00(MDOUBLE BranchLength);
	MDOUBLE gainProbGiven11(MDOUBLE BranchLength);
	MDOUBLE gainProbGiven10(MDOUBLE BranchLength);

	MDOUBLE lossProbGiven01(MDOUBLE BranchLength);
	MDOUBLE lossProbGiven00(MDOUBLE BranchLength);
	MDOUBLE lossProbGiven11(MDOUBLE BranchLength);
	MDOUBLE lossProbGiven10(MDOUBLE BranchLength);


	MDOUBLE gFunc_dr(MDOUBLE BranchLength, int startState);
    	
private:
	MDOUBLE m01(MDOUBLE BranchLength);
	MDOUBLE m00(MDOUBLE BranchLength);
	MDOUBLE m11(MDOUBLE BranchLength);
	MDOUBLE m10(MDOUBLE BranchLength);	

	MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d);

	MDOUBLE _Lambda1;
	MDOUBLE _Lambda2;
	int _maxNumOfChangesPerBranchSum;

	gFunc _gFuncStart0;
	gFunc _gFuncStart0MinusR;
	gFunc _gFuncStart1;
	gFunc _gFuncStart1MinusR;
};

#endif
