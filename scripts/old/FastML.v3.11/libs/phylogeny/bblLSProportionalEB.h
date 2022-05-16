#ifndef ___R4SP_BBL_LS
#define ___R4SP_BBL_LS

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"
#include "likelihoodComputation.h"
#include <vector>
using namespace std;

#define MAX_BRANCH_LENGTH 10.0

/*
This class optimize the branches using "naive" line search methodology.
go over each branch and optimize it using brent.
In one iteration it optimze seperatly all branches.
This procedure continues until convergence is reached or until the maximum number of iteration is reached.
*/
class bblLSProportionalEB {
public:
	
	explicit bblLSProportionalEB(tree& et, const vector<sequenceContainer>& sc, multipleStochasticProcess* msp, const gammaDistribution* pProportionDist, Vdouble& treeLikelihoodVec, const bool optimizeSelectedBranches=false, int maxIter=50, MDOUBLE epsilon=0.05);
	~bblLSProportionalEB() {};
	Vdouble getTreeLikelihoodVec() const {return _treeLikelihoodVec;}
private:
	Vdouble optimizeBranches(tree& et, const vector<sequenceContainer>& sc, multipleStochasticProcess* msp, const gammaDistribution* pProportionDist, Vdouble& treeLikelihoodVec, const bool optimizeSelectedBranches=false, int maxIter=50, MDOUBLE epsilon=0.05);

private:
	Vdouble _treeLikelihoodVec;
};

class evalR4SPBranch{
public:
	explicit evalR4SPBranch(tree::nodeP pNode, tree& et, const vector<sequenceContainer>& sc, multipleStochasticProcess* msp, const gammaDistribution* pProportionDist)
		:_pNode(pNode),_et(et), _sc(sc), _msp(msp), _pProportionDist(pProportionDist){};
private:
	tree::nodeP _pNode;
	tree& _et;
	const vector<sequenceContainer>& _sc;
	multipleStochasticProcess* _msp;
	const gammaDistribution* _pProportionDist;
public:
	MDOUBLE operator() (MDOUBLE bl) {
		_pNode->setDisToFather(bl);
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,_sc,_msp,_pProportionDist);
		MDOUBLE res = sumVdouble(likeVec);
		return -res;
	}
};

#endif
