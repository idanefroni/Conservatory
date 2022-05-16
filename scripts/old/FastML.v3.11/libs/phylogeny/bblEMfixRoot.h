// $Id: bblEM.h 4478 2008-07-17 17:09:55Z cohenofi $
#ifndef ___BBL_EM_GL__FIXED_ROOT
#define ___BBL_EM_GL__FIXED_ROOT

/********************************************************************************************
Class::bblEM (with variation: bblEMfixRoot, bblEM2codon)
compute_bblEM
allocatePlace (one more level for fixRoot - in computeDownAlg and countsTableVec)
bblEM_it (called at each iteration of BBL)
foreach pos{
	computeDown (use variants for fix root - fillComputeDownNonReversible
	   vector<suffStatGlobalGamPos> _cdown;		//_cdown[categ][letter@root][nodeid][letter][prob])
	addCounts
	addCountsFixedRoot (based on computeUp and computeDown... fill _computeCountsV)
	use class::computeCounts (but no duplicated class!!!)
}
optimizeBranches
foreach node{
		class::fromCountTableComponentToDistance (with variation: ...fixRoot, ...2Codon)
		computeDistance() + set - based on 
				class::likeDist (with variation: ...fixRoot, ...2Codon)
				giveDistance()
				giveDistanceBrent()
				C_evallikeDist and C_evallikeDist_d
				.... computation based on counts{alph1,alph2, root, rate(sp)}:      sumL+= _ctc.getCounts(alph1,alph2,rateCategor)*(log(   _sp.Pij_t(alph1,alph2,dist*rate)    )-log(_sp.freq(alph2)))
				

}
*********************************************************************************************/


#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "countTableComponent.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include "gainLossAlphabet.h"
#include "unObservableData.h"
#include <vector>

using namespace std;

class bblEMfixRoot {
public:
	explicit bblEMfixRoot(tree& et,
				const sequenceContainer& sc,
				const stochasticProcess& sp,
				const Vdouble * weights = NULL,
				const int maxIterations=50,
				const MDOUBLE epsilon=0.05,
				const MDOUBLE tollForPairwiseDist=0.001,
				unObservableData*  _unObservableData_p=NULL,
				const MDOUBLE* likelihoodLast=NULL);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEM(const int maxIterations,
					const MDOUBLE epsilon,
					const MDOUBLE tollForPairwiseDist,
					const MDOUBLE* likelihoodLast=NULL);
	void bblEM_it(const MDOUBLE tollForPairwiseDist);
	void computeDown(const int pos);
	void computeUp();
	void addCounts(const int pos);
	void addCountsFixedRoot(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig);


	void optimizeBranches(const MDOUBLE tollForPairwiseDist);
	void allocatePlace();



	MDOUBLE _treeLikelihood;
	tree& _et;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	vector< vector< countTableComponentGam > > _computeCountsV; // _computeCountsV [node] [letter@root] [rate][alph][alph]
	computePijGam _pij;
	suffStatGlobalGam _cup;						//_cup[pos][categ][nodeid][letter][prob]			
	//suffStatGlobalGamPos _cdown;				// for each pos: computeDown(pos);	addCounts(pos);  
	vector<suffStatGlobalGamPos> _cdown;		//_cdown[categ][letter@root][nodeid][letter][prob] - since fillComputeDownNonReversible uses this assumption
	const Vdouble * _weights;
	VdoubleRep _posLike;
	unObservableData*  _unObservableData_p;
};

#endif
