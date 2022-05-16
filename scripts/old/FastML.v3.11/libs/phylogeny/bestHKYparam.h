// $Id: bestHKYparam.h 9992 2011-11-08 03:57:29Z rubi $

#ifndef ___BEST_HKY_PARAM
#define ___BEST_HKY_PARAM

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "hky.h"
#include "multipleStochasticProcess.h"

class bestHkyParamFixedTree {
public:
	explicit bestHkyParamFixedTree(const tree& et,
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnHkyParam = 0.5,
					   const MDOUBLE epsilonHkyParamOptimization = 0.01);
	MDOUBLE getBestHkyParam() {return _bestHkyParam;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestHkyParam;
	MDOUBLE _bestL;
};

class bestHkyParamAndBBL {
public:
	explicit bestHkyParamAndBBL(tree& et, //find Best HkyParam and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnHkyParam = 5.0,
					   const MDOUBLE epsilonHkyParamOptimization= 0.01,
					   const MDOUBLE epsilonLikelihoodImprovment= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestHkyParam() {return _bestHkyParam;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestHkyParam;
	MDOUBLE _bestL;
};

class C_evalHkyParam{
public:
  C_evalHkyParam(	const tree& et,
				const sequenceContainer& sc,
				stochasticProcess& sp,
				const Vdouble * weights = NULL)
    : _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE HkyParam) {
		(static_cast<hky*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(HkyParam);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with HkyParam = "<<HkyParam<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalLocalHkyParam{
public:
  C_evalLocalHkyParam(	const tree& et,
				const sequenceContainer& sc,
				stochasticProcess& sp,
				const gammaDistribution* pProportionDist,
				const Vdouble * weights = NULL)
    : _et(et),_sc(sc),_weights(weights),_sp(sp),_pProportionDist(pProportionDist){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
	const gammaDistribution* _pProportionDist;
public:
	MDOUBLE operator() (MDOUBLE HkyParam) {
		(static_cast<hky*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(HkyParam);
		vector<sequenceContainer> tmpScVec;
		tmpScVec.push_back(_sc);
		vector<stochasticProcess> tmpSpVec;
		tmpSpVec.push_back(_sp);
		multipleStochasticProcess * tmpMsp = new multipleStochasticProcess();
		tmpMsp->setSpVec(tmpSpVec);	
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,tmpScVec,tmpMsp,_pProportionDist);
		MDOUBLE res = likeVec[0];
		delete(tmpMsp);
		LOG(5,<<" with HkyParam = "<<HkyParam<<" logL = "<<res<<endl);
		return -res;
	}
};

class bestHkyParamAlphaAndBBL {
public:
	explicit bestHkyParamAlphaAndBBL( //find best TrTv (=HkyParam), Alpha and best branch lengths
		tree& et,
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const int maxTotalIterations=5,
		const MDOUBLE epsilonLikelihoodImprovment= 0.05,
		const MDOUBLE epsilonHkyParamOptimization= 0.01,
		const MDOUBLE epsilonAlphaOptimization= 0.01,
		const MDOUBLE epsilonBBL= 0.01,
		const MDOUBLE upperBoundOnHkyParam = 5.0,
		const int maxBBLIterations=10,
		const MDOUBLE initAlpha = 1.5,
		const MDOUBLE upperBoundOnAlpha = 5.0);
	
	MDOUBLE getBestHkyParam() {return _bestHkyParam;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestHkyParam;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};


class bestHkyParamAlphaAndBBLProportional {
public:
	explicit bestHkyParamAlphaAndBBLProportional( //find best Kappa (=HkyParam), global Alpha, local Alpha, and best branch lengths
		tree& et,
		vector<sequenceContainer>& sc,
		multipleStochasticProcess* msp,
		gammaDistribution* pProportionDist,
		Vdouble initLocalAlphas,
		Vdouble initLocalKappas,
		const MDOUBLE upperBoundOnLocalAlpha,
		const MDOUBLE initGlobalAlpha,
		const MDOUBLE upperBoundOnGlobalAlpha,
		const MDOUBLE upperBoundOnHkyParam,
		const int maxTotalIterations,
		const int maxBBLIterations,
		const bool optimizeSelectedBranches=false,
		const bool optimizeTree = true,
		const string branchLengthOptimizationMethod="bblLS",
		const bool optimizeLocalParams = true,
		const bool optimizeGlobalAlpha = true,
		const Vdouble * weights=NULL,
		const MDOUBLE epsilonLikelihoodImprovment= 0.05,
		const MDOUBLE epsilonHkyParamOptimization= 0.01,
		const MDOUBLE epsilonLocalAlphaOptimization= 0.01,
		const MDOUBLE epsilonGlobalAlphaOptimization= 0.01,
		const MDOUBLE epsilonBBL= 0.01);
	
	MDOUBLE getBestHkyParam(int spIndex) {return _bestHkyParamVec[spIndex];}
	MDOUBLE getBestLocalAlpha(int spIndex) {return _bestLocalAlphaVec[spIndex];}
	MDOUBLE getBestGlobalAlpha(){return _bestGlobalAlpha;}
	Vdouble getBestL() {return _bestLvec;}
private:
	Vdouble _bestHkyParamVec;
	Vdouble _bestLocalAlphaVec;
	MDOUBLE _bestGlobalAlpha;
	Vdouble _bestLvec;
};





#endif


