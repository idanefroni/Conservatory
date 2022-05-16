// $Id: bestTamura92param.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___BEST_TAMURA92_PARAM
#define ___BEST_TAMURA92_PARAM

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "tamura92.h"


class bestTamura92ParamFixedTree {
public:
	explicit bestTamura92ParamFixedTree(const tree& et, // find best TrTv and theta
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights,
										const int maxTotalIterations = 5,
										const MDOUBLE epsilonLikelihoodImprovment = 0.05,
										const MDOUBLE epsilonLoglikelihoodForTrTvOptimization = 0.01,
										const MDOUBLE epsilonLoglikelihoodForThetaOptimization = 0.01,
										const MDOUBLE upperBoundOnTrTv = 5.0);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
};

class bestTamura92ParamAndBBL{
public:	
	explicit bestTamura92ParamAndBBL(tree& et, //find best TrTv, theta and best BBL
												 const sequenceContainer& sc,
												 stochasticProcess& sp,
												 const Vdouble * weights=NULL,
												 const int maxTotalIterations=5,
												 const MDOUBLE epsilonLikelihoodImprovment=0.05,
												 const MDOUBLE epsilonLoglikelihoodForTrTvOptimization=0.01,
												 const MDOUBLE epsilonLoglikelihoodForThetaOptimization=0.01,
												 const MDOUBLE epsilonLoglikelihoodForBBL=0.01,
												 const MDOUBLE upperBoundOnTrTv=5.0,
												 const int maxBBLIterations=10);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta(int spIndex) {return _bestTheta;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestL;
};


class bestTamura92ParamAlphaAndBBL {
public:
	explicit bestTamura92ParamAlphaAndBBL( //find best TrTv, theta, Alpha and best branch lengths
		tree& et,
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const int maxTotalIterations=5,
		const MDOUBLE epsilonLikelihoodImprovment= 0.05,
		const MDOUBLE epsilonLoglikelihoodForTrTvOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForThetaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.01,
		const MDOUBLE upperBoundOnTrTv = 5.0,
		const int maxBBLIterations=10,
		const MDOUBLE initAlpha = 1.5,
		const MDOUBLE upperBoundOnAlpha = 5.0);
	MDOUBLE getBestTrTv() {return _bestTrTv;}
	MDOUBLE getBestTheta() {return _bestTheta;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestTrTv;
	MDOUBLE _bestTheta;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestTamura92ParamAlphaAndBBLProportional {
public:
	explicit bestTamura92ParamAlphaAndBBLProportional( //find best TrTv, theta, loca Alpha for each gene, global Alpha and best branch lengths
		tree& et,
		vector<sequenceContainer>& sc,
		multipleStochasticProcess* msp,
 	    gammaDistribution* pProportionDist,
		Vdouble initLocalAlphas,
		Vdouble initLocalKappas,
		Vdouble initLocalThetas,
		const MDOUBLE upperBoundOnLocalAlpha,
		const MDOUBLE initGlobalAlpha,
		const MDOUBLE upperBoundOnGlobalAlpha,
		const MDOUBLE upperBoundOnTrTv,
		const int maxTotalIterations,
		const int maxBBLIterations,
		const bool optimizeSelectedBranches=false,
		const bool optimizeTree = true,
		const string branchLengthOptimizationMethod="bblLS",
		const bool optimizeLocalParams = true,
		const bool optimizeGlobalAlpha = true,
		const Vdouble * weights=NULL,
		const MDOUBLE epsilonLikelihoodImprovment= 0.05,
		const MDOUBLE epsilonLoglikelihoodForLocalTrTvOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForLocalThetaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForLocalAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForGlobalAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.01);
	MDOUBLE getBestTrTv(int spIndex) {return _bestTrTvVec[spIndex];}
	MDOUBLE getBestTheta(int spIndex) {return _bestThetaVec[spIndex];}
	MDOUBLE getBestLocalAlpha(int spIndex) {return _bestLocalAlphaVec[spIndex];}
	MDOUBLE getBestGlobalAlpha() {return _bestGlobalAlpha;}
	Vdouble getBestL() {return _bestLvec;}
private:
	Vdouble _bestTrTvVec;
	Vdouble _bestThetaVec;
	Vdouble _bestLocalAlphaVec;
	MDOUBLE _bestGlobalAlpha;
	Vdouble _bestLvec;
};


class C_evalTrTvParam{
public:
  C_evalTrTvParam( const tree& et,
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
	MDOUBLE operator() (MDOUBLE TrTv) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(TrTv);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with TrTv = "<<TrTv<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalLocalTrTvParam{
public:
  C_evalLocalTrTvParam( const tree& et,
				   const sequenceContainer& sc,
				   stochasticProcess& sp,
				   gammaDistribution* pProportionDist,
				   const Vdouble * weights = NULL)
	  : _et(et),_sc(sc),_weights(weights),_sp(sp),_pProportionDist(pProportionDist){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
	gammaDistribution* _pProportionDist;
public:
	MDOUBLE operator() (MDOUBLE TrTv) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(TrTv);
		vector<sequenceContainer> tmpScVec;
		tmpScVec.push_back(_sc);
		vector<stochasticProcess> tmpSpVec;
		tmpSpVec.push_back(_sp);
		multipleStochasticProcess * tmpMsp = new multipleStochasticProcess();
		tmpMsp->setSpVec(tmpSpVec);	
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,tmpScVec,tmpMsp,_pProportionDist);
		MDOUBLE res = likeVec[0];
		delete(tmpMsp);
		LOG(5,<<" with TrTv = "<<TrTv<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalLocalTheta{
public:
  C_evalLocalTheta( const tree& et,
				   const sequenceContainer& sc,
				   stochasticProcess& sp,
				   gammaDistribution* pProportionDist,
				   const Vdouble * weights = NULL)
	  : _et(et),_sc(sc),_weights(weights),_sp(sp),_pProportionDist(pProportionDist){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
	gammaDistribution* _pProportionDist;
public:
	MDOUBLE operator() (MDOUBLE theta) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTheta(theta);
		vector<sequenceContainer> tmpScVec;
		tmpScVec.push_back(_sc);
		vector<stochasticProcess> tmpSpVec;
		tmpSpVec.push_back(_sp);
		multipleStochasticProcess * tmpMsp = new multipleStochasticProcess();
		tmpMsp->setSpVec(tmpSpVec);	
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,tmpScVec,tmpMsp,_pProportionDist);
		MDOUBLE res = likeVec[0];
		delete(tmpMsp);
		LOG(5,<<" with Theta = "<<theta<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalTheta{
public:
  C_evalTheta(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE theta) {
		(static_cast<tamura92*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTheta(theta);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with theta = "<<theta<<" logL = "<<res<<endl);
		return -res;
	}
};

#endif


