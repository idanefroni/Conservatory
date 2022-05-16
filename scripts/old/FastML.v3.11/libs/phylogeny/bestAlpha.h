// $Id: bestAlpha.h 10000 2011-11-12 18:20:12Z rubi $

#ifndef ___BEST_ALPHA
#define ___BEST_ALPHA

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "logFile.h"

#ifndef VERBOS
#define VERBOS
#endif

class bestAlphaFixedTree {
public:
	explicit bestAlphaFixedTree(const tree& et,
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const MDOUBLE upperBoundOnAlpha = 15,
		const MDOUBLE epsilonAlphaOptimization = 0.01);
		MDOUBLE getBestAlpha() {return _bestAlpha;}
		MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestAlphaAndBBL {
public:
	explicit bestAlphaAndBBL(tree& et, //find Best Alpha and best BBL
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const MDOUBLE initAlpha = 1.5,
		const MDOUBLE upperBoundOnAlpha = 5.0,
		const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
		const int maxBBLIterations=10,
		const int maxTotalIterations=5);
		MDOUBLE getBestAlpha() {return _bestAlpha;}
		MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestAlphasAndBBLProportional {
public:
	explicit bestAlphasAndBBLProportional(tree& et, //find Best Alphas (per gene - local and proportional factors - global) and best BBL
		vector<sequenceContainer>& sc,
		multipleStochasticProcess* msp,
		gammaDistribution* pProportionDist,
		Vdouble initLocalRateAlphas,
		const MDOUBLE upperBoundOnLocalRateAlpha,
		const MDOUBLE initGlobalRateAlpha,
		const MDOUBLE upperBoundOnGlobalRateAlpha,
		const int maxBBLIterations,
		const int maxTotalIterations,
		const bool optimizeSelectedBranches=false,
		const bool optimizeTree = true,
		const string branchLengthOptimizationMethod="bblLS",
		const bool optimizeLocalAlpha = true,
		const bool optimizeGlobalAlpha = true,
		const Vdouble * weights=NULL,
		const MDOUBLE epsilonLoglikelihoodForLocalRateAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForGlobalRateAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.05);
		MDOUBLE getBestLocalAlpha(int spIndex){return _bestLocalAlphaVec[spIndex];}
		MDOUBLE getBestGlobalAlpha(){return _bestGlobalAlpha;}
		Vdouble getBestL() {return _bestLvec;}
private:
	Vdouble _bestLocalAlphaVec;
	MDOUBLE _bestGlobalAlpha;
	Vdouble _bestLvec;
};

class bestBetaAndBBL {
public:
	explicit bestBetaAndBBL(tree& et, //find Best Beta and best BBL
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const MDOUBLE initBeta = 1.5,
		const MDOUBLE upperBoundOnBeta = 5.0,
		const MDOUBLE epsilonLoglikelihoodForBetaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
		const int maxBBLIterations=10,
		const int maxTotalIterations=5);
		MDOUBLE getBestBeta() {return _bestBeta;}
		MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestBeta;
	MDOUBLE _bestL;
};

class bestAlphaAndBetaAndBBL {
public:
	explicit bestAlphaAndBetaAndBBL(tree& et, //find Best Alpha and best BBL
		const sequenceContainer& sc,
		stochasticProcess& sp,
		const Vdouble * weights=NULL,
		const MDOUBLE initAlpha = 1.5,
		const MDOUBLE initBeta = 1.5,
		const MDOUBLE upperBoundOnAlpha = 5.0,
		const MDOUBLE upperBoundOnBeta = 5.0,
		const MDOUBLE epsilonLoglikelihoodForAlphaOptimization= 0.01,
		const MDOUBLE epsilonLoglikelihoodForBetaOptimization = 0.01,
		const MDOUBLE epsilonLoglikelihoodForBBL= 0.05,
		const int maxBBLIterations=10,
		const int maxTotalIterations=5);
		MDOUBLE getBestAlpha() {return _bestAlpha;}
		MDOUBLE getBestBeta() {return _bestBeta;}
		MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestBeta;
	MDOUBLE _bestL;
};


class C_evalAlpha{
public:
	C_evalAlpha(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_sp.categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize alpha");
		}
		(static_cast<generalGammaDistribution*>(_sp.distr()))->setAlpha(alpha);

		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with alpha = "<<alpha<<" logL = "<<res<<endl);
#ifdef VERBOS
		LOG(7,<<" while in brent: with alpha = "<<alpha<<" logL = "<<res<<endl);
#endif
		return -res;
	}
};

class C_evalLocalAlpha{
public:
	C_evalLocalAlpha(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_sp.categories() == 1) {
			errorMsg::reportError("one category when trying to optimize local alpha");
		}
		(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
		vector<sequenceContainer> tmpScVec;
		tmpScVec.push_back(_sc);
		vector<stochasticProcess> tmpSpVec;
		tmpSpVec.push_back(_sp);
		multipleStochasticProcess * tmpMsp = new multipleStochasticProcess();
		tmpMsp->setSpVec(tmpSpVec);	
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,tmpScVec,tmpMsp,_pProportionDist);
		MDOUBLE res = likeVec[0];
		delete(tmpMsp);
		LOG(5,<<" with local alpha = "<<alpha<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalGlobalAlpha{
public:
	C_evalGlobalAlpha(	const tree& et,
		vector<sequenceContainer>& sc,
		multipleStochasticProcess* msp,
		gammaDistribution* pProportionDist,
		const Vdouble * weights = NULL)
		: _et(et),_sc(sc),_weights(weights),_msp(msp),_pProportionDist(pProportionDist){};
private:
	const tree& _et;
	vector<sequenceContainer>& _sc;
	const Vdouble * _weights;
	multipleStochasticProcess* _msp;
	gammaDistribution* _pProportionDist;
public:
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_pProportionDist->categories() < 1) {
			errorMsg::reportError(" less than one category when trying to optimize global alpha");
		}
		_pProportionDist->setAlpha(alpha);
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,_sc,_msp,_pProportionDist);
		MDOUBLE res = sumVdouble(likeVec);
		LOG(5,<<" with global alpha = "<<alpha<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalBeta{
public:
	C_evalBeta(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE beta) {
		if (_sp.categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize beta");
		}
		(static_cast<generalGammaDistribution*>(_sp.distr()))->setBeta(beta);

		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with alpha = "<<alpha<<" logL = "<<res<<endl);
#ifdef VERBOS
		LOG(7,<<" while in brent: with beta = "<<beta<<" logL = "<<res<<endl);
#endif
		return -res;
	}
};

#endif

