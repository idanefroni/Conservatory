// $Id: bestGtrModelparams.h 2008-28-04 15:13:34Z nimrod $

#ifndef ___BEST_GTRMODEL_PARAMS
#define ___BEST_GTRMODEL_PARAMS

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "generalGammaDistribution.h"
#include "tree.h"
#include "gtrModel.h"

typedef enum
  {
    Invalid = 0,
	a2c,
	a2g,
	a2t,
	c2g,
	c2t,
	g2t,
  }GTRParam;

#define maxBBLIt 10
#define epsilonLoglikeForBBL 0.01
#define inAlpha 1.5
#define epsilonLoglikeForAlphaOptimization 0.01
#define upperBoundForAlpha 5.0

class bestGtrModel {
public:
	explicit bestGtrModel(tree& et, // find best Gtr Model Params
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights=NULL,
										const int maxTotalIterations = 5,
										const MDOUBLE epsilonLikelihoodImprovment = 0.05,
										const MDOUBLE epsilonLoglikelihoodForGTRParam = 0.01,
										const MDOUBLE upperBoundGTRParam = 5.0,
										const bool optimizeTree = true,
                                        const bool optimizeAlpha = true);
	MDOUBLE getBesta2c() {return _best_a2c;}
	MDOUBLE getBesta2g() {return _best_a2g;}
	MDOUBLE getBesta2t() {return _best_a2t;}
	MDOUBLE getBestc2g() {return _best_c2g;}
	MDOUBLE getBestc2t() {return _best_c2t;}
	MDOUBLE getBestg2t() {return _best_g2t;}
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _best_a2c;
	MDOUBLE _best_a2g;
	MDOUBLE _best_a2t;
	MDOUBLE _best_c2g;
	MDOUBLE _best_c2t;
	MDOUBLE _best_g2t;
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestGtrModelProportional {
public:
	explicit bestGtrModelProportional(tree& et, // find best Gtr Model Params under a proportional model
										vector<sequenceContainer>& sc,
										multipleStochasticProcess* msp,
 										gammaDistribution* pProportionDist,
										Vdouble initLocalAlphas,
										Vdouble initLocala2cs,
										Vdouble initLocala2gs,
										Vdouble initLocala2ts,
										Vdouble initLocalc2gs,
										Vdouble initLocalc2ts,
										Vdouble initLocalg2ts,
										const MDOUBLE upperBoundOnLocalAlpha,
										const MDOUBLE initGlobalAlpha,
										const MDOUBLE upperBoundOnGlobalAlpha,
										const MDOUBLE upperBoundGTRParam,
										const int maxTotalIterations,
										const int maxBBLIterations,
										const bool optimizeSelectedBranches=false,
										const bool optimizeTree = true,
                                        const string branchLengthOptimizationMethod="bblLS",
										const bool optimizeLocalParams = true,
										const bool optimizeGlobalAlpha = true,
										const Vdouble * weights=NULL,
										const MDOUBLE epsilonLikelihoodImprovment = 0.05,
										const MDOUBLE epsilonLoglikelihoodForGTRParam = 0.01,
										const MDOUBLE epsilonLoglikelihoodForLocalAlphaOptimization= 0.01,
										const MDOUBLE epsilonLoglikelihoodForGlobalAlphaOptimization= 0.01,
										const MDOUBLE epsilonLoglikelihoodForBBL= 0.01);
	MDOUBLE getBesta2c(int spIndex) {return _best_a2cVec[spIndex];}
	MDOUBLE getBesta2g(int spIndex) {return _best_a2gVec[spIndex];}
	MDOUBLE getBesta2t(int spIndex) {return _best_a2tVec[spIndex];}
	MDOUBLE getBestc2g(int spIndex) {return _best_c2gVec[spIndex];}
	MDOUBLE getBestc2t(int spIndex) {return _best_c2tVec[spIndex];}
	MDOUBLE getBestg2t(int spIndex) {return _best_g2tVec[spIndex];}
	MDOUBLE getBestLocalAlpha(int spIndex) {return _bestLocalAlphaVec[spIndex];}
	MDOUBLE getBestGlobalAlpha() {return _bestGlobalAlpha;}
	Vdouble getBestL() {return _bestLvec;}
private:
	Vdouble _best_a2cVec;
	Vdouble _best_a2gVec;
	Vdouble _best_a2tVec;
	Vdouble _best_c2gVec;
	Vdouble _best_c2tVec;
	Vdouble _best_g2tVec;
	Vdouble _bestLocalAlphaVec;
	MDOUBLE _bestGlobalAlpha;
	Vdouble _bestLvec;
};

class C_evalGTRParam{
public:
  C_evalGTRParam(	const GTRParam param,
					const tree& et,
					const sequenceContainer& sc,
					stochasticProcess& sp,
					const Vdouble * weights = NULL)
		:_param(param), _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const GTRParam _param;
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE paramVal) {
		switch (_param){
			case a2c:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2c(paramVal);
				break;
			case a2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2g(paramVal);
				break;	
			case a2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2t(paramVal);
				break;
			case c2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2g(paramVal);
				break;
			case c2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2t(paramVal);
				break;
			case g2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_g2t(paramVal);
				break;
			default:
				errorMsg::reportError("Missing GTR parameter in C_evalGTRParam::operator ()");
				break;
		}
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		LOG(5,<<" with " + int2string(_param) + " = "<<paramVal<<" logL = "<<res<<endl);
		return -res;
	}
};

class C_evalGTRParamProportional{
public:
  C_evalGTRParamProportional(	const GTRParam param,
					const tree& et,
					const sequenceContainer& sc,
                    stochasticProcess& sp,
					const gammaDistribution* pProportionDist,
					const Vdouble * weights = NULL)
		:_param(param), _et(et),_sc(sc),_sp(sp),_pProportionDist(pProportionDist),_weights(weights){};
private:
	const GTRParam _param;
	const tree& _et;
	const sequenceContainer& _sc;
	const gammaDistribution* _pProportionDist;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE paramVal) {
		switch (_param){
			case a2c:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2c(paramVal);
				break;
			case a2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2g(paramVal);
				break;	
			case a2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_a2t(paramVal);
				break;
			case c2g:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2g(paramVal);
				break;
			case c2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_c2t(paramVal);
				break;
			case g2t:
				(static_cast<gtrModel*>(_sp.getPijAccelerator()->getReplacementModel()))->set_g2t(paramVal);
				break;
			default:
				errorMsg::reportError("Missing GTR parameter in C_evalGTRParamProportional::operator ()");
				break;
		}
		vector<sequenceContainer> tmpScVec;
		tmpScVec.push_back(_sc);
		vector<stochasticProcess> tmpSpVec;
		tmpSpVec.push_back(_sp);
		multipleStochasticProcess * tmpMsp = new multipleStochasticProcess();
		tmpMsp->setSpVec(tmpSpVec);	
		Vdouble likeVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(_et,tmpScVec,tmpMsp,_pProportionDist);
		MDOUBLE res = likeVec[0];
		delete(tmpMsp);
		LOG(5,<<" with " + int2string(_param) + " = "<<paramVal<<" logL = "<<res<<endl);
		return -res;
	}
};

#endif


