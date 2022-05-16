// $Id: likeDistProp.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___LIKE_DIST_PROP_EB
#define ___LIKE_DIST_PROP_EB

#include "definitions.h"
#include "countTableComponent.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"
#include "logFile.h"
#include <cmath>

class likeDistPropEB {
private:
	multipleStochasticProcess * _msp;
	const gammaDistribution* _pProportionDist;
	const MDOUBLE _maxPairwiseDistance;
	const MDOUBLE _minPairwiseDistance;
	const MDOUBLE _toll;
public:
	const MDOUBLE giveDistance(	const vector< vector<countTableComponentGamProportional> >& ctc,const int nodeID,
								MDOUBLE& resL,const MDOUBLE initialGuess= 0.03) const;
	explicit likeDistPropEB(multipleStochasticProcess * msp,
							const gammaDistribution* pProportionDist,
							const MDOUBLE toll =0.0001,
							const MDOUBLE maxPairwiseDistance = 5.0,
							const MDOUBLE minPairwiseDistance = 0.0000001) : _msp(msp) ,_pProportionDist(pProportionDist), _maxPairwiseDistance(maxPairwiseDistance), _minPairwiseDistance(minPairwiseDistance),_toll(toll){
	}
	likeDistPropEB(const likeDistPropEB & other)
		: _msp(other._msp),_pProportionDist(other._pProportionDist),_maxPairwiseDistance(other._maxPairwiseDistance),_minPairwiseDistance(other._minPairwiseDistance),_toll(other._toll){}
	virtual likeDistPropEB* clone() const {return new likeDistPropEB(*this);}
};



class C_evallikeDistPropEB_d{ // derivative.
public:
  C_evallikeDistPropEB_d(const vector< vector<countTableComponentGamProportional> >& ctc,
				 multipleStochasticProcess* msp,const gammaDistribution* pProportionDist,const int nodeID) : _ctc(ctc), _msp(msp), _pProportionDist(pProportionDist), _nodeID(nodeID) {};
private:
	const vector< vector<countTableComponentGamProportional> >& _ctc;
	multipleStochasticProcess* _msp;
	const gammaDistribution* _pProportionDist;
	const int _nodeID;
public:
	MDOUBLE operator() (MDOUBLE dist) {
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumDL = 0.0;
		for (int gene=0; gene < _msp->getSPVecSize(); ++gene) {
			for (int alph1=0; alph1 <  _ctc[gene][_nodeID].alphabetSize(); ++alph1){
				for (int alph2=0; alph2 <  _ctc[gene][_nodeID].alphabetSize(); ++alph2){
					for(int globalRateCategor = 0;globalRateCategor < _pProportionDist->categories();++globalRateCategor){
						_msp->getSp(gene)->setGlobalRate(_pProportionDist->rates(globalRateCategor));
						MDOUBLE globalRate = _pProportionDist->rates(globalRateCategor);
						for (int localRateCategor = 0; localRateCategor < _msp->getSp(gene)->categories(); ++localRateCategor) {
							MDOUBLE localRate = _msp->getSp(gene)->rates(localRateCategor);
							MDOUBLE pij= _msp->getSp(gene)->Pij_t(alph1,alph2,dist*globalRate*localRate);
							if (pij<epsilonPIJ) {
								pij = epsilonPIJ;
							}
							MDOUBLE dpij = _msp->getSp(gene)->dPij_dt(alph1,alph2,dist*globalRate*localRate);
							
							//sumDL+= _ctc[gene][_nodeID].getCounts(alph1,alph2,globalRateCategor,localRateCategor)*dpij*_pProportionDist->ratesProb(globalRateCategor)*sp->ratesProb(localRateCategor)
							//			*globalRate*localRate/pij;
							sumDL+= _ctc[gene][_nodeID].getCounts(alph1,alph2,globalRateCategor,localRateCategor)*dpij*globalRate*localRate/pij;
						}
					}
				}
			}
		}
		LOG(12,<<"check bl="<<dist<<" gives sumDL "<<sumDL<<endl);
		return -sumDL;
	};
};



class C_evallikeDistPropEB{
private:
	const vector< vector<countTableComponentGamProportional> >& _ctc;
	multipleStochasticProcess* _msp;
	const gammaDistribution* _pProportionDist;
	const int _nodeID;
public:
	C_evallikeDistPropEB(const vector< vector<countTableComponentGamProportional> >& ctc,
					multipleStochasticProcess* msp,const gammaDistribution* pProportionDist,const int nodeID):_ctc(ctc), _msp(msp), _pProportionDist(pProportionDist), _nodeID(nodeID) {};

	MDOUBLE operator() (MDOUBLE dist) {
		const MDOUBLE epsilonPIJ = 1e-10;
		MDOUBLE sumL = 0.0;
		for (int gene=0; gene < _msp->getSPVecSize(); ++gene) {
			for (int alph1=0; alph1 < _ctc[gene][_nodeID].alphabetSize(); ++alph1){
				for (int alph2=0; alph2 <  _ctc[gene][_nodeID].alphabetSize(); ++alph2){
					for(int globalRateCategor = 0;globalRateCategor < _pProportionDist->categories();++globalRateCategor){
						_msp->getSp(gene)->setGlobalRate(_pProportionDist->rates(globalRateCategor));
						MDOUBLE globalRate = _pProportionDist->rates(globalRateCategor);
						for (int localRateCategor = 0; localRateCategor < _msp->getSp(gene)->categories(); ++localRateCategor) {
							MDOUBLE localRate = _msp->getSp(gene)->rates(localRateCategor);
							MDOUBLE pij= _msp->getSp(gene)->Pij_t(alph1,alph2,dist*globalRate*localRate);
							if (pij<epsilonPIJ) {
								pij = epsilonPIJ;
							}
							sumL += _ctc[gene][_nodeID].getCounts(alph1,alph2,globalRateCategor,localRateCategor)*(log(pij)-log(_msp->getSp(gene)->freq(alph2)));//*_pProportionDist->ratesProb(globalRateCategor)*sp->ratesProb(localRateCategor);
						}
					}
				}
			}
		}
		LOG(8,<<"check bl="<<dist<<" gives sumL "<<sumL<<endl);
		return -sumL;
	};
};

#endif

