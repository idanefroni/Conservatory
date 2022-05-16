// $Id: fromCountTableComponentToDistanceProp.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_PROP_EB
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE_PROP_EB

#include "definitions.h"
#include "countTableComponent.h"
#include "multipleStochasticProcess.h"
#include "gammaDistribution.h"

class fromCountTableComponentToDistancePropEB {

public:
	explicit fromCountTableComponentToDistancePropEB(
		const vector< vector<countTableComponentGamProportional> >& ctc,
		const int nodeID,
		multipleStochasticProcess* msp,
		const gammaDistribution* pProportionDist,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess = 0.029);// =startingGuessForTreeBrLen

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	multipleStochasticProcess * _msp;
	const vector< vector<countTableComponentGamProportional> >& _ctc;
	const gammaDistribution* _pProportionDist;
	const int _nodeID;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	int alphabetSize() {return (_ctc.empty()?0:_ctc[0][_nodeID].alphabetSize());}
};

#endif

