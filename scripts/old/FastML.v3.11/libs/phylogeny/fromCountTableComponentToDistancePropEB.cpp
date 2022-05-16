// $Id: fromCountTableComponentToDistanceProp.cpp 962 2006-11-07 15:13:34Z privmane $

#include "fromCountTableComponentToDistancePropEB.h"
#include "likeDistPropEB.h"

fromCountTableComponentToDistancePropEB::fromCountTableComponentToDistancePropEB(
		const vector< vector<countTableComponentGamProportional> >& ctc,
		const int nodeID,
		multipleStochasticProcess *msp,
		const gammaDistribution* pProportionDist,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess ) : _msp(msp), _ctc(ctc), _nodeID(nodeID), _pProportionDist(pProportionDist){
	_distance =brLenIntialGuess;
	_toll = toll;
}

void fromCountTableComponentToDistancePropEB::computeDistance() {
	MDOUBLE maxPairwiseDistance = 10.0; // The default
	MDOUBLE minPairwiseDistance = 0.0000001; // The default
	likeDistPropEB likeDist1(_msp,_pProportionDist,_toll,maxPairwiseDistance,minPairwiseDistance);
	MDOUBLE initGuess = _distance;
	_distance = likeDist1.giveDistance(_ctc,_nodeID,_likeDistance,initGuess);
	assert(_distance>=0);
}
