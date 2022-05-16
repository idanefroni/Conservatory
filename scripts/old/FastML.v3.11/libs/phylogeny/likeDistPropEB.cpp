// $Id: likeDistProp.cpp 962 2006-11-07 15:13:34Z privmane $

#include "likeDistPropEB.h"
#include "numRec.h"

const MDOUBLE likeDistPropEB::giveDistance(	const vector< vector<countTableComponentGamProportional> >& ctc,const int nodeID,
								MDOUBLE& resL,const MDOUBLE initialGuess) const {
	const MDOUBLE ax = _minPairwiseDistance;
	const MDOUBLE bx = initialGuess;
	const MDOUBLE cx = _maxPairwiseDistance;
	const MDOUBLE tol = _toll;
	MDOUBLE dist=-1.0;
	resL = -dbrent(ax,bx,cx,
		  C_evallikeDistPropEB(ctc,_msp,_pProportionDist,nodeID),
		  C_evallikeDistPropEB_d(ctc,_msp,_pProportionDist,nodeID),
		  tol,&dist);
	return dist;
}

// the minus resL = -dbrent because C_evalDist return - value, because it is computing the min not the max...

