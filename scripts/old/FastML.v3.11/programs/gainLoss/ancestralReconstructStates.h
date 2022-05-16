/*
Copyright (C) 2011 Tal Pupko  TalP@tauex.tau.ac.il.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef ___ANCESTRAL_RECONSTRUCT_STATES_
#define ___ANCESTRAL_RECONSTRUCT_STATES_

#include "computeDownAlg.h"
#include "definitions.h"
#include "gainLossAlphabet.h"
#include "matrixUtils.h"
#include "seqContainerTreeMap.h"
#include "sequence.h"
#include "someUtil.h"
#include "treeIt.h"
#include "trivialAccelerator.h"
#include "logFile.h"


class ancestralReconstructStates {

public:
	explicit ancestralReconstructStates(const tree &tr, const sequenceContainer &sc, stochasticProcess *sp);
	virtual ~ancestralReconstructStates(){};

	void traverseUpML(VVVdouble &upL, VVVint &backtrack); // input as empty vector to be filled
	void traverseUpML(VVdouble &upL, VVint &backtrack, int pos);

	Vdouble traverseDownML(VVVdouble &upL, VVVint &backtrack,VVVint &transitionTypeCount); // input as already filled vector
	MDOUBLE traverseDownML(VVdouble &upL, VVint &backtrack,VVint &transitionTypeCount, int pos);
	VVint getStates() {return _statesV;}
	
	void computeAncestralPosterior(const VVVVdouble& jointPost); // posterior (marginal) reconstruction
	VVVdouble getAncestralProbs() {return _ancestralProbs;}


private:
	void initializeStatesVector(int pos);


	const tree &_tr;
	const sequenceContainer &_sc;
	stochasticProcess *_sp;
	VVint _statesV;
	VVVdouble _ancestralProbs; // VVVdouble[pos][node][state] _ancestralProbs
};

#endif
