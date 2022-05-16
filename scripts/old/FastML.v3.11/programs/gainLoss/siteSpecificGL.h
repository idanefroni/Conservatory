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


#ifndef ___SITE_SPECIFIC_GL__
#define ___SITE_SPECIFIC_GL__

#include "definitions.h"
#include "likelihoodComputation.h"
#include "seqContainerTreeMap.h"

// per all sites computation
void computeEB_EXP_siteSpecificGL(Vdouble & GainLossV,
								  Vdouble & stdV,
								  Vdouble & lowerBoundV,
								  Vdouble & upperBoundV,
								  VVdouble & posteriorsV,
								  const sequenceContainer& sc,
								  const vector<vector<stochasticProcess*> >& sp,
								  const tree& tr,
								  const distribution * gainDist,
								  const distribution * lossDist,
								  const distribution * distPrim,
								  const MDOUBLE alphaConf,
								  VVVdouble & postProbPerSpPerCatPerPos,	//2 fill (*postProbPerSpPerCatPerPos)[sp][pos]
								  unObservableData* unObservableData_p);

// per one site
void computeEB_EXP_siteSpecificGL(int pos,
									const sequenceContainer& sc,
									const vector<vector<stochasticProcess*> >& sp,
									//const computePijGam& cpg,
									const tree &tr,
									const distribution * gainDist,
									const distribution * lossDist,
									const distribution * distPrim,
									Vdouble & posteriorV,
									MDOUBLE& GainLossExpectation,
									MDOUBLE & stdForce,
									MDOUBLE & lowerConf,
									MDOUBLE & upperConf,
									const MDOUBLE alphaConf,
									VVVdouble & postProbPerSpPerCatPerPos,	//2 fill (*postProbPerSpPerCatPerPos)[sp][pos]
									unObservableData* unObservableData_p);








#endif
