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


#ifndef ___RATE_4_TRIAD
#define ___RATE_4_TRIAD

#include "definitions.h"



class rate4Triad {
public:
	explicit rate4Triad(const stochasticProcess* sp, const Vdouble& exp01V, const Vdouble& exp10V);
	virtual ~rate4Triad(){};

	//void rate4Triad::computePosteriorExpectationOfChangePerTriad();






private:
	const stochasticProcess* _sp;
	//const Vdouble  &_rateV;
	const Vdouble  &_exp01V;
	const Vdouble  &_exp10V;	

};

#endif
