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
#include "rate4Triad.h"
#include "errorMsg.h"
#include "logFile.h"
#include "gainLossOptions.h"

using namespace std;

/********************************************************************************************
rate4Triad
*********************************************************************************************/
rate4Triad::rate4Triad(const stochasticProcess* sp, const Vdouble& exp01V, const Vdouble& exp10V):
_sp(sp), _exp01V(exp01V), _exp10V(exp10V)
{
	if(!(_rateV.size()%3==0)){
		errorMsg::reportError("the length of the rates vector is not 'Triaded'");
	}
}

/********************************************************************************************
*********************************************************************************************/
//void rate4Triad::computePosteriorExpectationOfChangePerTriad(){
//	LOGnOUT(4,<<"Starting calculePosteriorExpectationOfChange for Triad..."<<endl);
//	
//	ofstream posteriorExpectationStreamTriad(gainLossOptions::_outFilePosteriorExpectationOfChangeTriad.c_str());
//	posteriorExpectationStreamTriad<<"POS"<<"\t"<<"000-001"<<"\t"<<"000-010"<<endl;
//
//
//	// printOut the final results
//	for (int pos = 0; pos <_sc.seqLen(); pos+=3){
//		posteriorExpectationStreamTriad<<pos+1<<"\t"<<expV01[pos]<<"\t"<<expV10[pos]<<endl;
//	}
//}
