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
#include "gainLoss.h"
#include "computePosteriorExpectationOfChange.h"
#include "gammaDistributionFixedCategories.h"
#include "mixtureDistribution.h"
#include "gainLossOptions.h"
#include "Parameters.h"

using namespace std;
int mainRunOptimize();

int main(int argc, char **argv){
	
	int pp = getpid();
	time_t t1,t2;
	time(&t1);
	if (argc == 1) {	
		printHelp();// here the -h option will be printed
		exit(0); 
	}
	long seed = static_cast<long>(t1) * pp;
	talRandom::setSeed(seed);	// set 1 for debug
	string paramStr = argv[1];
	gainLossOptions::initOptions(paramStr);

	myLog::setLog(gainLossOptions::_logFile, gainLossOptions::_logValue);
	LOG(4,<<"# Process_id= "<<getpid()<<endl);
	LOG(6,<<"# the seed = " <<seed<<endl);
	Parameters::dump(cout);
//enables other options...
	//string mainType = Parameters::getString("mainType");
	//if (mainType == "mainRunOptimize")	
		mainRunOptimize();	

	time(&t2);
	LOGnOUT(4,<<endl<<"TOTAL RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	return 0;
}

/********************************************************************************************
*********************************************************************************************/
int mainRunOptimize(){
	gainLoss gl;
	gl.run();
	return 0;
}
