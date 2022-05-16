#include "indelCoder.h"
#include "indelCoderOptions.h"
#include "indelCoderUtils.h"
#include "Parameters.h"



using namespace std;

int main(int argc, char **argv){

	//printICProgramInfo();
	//time_t t1,t2;
	//time(&t1);
	if (argc == 1) {printICHelp();// here the -h option will be printed
		return 0; 
	}
	string paramStr = argv[1];
	indelCoderOptions::initOptions(paramStr);

	myLog::setLog(indelCoderOptions::_logFile, indelCoderOptions::_logValue);

	//Parameters::dump(cout);

	indelCoder gl;
	gl.run();

	//time(&t2);
	//LOGnOUT(4,<<endl<<"TOTAL RUNNING TIME = "<<(t2-t1)/60.0<<" minutes"<<endl);
	return 0;
}

