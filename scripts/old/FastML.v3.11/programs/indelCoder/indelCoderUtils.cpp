
#include "indelCoderUtils.h"
#include "indelCoder.h"




void printICHelp(){
	cout <<"+-------------------------------------------+"<<endl;
	cout <<"*** The indelCoder project.					"<<endl;
	cout <<"use a parameter file with these options:    "<<endl;
	cout <<"+-------------------------------------------+"<<endl;
	cout <<"_seqFile									"<<endl;
	cout <<"|------------------------------------------|"<<endl;
	cout <<"_logFile									"<<endl;
	cout <<"_logValue									"<<endl;
	//cout <<"_outDir										"<<endl;
	cout <<"...(a partial list)							"<<endl;
	cout <<"+------------------------------------------+"<<endl;
}

void printICProgramInfo(){
	LOGnOUT(3,<<"+=================================================================+"<<endl);
	LOGnOUT(3,<<"+  The indelCoder project:									        "<<endl);
	LOGnOUT(3,<<"+  Transforming a multiple sequence alignment (MSA)                "<<endl);
	LOGnOUT(3,<<"+  of amino acids into 0/1 characters                              "<<endl);
	LOGnOUT(3,<<"+  Implementation of Indel Coding scheme SIC (Simmons et al. 2002) "<<endl);
	LOGnOUT(3,<<"+  "<<PROG_INFO<<"											        "<<endl);
	LOGnOUT(3,<<"+  Ofir Cohen - ofircohe@tau.ac.il	  							    "<<endl);
	LOGnOUT(3,<<"+  Tal Pupko  - talp@post.tau.ac.il 							    "<<endl);
	LOGnOUT(3,<<"+  Dorothee Huchon  - huchondp@post.tau.ac.il 			    	    "<<endl);
	LOGnOUT(3,<<"+=================================================================+"<<endl);
}
