#ifndef __indelCoderOptionsParams_OPTION
#define __indelCoderOptionsParams_OPTION

#include "definitions.h"
#include <string>
#include <fstream>

using namespace std;


class indelCoderOptions{
public:
	enum codingType {SIC, MCIC, MCIC2};

public:
	virtual ~indelCoderOptions();

	static void initOptions(const string& paramFileName);
	static void initDefault();
	static void readParameters(const string& paramFileName);
	static void getParamsFromFile(const string& paramFileName);
	static void getOutDirFromFile(const string& paramFileName);
	static void verifyConsistParams();

	static string getCodingType(codingType type);
	static codingType getCodingType(const string& str);


public:
//################### Basic parameters:
// input (general)
	static string _seqFile;		// essential - fasta file with presence(1)/absence(0) for each species over all gene families (positions)
	static string _indelOutputInfoFile;		// a file in which all the indel information is given (not just the 0/1 codes)
	static string _indelOutputFastaFile;	// a file in which ajust the 0/1 coding is given
	static string _nexusFileName; // a file in which the 0/1 coding is given in nexus format
	//static string _outDir;		// _outDir = "RESULTS", concatenated after current dir location 'pwd'
	static string _logFile;		// print-outs of the running progress including the estimated parameters optimization
	static int _logValue;		// verbosity level - ~4 - normal, >7 - load of info

	//static bool _isMCIC2;
	static codingType _codingType;		// SIC, MCIC, MCIC2

	static bool _isCheckForTriangleInequality;
	static bool _isOmitLeadingAndEndingGaps;	// ignore gaps that either start at 5' or end at 3'

private:

};
#endif
