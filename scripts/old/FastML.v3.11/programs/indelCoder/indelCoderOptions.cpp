/********************************************************************************************
indelCoderOptions - a class that contains all the parameters for the indelCoderProjest as static
use the 'Parameters' class to read info from txt file.
initDefault. (+Parameters::addParameter)
getParamsFromFile. ->with alterations of defults for consistancy
verifyConsistParams.
*********************************************************************************************/
#include "indelCoderOptions.h"
#include "errorMsg.h"
#include "someUtil.h"
#include "Parameters.h"
#include <iostream>
#include <cmath>

using namespace std;

// recognize all the static members defined at .h
string indelCoderOptions::_seqFile;
string indelCoderOptions::_logFile;
int indelCoderOptions::_logValue;
//string indelCoderOptions::_outDir;
string indelCoderOptions::_indelOutputInfoFile;	
string indelCoderOptions::_indelOutputFastaFile;
string indelCoderOptions::_nexusFileName;

indelCoderOptions::codingType indelCoderOptions::_codingType;
//bool indelCoderOptions::_isMCIC2;

bool indelCoderOptions::_isCheckForTriangleInequality;
bool indelCoderOptions::_isOmitLeadingAndEndingGaps;




/********************************************************************************************
*********************************************************************************************/
void indelCoderOptions::initOptions(const string& paramFileName)
{	
	//getOutDirFromFile(paramFileName);	// first set _outDir to be used next
	//createDir("", indelCoderOptions::_outDir);
	ifstream params(paramFileName.c_str());
	if(params.good())
		Parameters::readParameters(params);
	params.close();
	initDefault();
	getParamsFromFile(paramFileName);
	//verifyConsistParams();
}



/********************************************************************************************
*********************************************************************************************/
//void indelCoderOptions::getOutDirFromFile(const string& paramFileName)
//{
//	_outDir = "INDEL_CODER_RES";
//	Parameters::addParameter("_outDir", _outDir);
//	
//	_outDir = Parameters::getString("_outDir");
//}

/********************************************************************************************
initDefault
*********************************************************************************************/
void indelCoderOptions::initDefault()
{
// all the default values are stored in the gainLossOptions:: static members
//################### Basic parameters:
// input (general)
	_seqFile = "";				// essential - fasta file with presence(1)/absence(0) for each species over all gene families (positions)
	_indelOutputInfoFile= "";
	_indelOutputFastaFile="";
	_nexusFileName="";
// output
//_outDir = "RESULTS";		// concatenated after current dir location 'pwd'
	_logFile = "log.txt";		// print-outs of the running progress including the estimated parameters optimization
	_logValue = 4;				// verbosity level - ~4 - normal, >7 - load of info
	//_isMCIC2 = true;
	_codingType =SIC;
	_isCheckForTriangleInequality = false;
	_isOmitLeadingAndEndingGaps = true;	// The typical approach is to omit (SeqState)

	Parameters::addParameter("_seqFile", _seqFile);
	Parameters::addParameter("_logFile", _logFile);
	Parameters::addParameter("_indelOutputInfoFile", _indelOutputInfoFile);
	Parameters::addParameter("_indelOutputFastaFile", _indelOutputFastaFile);
	Parameters::addParameter("_nexusFileName", _nexusFileName);

	Parameters::addParameter("_logValue", _logValue);
	Parameters::addParameter("_codingType", getCodingType(_codingType));
	//Parameters::addParameter("_isMCIC2", (_isMCIC2 == true) ? 1 : 0);
	Parameters::addParameter("_isCheckForTriangleInequality", (_isCheckForTriangleInequality == true) ? 1 : 0);
	Parameters::addParameter("_isOmitLeadingAndEndingGaps", (_isOmitLeadingAndEndingGaps == true) ? 1 : 0);

}


/********************************************************************************************
getParamsFromFile
*********************************************************************************************/
void indelCoderOptions::readParameters(const string& paramFileName)
{
	ifstream params(paramFileName.c_str());
	if(params.good())
		Parameters::readParameters(params);		// only place where params are read, updateParameter(paramName, param.c_str()) used
	params.close();
}
/********************************************************************************************
getParamsFromFile
*********************************************************************************************/
void indelCoderOptions::getParamsFromFile(const string& paramFileName)
{	
	readParameters(paramFileName);
	_logFile = Parameters::getString("_logFile");
	_seqFile = Parameters::getString("_seqFile");
	_indelOutputFastaFile = Parameters::getString("_indelOutputFastaFile");
	_nexusFileName = Parameters::getString("_nexusFileName");
	_indelOutputInfoFile = Parameters::getString("_indelOutputInfoFile");
	if(_seqFile=="") errorMsg::reportError("_seqFile is needed");
	if(_indelOutputFastaFile=="") errorMsg::reportError("_indelOutputFastaFile is needed");
	if(_nexusFileName=="") errorMsg::reportError("_nexusFileName is needed");

	if(_indelOutputInfoFile=="") errorMsg::reportError("_indelOutputInfoFile is needed");
	//_isMCIC2 = (Parameters::getInt("_isMCIC2") == 1) ? true : false;
	_codingType = getCodingType(Parameters::getString("_codingType"));

	_isCheckForTriangleInequality = (Parameters::getInt("_isCheckForTriangleInequality") == 1) ? true : false;
	_isOmitLeadingAndEndingGaps = (Parameters::getInt("_isOmitLeadingAndEndingGaps") == 1) ? true : false;
	_logValue = Parameters::getInt("_logValue");
}

/********************************************************************************************
enum distributionType {SIC, MCIC, MCIC2};
*********************************************************************************************/
string indelCoderOptions::getCodingType(codingType type) 
{
	string res = "";
	switch (type)
	{
	case SIC:
		res = "SIC";
		break;
	case MCIC:
		res = "MCIC";
		break;
	case MCIC2:
		res = "MCIC2";
		break;
	default:
		errorMsg::reportError("unknown type in codingType - {SIC, MCIC, MCIC2}");
	}
	return res;
}
//////////////////////////////////////////////////////////////////////////
indelCoderOptions::codingType indelCoderOptions::getCodingType(const string& str) 
{
	if (str == "SIC")
		return SIC;
	if (str == "MCIC")
		return MCIC;
	if (str == "MCIC2")
		return MCIC2;
	else
		errorMsg::reportError("unknown type in codingType - {SIC, MCIC, MCIC2}");
	return SIC;
}



