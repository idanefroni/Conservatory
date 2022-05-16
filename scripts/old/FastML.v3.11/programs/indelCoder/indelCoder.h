#ifndef ___INDELCODER_
#define ___INDELCODER_

#include "gaps.h"
#include "character.h"
#include "indelCoderOptions.h"

#include "amino.h"
#include "sequenceContainer.h"
#include "recognizeFormat.h"

#include "logFile.h"
#include "talRandom.h"
#include <ctime>
#include <vector>

using namespace std;
// The implementation of IndelCoding scheme MCIC (Muller)
// for "Large-scale parsimony analysis of metazoan indels in protein-coding genes"
// (Parsimony tree reconstruction using indel information - supporting the Ecdysozoa hypothesis - sponges evolutionary close to animals)
// coded as gaps only those gaps in the alignments that were shorter than 50 amino-acids and those which did not start at the N-terminus or end at the C-terminus of the alignment.
//
//
// Simple Indel Coding – SIC (Simmons and Ochoterena 2000) 
// each indel receives a separate 2-state character of presence/absence. 
// Any overlapping indels that exceed the boundaries of this indel are scored as missing data for that indel character.
//
// Modified Complex Indel Coding – MCIC (Muller 2006). 
// MCIC differs from SIC only in the treatment of overlapping indels. 
// uses multistate characters to code overlapping indels and assigns a distinct symmetrical step matrix to those gaps.

// Note:
// (*) implemented for Amino Acids seq. (Later, we can add parameter with Alphabet type)


class indelCoder {

public:
	explicit indelCoder(){};
	virtual ~indelCoder(){};

	void startSequenceContainer();
	void readSequenceIntoGaps();
	//void delimitationOfCharacters();
	void delimitationOfCharacters(indelCoderOptions::codingType type);
	void delimitationOfCharactersSIC();
	//void delimitationOfCharactersMCIC2();

	void determinationCharacterState();
	void determinationCharacterStateSIC();

	void determinationStepsMatrix();
	void printCharacters();
	void printNexus();
	void printFasta();
	void printGapsInfo();
	void printIndelSummary(); 
	void run();

private:
	sequenceContainer _sc;
	vector< vector<int> > _matrix;
	gaps _gaps;
	gaps _unknowns;
	vector<gaps> _gapsVperSc;
	vector<character*> _characters;
};






#endif

