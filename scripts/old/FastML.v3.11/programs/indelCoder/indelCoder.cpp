#include "indelCoder.h"
#include "indelCoderUtils.h"


using namespace std;


/********************************************************************************************
run
*********************************************************************************************/
void indelCoder::run(){
	startSequenceContainer();	// note: only Amino seq is implemented
	readSequenceIntoGaps();		// Find and sort all gaps in MSA
	printGapsInfo();
	switch (indelCoderOptions::_codingType)
	{
		case (indelCoderOptions::SIC):
			delimitationOfCharactersSIC();
			break;
		case (indelCoderOptions::MCIC):
			LOGnOUT(2,<<endl<< "WARNING: The MCIC implementation is incomplete.\n Please re-run using SIC coding\n"<<endl);
			return;
			delimitationOfCharacters(indelCoderOptions::_codingType);
			break;
		case (indelCoderOptions::MCIC2):
			LOGnOUT(2,<<endl<< "WARNING: The MCIC2 implementation is incomplete.\n Please re-run using SIC coding\n"<<endl);
			return;
			delimitationOfCharacters(indelCoderOptions::_codingType);
			break;
		default:
			errorMsg::reportError("unknown type in codingType - {SIC, MCIC, MCIC2}");
	}	
	
	//if(indelCoderOptions::_isMCIC2)
	//	delimitationOfCharactersMCIC2();
	//else
	//	delimitationOfCharacters();
	
	resizeMatrix(_matrix,_sc.numberOfSeqs(),_characters.size());
	if(indelCoderOptions::_codingType==indelCoderOptions::SIC)
		determinationCharacterStateSIC();
	else{
		determinationCharacterState();
		determinationStepsMatrix();		
	}
	//printCharacters(); //DEBUG	
	printFasta();
	printNexus();
	printIndelSummary(); // Tal's own format in which all indel information is provided.
}

/********************************************************************************************
printCharacters
*********************************************************************************************/
void indelCoder::printCharacters(){
	for(int i = 0; i < _characters.size(); ++i)
	{
		cout<<"Character  "<<_characters[i]->getCoord5()<<" "<<_characters[i]->getCoord3()<<" "<<_characters[i]->getNumOfGaps()<<" "<<_characters[i]->getNumOfStates()<<endl;
		//_characters[i]->printScStates(); //DEBUG
		_characters[i]->printStates();
		_characters[i]->printStepsMatrix();
	}
}

/********************************************************************************************
startSequenceContainer
*********************************************************************************************/
void indelCoder::startSequenceContainer(){
	amino alph;	// note: we can add parameter with Alphabet type
	ifstream in(indelCoderOptions::_seqFile.c_str());
	_sc = recognizeFormat::read(in,&alph);
	_gapsVperSc.resize(_sc.numberOfSeqs());
	LOGnOUT(4,<<"Seq "<<indelCoderOptions::_seqFile.c_str()<<endl);	
	LOGnOUT(4,<<"numberOfSeqs="<<_sc.numberOfSeqs()<<endl);	
	LOGnOUT(4,<<"seqLen="<<_sc.seqLen()<<endl);
}

/********************************************************************************************
readSequenceIntoGaps
// Sort the vector containing all indels by  I =(i1,i2), K =(k1,k2), I<K iff i1<k1 or i1=k1 and i2<k2
*********************************************************************************************/
void indelCoder::readSequenceIntoGaps(){
	LOGnOUT(4,<<endl<< "Step (1) readSequenceIntoGaps..."<<endl);	
	LOGnOUT(5,<< " All MSA gaps are sorted by coordinates"<<endl);	
	LOGnOUT(5,<< " Sort by: I =(i1,i2), K =(k1,k2), I<K iff i1<k1 or i1=k1 and i2<k2"<<endl);
	int gapSign = -1;
	int UnknownSign = _sc.getAlphabet()->unknown(); // Note that within amino class, 'X' is also coded as unknown
	int coord5=0;
	int coord3=0;
	int seqID=0;
	int coord5abs=0; //coord5MinusNumOfGapPositionsFromGenomeStart
	for(int s=0; s<_sc.numberOfSeqs(); ++s){
		cout<<_sc[s].id()<<" "<<_sc[s].name()<<"\n";
		int numOfGapPositionsFromGenomeStart = 0;
		int seqLength =_sc.seqLen();
		for(int pos=0; pos<seqLength; ++pos){
			if(_sc[s][pos] == gapSign){				
				coord5 = pos;
				coord5abs = coord5-numOfGapPositionsFromGenomeStart;
				++numOfGapPositionsFromGenomeStart;
				while(pos<(seqLength-1) && _sc[s][pos+1] == gapSign){
					++pos;
					++numOfGapPositionsFromGenomeStart;
				}					
				coord3 = pos;
				seqID = _sc[s].id();
				//cout<<"new gap found "<<seqID<<" "<<coord5<<" "<<coord3<<endl;
				if(indelCoderOptions::_isOmitLeadingAndEndingGaps && (coord5abs==0 || coord3==seqLength-1)){
					LOGnOUT(4,<< "Skip Leading/Ending Gap. seq="<< s<<" coord5="<<coord5abs<<" coord3="<<coord3<<endl);
					_unknowns.insertNewGap(coord5, coord3,seqID, coord5abs);
				}else{
					_gaps.insertNewGap(coord5, coord3,seqID, coord5abs);
					_gapsVperSc[seqID].insertNewGap(coord5, coord3,seqID, coord5abs);	// used additionally were gaps are pre-sorted by seq
				}				
			}
			if(_sc[s][pos] == UnknownSign){
				coord5 = pos;
				coord5abs = coord5-numOfGapPositionsFromGenomeStart;
				while(pos<(seqLength-1) && _sc[s][pos+1] == UnknownSign){
					++pos;
				}					
				coord3 = pos;
				seqID = _sc[s].id();
				_unknowns.insertNewGap(coord5, coord3,seqID, coord5abs);
			}
		}
	}
	LOGnOUT(4,<<endl<< "There are "<<_gaps.numOfGaps()<<" gaps"<<endl);
	//_gaps.printGaps();	//DEBUG
}

/********************************************************************************************
delimitationOfCharacters
// 1) delimitation of the characters. 
// Each character is a region of the alignment that is fully represented by one indel, 
// and this indel is the longest one in these coordinates.

// Start with the first gap(indel) in the sorted vector to define the first character, following characters 
// character_1 is defined by gap_1
// For i in gap_1:gap_N (N gaps in the sorted vector)
//  while gap_i(3') < character_j(3')
//   gap_i is within character_j
//  else j++
*********************************************************************************************/
void indelCoder::delimitationOfCharacters(indelCoderOptions::codingType type){
	LOGnOUT(4,<<endl<< "Step (2) delimitationOfCharacters... Complex Coding "<<indelCoderOptions::getCodingType(type)<<endl);	
	LOGnOUT(5,<< " Finding the required positions(=characters) in the coded sequence"<<endl);	
	LOGnOUT(5,<< " The number of characters <= the number of found gaps (each character may include several gaps)"<<endl);

	if(type == indelCoderOptions::MCIC)
		LOGnOUT(5,<< " gap_i is extending previous character if it's start is the same but ends further"<<endl);
	if(type == indelCoderOptions::MCIC2)
		LOGnOUT(5,<< " gap_i is extending previous character if it's start is included in the previous character but ends further"<<endl);

	int i=0;
	character* character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
	_characters.push_back(character_p);
	while( i<_gaps.numOfGaps())
	{
		// gap_i is included in previous character if it's start(5) & end(3) coord are within the start & end coord of the character
		while(  _gaps[i]->getCoord5()>=character_p->getCoord5() && _gaps[i]->getCoord3()<=character_p->getCoord3()){
			character_p->addGap(_gaps[i],i+1);
			i++;
			if(i>=_gaps.numOfGaps())
				break;
		}
		if(i>=_gaps.numOfGaps())
			break;

		bool condition;
		if(type == indelCoderOptions::MCIC)
			condition = _gaps[i]->getCoord5()==character_p->getCoord5(); // gap_i is extending previous character if it's start is the same but ends(3) further
		if(type == indelCoderOptions::MCIC2)
			condition = _gaps[i]->getCoord5()<=character_p->getCoord3(); // gap_i is extending previous character if it's start is included in the previous character but ends(3) further		
		while(condition && _gaps[i]->getCoord3()>character_p->getCoord3() ){
			character_p->setCoord3(_gaps[i]->getCoord3());
			character_p->addGap(_gaps[i],i+1);
			i++;
			if(i>=_gaps.numOfGaps())
				break;
		}
		
		// new character is required for this gap
		if(i<_gaps.numOfGaps() &&  _gaps[i]->getCoord5() > character_p->getCoord5() && _gaps[i]->getCoord3() > character_p->getCoord3()){
			character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
			_characters.push_back(character_p);
			character_p->addGap(_gaps[i],i+1);
			i++;
			//cout<<" Character  "<<i<<" "<< character_p->getCoord5()<<" "<< character_p->getCoord3()<<" " <<endl;
			//break;
		}
	}
	LOGnOUT(4,<<endl<< "There were "<<_characters.size()<<" characters"<<endl);
}

/********************************************************************************************
delimitationOfCharactersMCIC2
*********************************************************************************************/
//void indelCoder::delimitationOfCharactersMCIC2(){
//	LOGnOUT(4,<<endl<< "Step (2) delimitationOfCharacters... Complex Coding (MCIC2)"<<endl);	
//	LOGnOUT(5,<< " Finding the required positions(=characters) in the coded sequence"<<endl);	
//	LOGnOUT(5,<< " The number of characters <= the number of found gaps (each character may include several gaps)"<<endl);
//
//	int i=0;
//	character* character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
//	_characters.push_back(character_p);
//	while( i<_gaps.numOfGaps())
//	{
//		//coord5_c = _gaps[i]->getCoord5();
//		//coord3_c = _gaps[i]->getCoord3();		
//		while(  _gaps[i]->getCoord5()>=character_p->getCoord5() && _gaps[i]->getCoord3()<=character_p->getCoord3()){
//			character_p->addGap(_gaps[i],i+1);
//			i++;
//			if(i>=_gaps.numOfGaps())
//				break;
//		}
//		// gap_i is extending previous character it's start is included in the previous character
//		while(i<_gaps.numOfGaps() 
//				&& _gaps[i]->getCoord5()<=character_p->getCoord3() && _gaps[i]->getCoord3()>character_p->getCoord3() ){
//			character_p->setCoord3(_gaps[i]->getCoord3());
//			character_p->addGap(_gaps[i],i+1);
//			i++;
//			if(i>=_gaps.numOfGaps())
//				break;
//		}
//		// new character is required for this gap
//		if(i<_gaps.numOfGaps() &&  _gaps[i]->getCoord5() > character_p->getCoord5() && _gaps[i]->getCoord3() > character_p->getCoord3()){
//			character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
//			_characters.push_back(character_p);
//			character_p->addGap(_gaps[i],i+1);
//			i++;
//			//cout<<" Character  "<<i<<" "<< character_p->getCoord5()<<" "<< character_p->getCoord3()<<" " <<endl;
//			//break;
//		}
//	}
//	LOGnOUT(4,<<endl<< "There were "<<_characters.size()<<" characters"<<endl);
//}


/********************************************************************************************
delimitationOfCharactersSIC
*********************************************************************************************/
void indelCoder::delimitationOfCharactersSIC(){
	LOGnOUT(4,<<endl<< "Step (2) delimitationOfCharacters... Simple Coding (SIC)"<<endl);	
	LOGnOUT(5,<< " Finding the required positions(=characters) in the coded sequence"<<endl);	
	LOGnOUT(5,<< " The number of characters <= the number of found gaps (each character may include several gaps)"<<endl);

	int i=0;
	character* character_p=NULL;
	if (_gaps.numOfGaps()>0) {
        character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
		_characters.push_back(character_p);
	}
	while( i<_gaps.numOfGaps())
	{
		while(  _gaps[i]->getCoord5()==character_p->getCoord5() && _gaps[i]->getCoord3()==character_p->getCoord3()){
			character_p->addGap(_gaps[i],i+1);
			i++;
			if(i>=_gaps.numOfGaps())
				break;
		}
		// new character is required for this gap
		if(i<_gaps.numOfGaps() ){
			character_p = new character(_gaps[i]->getCoord5(), _gaps[i]->getCoord3(),_sc.numberOfSeqs());
			_characters.push_back(character_p);
			character_p->addGap(_gaps[i],i+1);
			i++;
			//cout<<" Character  "<<i<<" "<< character_p->getCoord5()<<" "<< character_p->getCoord3()<<" " <<endl;
			//break;
		}
	}
	LOGnOUT(4,<<endl<< "There were "<<_characters.size()<<" characters"<<endl);
}


/********************************************************************************************
determinationCharacterState
// 2) determination of the character state of each character. 
// Each sequence presenting a different indel pattern at the corresponding character region is coded as a different state.
// state_0 is defined by no-gaps in this region 
// Foreach j in character_1:character_M (M where found in previous step)
//  Foreach s in seq1:seqS (S taxons in the MSA)
//   if s gaps coordinates are equal to one of the gaps coordinates of previous states
//    next;
//   else 
//    st++, where state_st is defined by X set of 5'-3' coordinates of the X gaps within s  
*********************************************************************************************/
void indelCoder::determinationCharacterState(){
	LOGnOUT(4,<<endl<< "Step (3) determinationCharacterState... "<<endl);	

	// loop over characters
	for(int c = 0; c < _characters.size(); ++c)
	{
		int coord_5p = _characters[c]->getCoord5();
		int coord_3p = _characters[c]->getCoord3();
		_characters[c]->addZeroState();	// the default state - ones in length of the character
		//cout<<" Char "<<" "<<coord_5p<<" "<<  coord_3p<<endl;
		_characters[c]->resizeSc_states();	// the _sc_states matrix (#species X length) is init with ones
		
		// loop over taxa
		for(int s = 0; s < _sc.numberOfSeqs(); ++s){
			int seqID = _sc[s].id();
			if(seqID != s)
				cout<<"error: seqID not eq s";
			//cout<<"SeqID vs. s "<<seqID<<" "<<s<<endl; // DEBUG			
			
			// loop over gaps - ToDo - highly wasteful! - fix
			for(int g = 0; g <_gapsVperSc[seqID].numOfGaps(); ++g){
				int coord_5_gap = _gapsVperSc[seqID][g]->getCoord5();
				int coord_3_gap = _gapsVperSc[seqID][g]->getCoord3();
				if(coord_5_gap>=coord_5p && coord_3_gap<= coord_3p){
					// if the gap of the species is included in character - update _sc_state with zeros to designate this gap
					_characters[c]->setGapsInSc_states( seqID, coord_5_gap, coord_3_gap); 
				}
			}

			bool isNewState = true;
			if(_characters[c]->getScStates()[s]==_characters[c]->getStates()[0]  ){
				isNewState = false;
				_matrix[s][c] = 0;	// no gaps for species s in character c
			}
			else{
				for(int sq = s-1; sq>=0; --sq){
					if(_characters[c]->getScStates()[s] == _characters[c]->getScStates()[sq] ){
						isNewState = false;				// this state was already found, no need for new
						_matrix[s][c] = _matrix[sq][c];	// gap in species s, same state as previously found state in species sq
					}
				}
			}
			if(isNewState){	// state was not found in previous species, need new
				_characters[c]->addState(_characters[c]->getScStates()[s]);
				_matrix[s][c] = _characters[c]->getNumOfStates()-1;	// gap in species s, new state type
			}
		}
	}
}

/********************************************************************************************
determinationCharacterState
// 2) determination of the character state of each character. 
// Each sequence presenting a different indel pattern at the corresponding character region is coded as a different state.
// state_0 is defined by no-gaps in this region 
// state_1 for species with this (exact) gap
// state_? for species with gap overlapping this one
*********************************************************************************************/
void indelCoder::determinationCharacterStateSIC(){
	LOGnOUT(4,<<endl<< "Step (3) determinationCharacterState... "<<endl);	

	resizeMatrix(_matrix,_sc.numberOfSeqs(),_characters.size());	// all zeroes
	// loop over characters
	for(int c = 0; c < _characters.size(); ++c)
	{
		int coord_5_char = _characters[c]->getCoord5();
		int coord_3_char = _characters[c]->getCoord3();

		// loop over all gaps
		for(int g = 0; g <_gaps.numOfGaps(); ++g){
			int coord_5_gap = _gaps[g]->getCoord5();
			int coord_3_gap = _gaps[g]->getCoord3();
			int s = _gaps[g]->getSeqID();	// ????
			string nameG = _sc.name(_gaps[g]->getSeqID());
			
			if(coord_5_gap==coord_5_char && coord_3_gap==coord_3_char){
				_matrix[s][c] = 1;				
			}				
			else if(	//(coord_5_gap>=coord_5_char && coord_3_gap<= coord_3_char) ||		//  gap (in genome 'g') is within the char (indel=c)
						(coord_5_gap<=coord_5_char && coord_3_gap>= coord_3_char ) //||      // char (indel=c) is within the gap (in genome 'g')
						//(coord_5_gap<=coord_5_char && coord_3_gap>=coord_5_char ) ||      // 5' of char is within gap (partial overlap)  
						//(coord_5_gap<=coord_3_char && coord_3_gap>=coord_3_char )        // 3' of char is within gap (partial overlap)  
				)
				_matrix[s][c] = 2;	// same as '?', need to find&replace			
		}
		for(int g = 0; g <_unknowns.numOfGaps(); ++g){
			int coord_5_gap = _unknowns[g]->getCoord5();
			int coord_3_gap = _unknowns[g]->getCoord3();
			int s = _unknowns[g]->getSeqID();	// ????
			if(coord_5_gap<=coord_5_char && coord_3_gap>= coord_3_char )
				_matrix[s][c] = 2;	// same as '?', need to find&replace

			if(_matrix[s][c] == 1 && (coord_5_char==(coord_3_gap+1) || coord_3_char==(coord_5_gap-1)) ) // the indel is flaked by unKnown, thus it is ?
				_matrix[s][c] = 2;	// same as '?', need to find&replace	
		}

	}

	
}



/********************************************************************************************
determinationStepsMatrix
// 3) determination of the number of steps between every 2 character states. 
// Each pair of sequences is compared separately for the corresponding character area and the minimum number of steps between every 2 character states is then determined

// cost_c_x_y is initiated with lenght of character c
// Foreach c in character_1:character_M
//  Foreach st_x and st_y in state_0:state_c_ST (There are ST states in character c) (go over all state combinations)
//   Do A to E steps for the pair st_x and st_y:
//   A) translate into 01 to X set of 5'-3' coordinats of the X gaps within st_x and st_y
//   B) ignore 0-0 (cost_c_x_y =- #0-0 colomns)
//   C) merge adjacent 0-1 and 1-0 ((cost_c_x_y =- #adjacent 0-1 and 1-0 colomns)
//   D) ignore 1-1 (cost_c_x_y =- #1-1 colomns)
*********************************************************************************************/
void indelCoder::determinationStepsMatrix(){
	LOGnOUT(4,<<endl<< "determinationStepsMatrix..."<<endl);
	for(int c = 0; c < _characters.size(); ++c)
	{
        _characters[c]->determinationStepsMatrix();
	}
}
/********************************************************************************************
//  print to out file the required data
*********************************************************************************************/
void indelCoder::printGapsInfo(){
	string fileGapsString = "gapsInfo.txt";
	ofstream fileGapsStream(fileGapsString.c_str());

	fileGapsStream<<"# Start coordinate are with the genome as reference (not MSA)."<<endl;
	fileGapsStream<<"# Count starts from zero (first position = 1)."<<endl;
	fileGapsStream<<"seqID"<<"\t"<<"seqName"<<"\t"<<"start"<<"\t"<<"length"<<endl;
	int gapNum = 1;
	for(int s = 0; s < _sc.numberOfSeqs(); ++s){
		int seqID = _sc[s].id();
		string seqName = _sc[s].name();
		for(int g = 0; g <_gapsVperSc[seqID].numOfGaps(); ++g){
			fileGapsStream<<seqID<<"\t"<<seqName<<"\t"<<_gapsVperSc[seqID][g]->getCoord5Abs()+1<<"\t"<<_gapsVperSc[seqID][g]->getLength()<<endl;
		}
	}
}


/********************************************************************************************
// print to out file the required data as fasta file
*********************************************************************************************/
void indelCoder::printFasta(){
	//string fileString = indelCoderOptions::_outDir + "//" + "outFileCodedSeq.fa";
	//ofstream fileStream(fileString.c_str());
	ofstream fileStream(indelCoderOptions::_indelOutputFastaFile.c_str());
	bool isSIC  = indelCoderOptions::_codingType == indelCoderOptions::SIC;
	for(int s=0; s<_sc.numberOfSeqs();++s){
		fileStream<<">"<<_sc.name(s)<<"\n";
		for(int c=0; c<_matrix[0].size(); ++c ){
			if(isSIC && _matrix[s][c]==2)
				fileStream<<'?';
			else
				fileStream<<_matrix[s][c];
		}
		fileStream<<endl;	// prev- 2 endl
	}
	fileStream.close();
}


/********************************************************************************************
// 4) print to out file the required data
// 4.1) MATRIX of S species over M characters
// 4.2) foreach character of more than 2 states print the transition costs stepmatrix
*********************************************************************************************/
void indelCoder::printNexus()  {
	bool isSIC  = indelCoderOptions::_codingType == indelCoderOptions::SIC;
	string fileNexusString = indelCoderOptions::_nexusFileName;
	ofstream fileNexusStream(fileNexusString.c_str());
	fileNexusStream<<"#NEXUS"<<endl<<endl;
	fileNexusStream<<"[!matrix with indels coded according to "<<indelCoderOptions::getCodingType(indelCoderOptions::_codingType)<<" coding]"<<endl;
		//if(indelCoderOptions::_isMCIC2){fileNexusStream<<"2";}
		//fileNexusStream<<"]"<<endl;
	fileNexusStream<<"[! "<< PROG_INFO <<" ]"<<endl<<endl;

	fileNexusStream<<"BEGIN CHARACTERS;"<<endl;
	fileNexusStream<<"DIMENSIONS newtaxa ntax="<<_sc.numberOfSeqs()<<" NCHAR="<<_sc.seqLen()+_characters.size()<<";"<<endl<<endl;

	fileNexusStream<<"FORMAT "<<endl;
	fileNexusStream<<"	DATATYPE = standard"<<endl; 
	fileNexusStream<<"	GAP = - "<<endl;
	fileNexusStream<<"	MISSING = ?"<<endl;
	fileNexusStream<<"	SYMBOLS="<<'"'<<"0123456789A#C$EFG.IJ&L%>OPQ/'TU:*X<Z"<<'"'<<endl;
	fileNexusStream<<"	EQUATE="<<'"'<<"R={AG} Y={CT} M={AC} K={GT} S={CG} W={AT} H={ACT} B={CGT} V={ACG} D={AGT} N={ACGT} r={AG} y={CT} m={AC} k={GT} s={CG} w={AT} h={ACT} b={CGT} v={ACG} d={AGT} n={ACGT}"<<'"'<<endl;
	fileNexusStream<<"INTERLEAVE;"<<endl;
	
	
	//string fileGapString = indelCoderOptions::_outDir + "//" + "outGapInfoCHARSTATELABELS.txt";
	//ofstream gapStream(fileGapString.c_str());	
	//gapStream<<"# each character is listing all included gaps.\n";
	//gapStream<<"# date for each gap: seqID(num start from 0), coord5Abs(not by MSA), length.\n";
	
	fileNexusStream<<"CHARSTATELABELS"<<endl;
	int seqLeng = _sc.seqLen();
	int c=0;
	for(c=0; c<_characters.size();++c){
		fileNexusStream<<"\t"<<c+1+seqLeng<<" "<<"ind_pos_"<<_characters[c]->getCoord5()+1<<"_to_"<<_characters[c]->getCoord3()+1<<" ";
		fileNexusStream<<"/absent ";	
//		gapStream<<c<<" character\t"<<_characters[c]->getCoord5()+1<<" to "<<_characters[c]->getCoord3()+1<<"\tincluding gaps:";		
		for(int g=0; g<_characters[c]->getGapsIndices().size(); ++g){
				int gapNum = _characters[c]->getGapsIndices()[g];
				fileNexusStream<<" indel_"<<gapNum;
				//gapStream<<"\tgap "<<gapNum<<": "<<_gaps[gapNum-1]->getSeqID();
				//gapStream<<", "<<_gaps[gapNum-1]->getCoord5Abs()<<", "<<_gaps[gapNum-1]->getLength()<<";";
			}
			fileNexusStream<<endl;
			//gapStream<<endl;
	}

	fileNexusStream<<"MATRIX"<<endl<<endl;
	for(int s=0; s<_sc.numberOfSeqs();++s){
		fileNexusStream<<""<<_sc.name(s)<<"";	//prev- "Species_"
		for(c=0; c<_matrix[0].size(); ++c ){
			if(isSIC && _matrix[s][c]==2)
				fileNexusStream<<'?';
			else
				fileNexusStream<<_matrix[s][c];
		}
		fileNexusStream<<endl;	// prev- 2 endl
	}

	fileNexusStream<<";"<<endl<<endl;
	fileNexusStream<<"END;"<<endl<<endl<<endl;
	fileNexusStream<<"BEGIN ASSUMPTIONS; [below are the cost matrices of character change between the indel character state, these matrix exist only for characters that have more than two states]"
		<<endl<<endl<<endl<<endl;
	
	if(!isSIC){
		for(c=0; c<_characters.size();++c){
			int numOfStates = _characters[c]->getNumOfStates();
			if(numOfStates>2){
				fileNexusStream<<"[char "<<c+1+seqLeng<<", indel char "<<c+1<<" "<<"("<<_characters[c]->getCoord5()+1<<"-"<<_characters[c]->getCoord3()+1<<"):";
				fileNexusStream<<"0 (absent)";
				for(int st=1; st<numOfStates; ++st ){
					fileNexusStream<<", "<<st<<" (indel_"<<_characters[c]->getGapsIndices()[st] <<")";
				}
				fileNexusStream<<" ]"<<endl;

				fileNexusStream<<"usertype stepmatrix"<<c+1+seqLeng<<" (stepmatrix)="<<numOfStates<<endl<<endl;
				if(indelCoderOptions::_isCheckForTriangleInequality){
					if(_characters[c]->isTriangleInequalityCorrectionNeeded()){
						_characters[c]->printStepsMatrix(fileNexusStream,_characters[c]->isTriangleInequalityCorrectionNeeded());
						fileNexusStream<<"\n[prior to adjustment to satisfy triangle inequality:]\n";
						fileNexusStream<<"[";
						_characters[c]->printStepsMatrix(fileNexusStream);
						fileNexusStream<<"]\n";
					}
					else{
						_characters[c]->printStepsMatrix(fileNexusStream);
					}
				}
				else{
					_characters[c]->printStepsMatrix(fileNexusStream);
				}
				fileNexusStream<<";"<<endl<<endl<<endl;
			}
		}
		fileNexusStream<<endl<<endl<<endl;

		fileNexusStream<<"[below is the line that says to which indel correspond which matrix]"<<endl<<endl;

		fileNexusStream<<"typeset complexIndelCoding = ";
		for(c=0; c<_characters.size()-1;++c){
			int numOfStates = _characters[c]->getNumOfStates();
			if(numOfStates>2){
				fileNexusStream<<"stepmatrix"<<c<<" :"<<c<<",";
			}
		}
		fileNexusStream<<"stepmatrix"<<c<<" :"<<c<<";"<<endl<<endl;
		fileNexusStream<<";"<<endl<<endl;
		fileNexusStream<<"END;"<<endl<<endl<<endl<<endl;
	}

	
	fileNexusStream<<"BEGIN SETS;"<<endl<<endl;
	fileNexusStream<<"CHARSET indels = 1-"<<_characters.size()<<";"<<endl<<endl;
	fileNexusStream<<"END;"<<endl<<endl<<endl;
	fileNexusStream<<"BEGIN PAUP;"<<endl<<endl;
	fileNexusStream<<"assume typeset=complexIndelCoding;"<<endl<<endl;
	fileNexusStream<<"END;"<<endl<<endl;

	fileNexusStream<<"[Indels:"<<endl;
	fileNexusStream<<"No.	extension"<<endl;
	//for(int c=0; c<_gaps.numOfGaps();++c){
	//	fileNexusStream<<c+1<<"\t"<<_gaps[c]->getCoord5()<<"-"<<_gaps[c]->getCoord3()<<endl;		
	//}
	for(c=0; c<_characters.size();++c){
		fileNexusStream<<c+1<<"\t"<<_characters[c]->getCoord5()+1<<"-"<<_characters[c]->getCoord3()+1<<endl;		
	}
	fileNexusStream<<"]"<<endl;
}

void indelCoder::printIndelSummary()  {
	bool isSIC  = indelCoderOptions::_codingType == indelCoderOptions::SIC;
	//string fileGapString = indelCoderOptions::_outDir + "//" + "outGapInfoCHARSTATELABELS.txt";
	//ofstream gapStream(fileGapString.c_str());	
	ofstream gapStream(indelCoderOptions::_indelOutputInfoFile.c_str());

	
	gapStream<<"# each character is listing all included gaps.\n";
	gapStream<<"# Character: Start position relative to MSA (first pos of gap, count from 0); End position relative to MSA (+1);  length.\n";
	gapStream<<"# Gap: seqID(num start from 0); coord5Abs (relative to genome. not by MSA,first pos of gap, count from 0); length.\n";
	
	int seqLeng = _sc.seqLen();
	for(int c=0; c<_characters.size();++c){
		gapStream<<"character number: "<<c<<endl;
		gapStream<<"Start position relative to MSA: "<<_characters[c]->getCoord5()<<endl;
		gapStream<<"End position relative to MSA: "<<_characters[c]->getCoord3()+1<<endl;
		gapStream<<"Length: "<<_characters[c]->getCoord3()-_characters[c]->getCoord5()+1<<endl;

		vector<bool> isSpeciesWithGap(_sc.numberOfSeqs(),false);

		for(int g=0; g<_characters[c]->getGapsIndices().size(); ++g){
				int gapNum = _characters[c]->getGapsIndices()[g];
				int speciesSeqID = _gaps[gapNum-1]->getSeqID();
				isSpeciesWithGap[speciesSeqID] = true;
				gapStream<<"Found in species: "<<_sc.name(speciesSeqID);
				gapStream<<" Start position relative to genome: "<<_gaps[gapNum-1]->getCoord5Abs();
				gapStream<<" Length: "<<_gaps[gapNum-1]->getLength()<<endl;
				//gapStream<<"\tgap "<<gapNum<<": "<<_gaps[gapNum-1]->getSeqID();
				//gapStream<<", "<<_gaps[gapNum-1]->getCoord5Abs()<<", "<<_gaps[gapNum-1]->getLength()<<";";

		}
		gapStream<<"NOT FOUND in species: ";
		for(int i=0; i<_sc.numberOfSeqs(); ++i){
			if(!isSpeciesWithGap[i] && _matrix[i][c]!=2)
				gapStream<<_sc.name(i)<<",";
		}
		gapStream<<"\n";
		gapStream<<"ENDCHARACTER"<<endl<<endl;
	}
}

