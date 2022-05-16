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
#include "sankoffReconstructGL.h"
#include "threeStateAlphabet.h"
#include "treeIt.h"
#include "matrixUtils.h"
#include "sequence.h"
#include "someUtil.h"
#include "recognizeFormat.h"
#include "seqContainerTreeMap.h"
#include "treeUtil.h"
#include "amino.h"
#include "nucleotide.h"
#include "integerAlphabet.h"
#include "logFile.h"
#include "gainLossOptions.h"

sankoffReconstructGL::sankoffReconstructGL(sequenceContainer& sc, tree& tr, string& outDir, MDOUBLE costMatrixGainLossRatio, MDOUBLE distanceFromRootForRecent):
_tr(tr),_sc(sc),_outDir(outDir),_costMatrixGainLossRatio(costMatrixGainLossRatio),_distanceFromRootForRecent(distanceFromRootForRecent)
{
	initialize();
	//myLog::setLog(MPoptions::_logfile, 5);
	run();
}

sankoffReconstructGL::~sankoffReconstructGL(){
	//if (_alph)
	//	delete _alph;
}

/********************************************************************************************
initialize
*********************************************************************************************/
void sankoffReconstructGL::initialize() {
	//string paramStr = argv[1];
	//MPoptions::initOptions(paramStr);
	//startTree();
	//startSequenceContainer();
	_states.resize(_tr.getNodesNum(),-1000);
	_gainMPPerPos.resize(_sc.seqLen());
	_lossMPPerPos.resize(_sc.seqLen());
	 resizeVVV(_sc.seqLen(),_sc.alphabetSize(),_sc.alphabetSize(), _MPPerPos);
	 resizeVVV(_tr.getNodesNum(),_sc.alphabetSize(),_sc.alphabetSize(), _MPPerBranch);
	 resizeVVVV(_sc.seqLen(),_tr.getNodesNum(),_sc.alphabetSize(),_sc.alphabetSize(), _MPPerPosPerNode);	
	 startCostMatrix();
	_costOfTree = 0.0;
	_numOfGains = 0;
	_numOfLosses = 0;
}

/********************************************************************************************
More functions
*********************************************************************************************/
//void sankoffReconstructGL::startTree(){
//	tree t(MPoptions::_treefile);
//	_tr = t;
//	if (!(MPoptions::_rootAt =="")){
//		tree::nodeP myroot = _tr.findNodeByName(MPoptions::_rootAt); //returns NULL if not found
//		if (myroot){
//			_tr.rootAt(myroot);
//		}
//		else {
//			errorMsg::reportError("Requested root name is not found");
//		}
//	}
//	else {
//		LOGnOUT(5,<<"Default rooting used, root name is "<<_tr.getRoot()->name()<<endl);
//	}
//	LOGnOUT(5,<<"sons of root are "<<endl);
//	for (int son=0; son<_tr.getRoot()->getNumberOfSons();++son){
//		LOGnOUT(5,<<_tr.getRoot()->getSon(son)->name()<<endl);
//	}
//}

/********************************************************************************************
*********************************************************************************************/
//void sankoffReconstructGL::startSequenceContainer(){
//	switch (MPoptions::_alphabetType) {
//		case (MPoptions::amino): 
//			_alph = new amino; break;
//		case (MPoptions::nuc):
//			_alph = new nucleotide; break;
//		case (MPoptions::integer):
//			_alph = new integerAlphabet(MPoptions::_alphabetSize); break;
//		case (MPoptions::threeState): default:
//			_alph = new threeStateAlphabet; break;
//
//
//	}
//	string strFile = MPoptions::_seqfile;
//	ifstream in(strFile.c_str());
//	_sc = recognizeFormat::read(in,_alph);
//	_states.resize(_tr.getNodesNum(),-1000);
//	checkThatNamesInTreeAreSameAsNamesInSequenceContainer(_tr,_sc);
//
//}

/********************************************************************************************
*********************************************************************************************/
void sankoffReconstructGL::startCostMatrix(){
	switch (gainLossOptions::_costMatrixType) {
		case (gainLossOptions::file):
			// if specified an input cost matrix:
			if (gainLossOptions::_costMatrixfile != "") {
				readMatrixFromFile(_costMatrix,gainLossOptions::_costMatrixfile);
				if (_costMatrix.size() != _sc.alphabetSize()) {
					errorMsg::reportError("error in sankoff::startCostMatrix, the cost matrix must be the same size as the alphabet");
				}
			} else 
				errorMsg::reportError("error in sankoff::startCostMatrix, the cost matrix file is not specified after the -mf flag");
			break;
		case (gainLossOptions::diff) :
			resizeMatrix(_costMatrix,_sc.alphabetSize(),_sc.alphabetSize());
			for (int row=0; row<_costMatrix.size(); ++row){
				for (int col=0; col<_costMatrix.size(); ++col){
					_costMatrix[row][col]=(row-col);
				}
			}
			break;
		case (gainLossOptions::diffSquare) :
			resizeMatrix(_costMatrix,_sc.alphabetSize(),_sc.alphabetSize());
			for (int row=0; row<_costMatrix.size(); ++row){
				for (int col=0; col<_costMatrix.size(); ++col){
					_costMatrix[row][col]=(row-col)*(row-col);
				}
			}
			break;
		case (gainLossOptions::gainLossCost) :
			resizeMatrix(_costMatrix,_sc.alphabetSize(),_sc.alphabetSize());
			for (int row=0; row<_costMatrix.size(); ++row){
				for (int col=0; col<_costMatrix.size(); ++col){
					MDOUBLE cost = (row==col? 0: (row<col?_costMatrixGainLossRatio:1)); //gain 2(or other set value), loss 1
					_costMatrix[row][col]=cost;
				}
			}
			break;
		case (gainLossOptions::fitch) : default:
			//default: equal cost matrix (Fitch)
			resizeMatrix(_costMatrix,_sc.alphabetSize(),_sc.alphabetSize());
			for (int row=0; row<_costMatrix.size(); ++row){
				for (int col=0; col<_costMatrix.size(); ++col){
					_costMatrix[row][col]=(row==col?0.0:1.0);
				}
			}
			break;
	}	
}


/********************************************************************************************
*********************************************************************************************/
void sankoffReconstructGL::run(){
	//ofstream oStream((MPoptions::_outfile).c_str());
	//oStream<<"Maximum parsimony reconstruction"<<endl;
	//oStream<<"For each position, the reconstructed tree is presetned with reconstructed data as BP at each node, followed by a matrix specifying the number of each transition found"<<endl<<endl;
	//oStream.close();
	LOGnOUT(4,<<" MaxParsimony with costMatrix - gainLossRatio 1:"<<_costMatrixGainLossRatio<<endl);
	
	string MPprints =  _outDir + "//" + "MPprints." + double2string(_costMatrixGainLossRatio)+ ".txt"; 
	ofstream MPprintsStream(MPprints.c_str());
	MPprintsStream<<"# Various MP prints: "<<endl;	
	
	string gainLossMPPerPosPerBranch =  _outDir + "//" + "gainLossMP." + double2string(_costMatrixGainLossRatio)+ ".PerPosPerBranch.txt"; 
	ofstream gainLossMPPerPosPerBranchStream(gainLossMPPerPosPerBranch.c_str());
	gainLossMPPerPosPerBranchStream<<"# print with MP based on the cost matrix: "<<endl;

	// per pos
	string gainLossMPPerPos =  _outDir + "//" + "gainLossMP." + double2string(_costMatrixGainLossRatio)+ ".PerPos.txt"; 
	ofstream gainLossMPPerPosStream(gainLossMPPerPos.c_str());
	gainLossMPPerPosStream<<"# print with MP based on the cost matrix: "<<endl;	

	// per branch
	string gainLossMPPerBranch =  _outDir + "//" + "gainLossMP." + double2string(_costMatrixGainLossRatio)+ ".PerBranch.txt"; 
	ofstream gainLossMPPerBranchStream(gainLossMPPerBranch.c_str());
	gainLossMPPerBranchStream<<"# print with MP based on the cost matrix: "<<endl;

	// state per node, Sankoff reconstruction
	string gainLossMPAncestralReconstruct =  _outDir + "//" + "gainLossMP." + double2string(_costMatrixGainLossRatio)+ ".AncestralReconstructSankoff.txt"; 
	ofstream gainLossMPAncestralReconstructStream(gainLossMPAncestralReconstruct.c_str());
	gainLossMPAncestralReconstructStream<<"# print with MP based on the cost matrix: "<<endl;	
	
	for (int row=0; row<_costMatrix.size(); ++row){
		for (int col=0; col<_costMatrix.size(); ++col){
			MPprintsStream<<"#  "<<row<<"->"<<col<<" ="<<_costMatrix[row][col]<<endl;
			gainLossMPPerPosPerBranchStream<<"#  "<<row<<"->"<<col<<" ="<<_costMatrix[row][col]<<endl; 
			gainLossMPPerPosStream<<"#  "<<row<<"->"<<col<<" ="<<_costMatrix[row][col]<<endl; 
			gainLossMPPerBranchStream<<"#  "<<row<<"->"<<col<<" ="<<_costMatrix[row][col]<<endl;
			gainLossMPAncestralReconstructStream<<"#  "<<row<<"->"<<col<<" ="<<_costMatrix[row][col]<<endl;
		}
	}
	
	gainLossMPPerPosPerBranchStream<<"G/L"<<"\t"<<"POS"<<"\t"<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2Root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"probability"<<"\t"<<"expectation"<<endl;
	gainLossMPAncestralReconstructStream<<"POS"<<"\t"<<"Node"<<"\t"<<"State"<<endl;
	for (int pos = 0 ; pos<_sc.seqLen(); ++pos) {
		LOGnOUT(7,<<"Running position "<<pos+1<<endl<<"=========================="<<endl);	
		_costOfTree += runPosition(pos, gainLossMPPerPosPerBranchStream, MPprintsStream,gainLossMPAncestralReconstructStream);
	}
	LOGnOUT(4,<<"Cost of tree is "<<_costOfTree<<" (with "<<_numOfGains+_numOfLosses<<" events)"<<endl);
	LOGnOUT(4,<<" Gain="<<_numOfGains<<endl);
	LOGnOUT(4,<<" Losses="<<_numOfLosses<<endl);

	// Per Branch
	printMPPerBranch(gainLossMPPerBranchStream);

	// Per Pos
	printMPPerPos(gainLossMPPerPosStream);

}

/********************************************************************************************
*********************************************************************************************/
// transitionTypeCount is printed to the outfile only if the _costMatrixType != diffSquare or diff
// totalCosts is printed to the outfile only if _costMatrixType == diffSquare or diff
MDOUBLE sankoffReconstructGL::runPosition(int pos, ofstream& gainLossMPPerPosPerBranchStream, ofstream& MPprints, ofstream& gainLossMPAncestralReconstructStream){
	// intialize _states veactor with values of leaves
	seqContainerTreeMap scTreeMap(_sc,_tr);	
	vector <tree::nodeP> leaves;
	_tr.getAllLeaves(leaves,_tr.getRoot());
	for (int i=0; i< leaves.size();i++){
		int myleafId = (leaves[i])->id();
		int mySeqId = scTreeMap.seqIdOfNodeI(myleafId);
		_states[myleafId] = _sc[mySeqId][pos];
	}
	VVdouble upcosts;
	vector <VVint> backtrack;
	traverseUpMP(upcosts, backtrack);
	//ofstream oStream((MPoptions::_outfile).c_str(),ios::app);
	//oStream<<"======================================"<<endl;
	//oStream<<"POSITION "<<pos<<endl;
	//oStream<<"======================================"<<endl;
	VVint transitionTypeCount;
	VVdouble totalCosts;

	MDOUBLE costoftree = traverseDownMP(upcosts, backtrack, transitionTypeCount,totalCosts);
	Vstring data;
	preparePrintData(data);
	//printDataOnTreeAsBPValues(oStream,data,_tr);
	//oStream<<endl<<"Cost of tree is "<<costoftree<<endl<<"Transition type count:"<<endl;
	LOGnOUT(7,<<"Cost of position "<< pos+1 <<" is "<<costoftree<<endl<<endl);


	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		gainLossMPAncestralReconstructStream<<pos+1<<"\t"<<mynode->name()<<"\t"<<_states[mynode->id()]<<endl;
		if(mynode->isRoot())
			continue;
		int stateAtNode = _states[mynode->id()];
		int stateAtFather = _states[mynode->father()->id()];
		if(stateAtNode > stateAtFather){
			gainLossMPPerPosPerBranchStream<<"gain"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<"1"<<"\t"<<stateAtNode-stateAtFather<<endl;
			_gainMPPerPos[pos]++;
			_MPPerPos[pos][0][1]++;
			_MPPerBranch[mynode->id()][0][1]++;
			_MPPerPosPerNode[pos][mynode->id()][0][1]++;
			_numOfGains++;
		}
		if(stateAtNode < stateAtFather){
			gainLossMPPerPosPerBranchStream<<"loss"<<"\t"<<pos+1<<"\t"<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<"1"<<"\t"<<-(stateAtNode-stateAtFather)<<endl;
			_lossMPPerPos[pos]++;			
			_MPPerPos[pos][1][0]++;
			_MPPerBranch[mynode->id()][1][0]++;
			_MPPerPosPerNode[pos][mynode->id()][1][0]++;
			_numOfLosses++;
		}
	}

	if ((gainLossOptions::_costMatrixType  != gainLossOptions::diffSquare) &&
		(gainLossOptions::_costMatrixType  != gainLossOptions::diff) )		
	{
		for (int i = 0; i < transitionTypeCount.size(); i++) {
			for (int j = 0; j < transitionTypeCount[i].size(); j++) {
				MPprints<<transitionTypeCount[i][j]<<" ";
			}
			MPprints<<endl;
		}
	} else {
		for (int i=0; i< totalCosts.size();++i) {
			MPprints << "node " << i ;
			if (_tr.findNodeById(i)->isLeaf())
				MPprints << " (leaf)" ;
			if (_tr.findNodeById(i)->isRoot())
				MPprints << " (root)" ;
			MPprints <<" :" << endl ;
			for (int j=0; j < _costMatrix.size();++j)
				MPprints<< totalCosts[i][j] << " ";
			MPprints << endl;
		}
	}
	return costoftree;
}

/********************************************************************************************
*********************************************************************************************/
void sankoffReconstructGL::traverseUpMP(VVdouble &upCosts, vector <VVint> &backtrack) { 
	// upCosts[i][j] i for node, j for size of cost matrix
	// backtrack[i][j][k] remembers the state for which a min was obtained for node i, state j, from both sons (k=0 and k=1)
	int i;
	gainLossAlphabet alph;
	upCosts.resize(_tr.getNodesNum());
	for (i = 0; i < upCosts.size(); i++) 
		upCosts[i].resize(_costMatrix.size(),0.0);
	backtrack.resize(_tr.getNodesNum());
	for (i = 0; i < backtrack.size(); i++) {
		backtrack[i].resize(_costMatrix.size());
	}
	
	// fill upCosts, starting with leafs (0,Inf) according to the observed character
	treeIterDownTopConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf()) {
			for (int j = 0; j < _costMatrix.size(); j++) {
				upCosts[mynode->id()][j] = ( (_states[mynode->id()] == j || _states[mynode->id()] == alph.unknown() ) ? 0 : VERYBIG);
			}
		}
		else {
			for (int k = 0; k < _costMatrix.size(); k++) { // this loop fills each cell in the vector for node mynode
				for (int son=0; son<mynode->getNumberOfSons(); ++son) {	// go over all sons
					MDOUBLE minSon = VERYBIG;
					int argMinSon=-1;	// for backtrack
					int idSon = (mynode->getSon(son))->id();
					
					//for (int l = _costMatrix.size()-1; l >= 0; l--) { // loop to find the min, 1 is preferred
					for (int l = 0; l < _costMatrix.size(); l++) { // loop to find the min, 0 is preferred
						MDOUBLE sumSon = upCosts[idSon][l]+_costMatrix[k][l];
						if ( sumSon < minSon) {
							minSon = sumSon;
							argMinSon = l;
						}
					}
					if ((argMinSon==-1) || (minSon==VERYBIG)){
						errorMsg::reportError("Error in sankoff::traverseUpMP, unknown reason");
					}

					upCosts[mynode->id()][k]+=minSon;
					backtrack[mynode->id()][k].push_back(argMinSon);
				}
			}
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
// totalCosts is only filled for _costMatrixType==diffSquare or diff
MDOUBLE sankoffReconstructGL::traverseDownMP(VVdouble &upCosts, vector <VVint> &backtrack, 
											 VVint &transitionTypeCount,VVdouble &totalCosts) {
	if (upCosts.size() == 0) 
		errorMsg::reportError("error in sankoff::traverseDownMP, input vector upCosts must be filled (call traverseUpMP() first)");
	if (backtrack.size() == 0) 
		errorMsg::reportError("error in sankoff::traverseDownMP, input vector backtrack must be filled (call traverseUpMP() first)");
	int sizeOfCosts = upCosts[0].size();
	totalCosts.resize(_tr.getNodesNum());
	for (int i = 0; i < totalCosts.size(); i++) {
		totalCosts[i].resize(_costMatrix.size(),0.0);
	}

	MDOUBLE costOfTree = 0;
	int stateOfRoot;
	findMinInVector(upCosts[(_tr.getRoot())->id()], costOfTree, stateOfRoot);	// first, reconstruct Root
	_states[(_tr.getRoot())->id()] = stateOfRoot;
	
	transitionTypeCount.resize(sizeOfCosts);
	for (int i = 0; i < transitionTypeCount.size(); i++) 
		transitionTypeCount[i].resize(sizeOfCosts,0);
	
	treeIterTopDownConst tIt(_tr);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf()) continue;
		int myId = mynode->id();
		for (int j=0; j<mynode->getNumberOfSons(); ++j) {
			int idSon = (mynode->getSon(j))->id();
			_states[idSon] = backtrack[myId][_states[myId]][j];
			transitionTypeCount[_states[myId]][_states[idSon]]++;
			if ((gainLossOptions::_costMatrixType  == gainLossOptions::diffSquare) ||
				(gainLossOptions::_costMatrixType  == gainLossOptions::diff)){
					for (int z=0; z <_costMatrix.size(); ++z) // go over all the states
						totalCosts[idSon][z] = upCosts[idSon][z] + _costMatrix[_states[myId]][z];
				}
		}
		// fill totalCosts of the root
		if (mynode->isRoot()) {
			if ((gainLossOptions::_costMatrixType  == gainLossOptions::diffSquare) ||
				(gainLossOptions::_costMatrixType  == gainLossOptions::diff)){
					for (int z=0; z <_costMatrix.size(); ++z) // go over all the states
						totalCosts[myId][z] = upCosts[myId][z];
				}
		}
	}
	return costOfTree;
}

/********************************************************************************************
*********************************************************************************************/
//prepares the data to be printed as BP data on the tree
void sankoffReconstructGL::preparePrintData(Vstring &data){
	data.resize(_tr.getNodesNum());
	for (int i=0; i< data.size(); ++i) {
		data[i] = double2string(_states[i]);
		data[i]+="[";
		data[i]+=_tr.findNodeById(i)->name();
		data[i]+="]";
	}
}


/********************************************************************************************
*********************************************************************************************/
void sankoffReconstructGL::printMPPerBranch(ostream& out)
{
	treeIterTopDownConst tIt(_tr);
	out<<"# MP Gain and Loss counts"<<"\n";
	out<<"branch"<<"\t"<<"branchLength"<<"\t"<<"distance2root"<<"\t"<<"distance2NearestOTU"<<"\t"<<"numOfNodes2NearestOTU"<<"\t"<<"exp01"<<"\t"<<"exp10"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		out<<mynode->name()<<"\t"<<mynode->dis2father()<<"\t"<<mynode->getDistance2ROOT()<<"\t"<<mynode->getMinimalDistance2OTU()<<"\t"<<mynode->getMinimalNumOfNodes2OTU()<<"\t"<<_MPPerBranch[mynode->id()][0][1]<<"\t"<<_MPPerBranch[mynode->id()][1][0]<<endl;
	}
}


/********************************************************************************************
printProbExp()
print perPos (over all branches)
use the members _expV01, _expV10 for basic 
*********************************************************************************************/
void sankoffReconstructGL::printMPPerPos(ostream& out)
{
	out<<"POS"<<"\t"<<"MP01"<<"\t"<<"MP10"<<endl;
	for (int pos = 0; pos <_sc.seqLen(); ++pos){
		out<<pos+1<<"\t"<<_MPPerPos[pos][0][1]<<"\t"<<_MPPerPos[pos][1][0]<<endl;
	}
}

