// $Id: bblEM.cpp 4478 2008-07-17 17:09:55Z cohenofi $
#include "bblEMfixRoot.h"
#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistancefixRoot.h"
#include <ctime>

bblEMfixRoot::bblEMfixRoot(tree& et,
				const sequenceContainer& sc,
				const stochasticProcess& sp,
				const Vdouble * weights,
				const int maxIterations,
				const MDOUBLE epsilon,
				const MDOUBLE tollForPairwiseDist,
				unObservableData*  unObservableData_p,
				const MDOUBLE* likelihoodLast) :
_et(et),_sc(sc),_sp(sp),_weights (weights),_unObservableData_p(unObservableData_p) 
{
	//if(!plogLforMissingData){
	//	_plogLforMissingData = NULL;
	//}
	_treeLikelihood = compute_bblEM(maxIterations,epsilon,tollForPairwiseDist,likelihoodLast);
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE bblEMfixRoot::compute_bblEM(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist,
			const MDOUBLE* likelihoodLast){
	allocatePlace();
	MDOUBLE oldL=VERYSMALL;
	MDOUBLE currL = VERYSMALL;
	tree oldT = _et;
	for (int i=0; i < maxIterations; ++i) {
		//if(_unObservableData_p)
		//	_unObservableData_p->setLforMissingData(_et,&_sp);

		computeUp();
		currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights,_unObservableData_p);
		LOGnOUT(4,<<"--- Iter="<<i<<" logL="<<currL<<endl);
		//if(_unObservableData_p){
			//if(!likelihoodLast)
			//	LOGnOUT(4,<<" WARNING!!! likelihoodLast was not sent to bblEM with unObservableData prog"<<endl);				
		if(oldL<=currL){	// make sure not to use tree with lower likelihood then last computed likelihood (before BBL-EM)
				if(likelihoodLast){
					if(*likelihoodLast<=currL)
						oldT = _et;	 // L didn't go down
					else
						LOGnOUT(4,<<"Likelihood went down compared pre-BBL oldL="<<*likelihoodLast<<" newL="<<currL<<" Do not update tree"<<endl);
				}
				else{
					oldT = _et;	 // L didn't go down
					LOGnOUT(7,<<"LikelihoodLast was not sent to bblEM"<<endl);
				}
		}
		else
			LOGnOUT(4,<<"Likelihood went down oldL="<<oldL<<" newL="<<currL<<" Do not update tree"<<endl);				
		//}
		//else
			//oldT = _et;	// all application that don't need to correct for unObservableData, update tree
		
		
		if (currL < oldL + epsilon) { // need to break
			if (currL<=oldL) {
				_et = oldT;
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_et,&_sp);
				return oldL; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				return currL;
			}
		}
        bblEM_it(tollForPairwiseDist);
		oldL = currL;
	}
	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUp();
	if(_unObservableData_p)
		_unObservableData_p->setLforMissingData(_et,&_sp);
	currL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights, _unObservableData_p);

	if (currL<=oldL) 
	{
		_et = oldT;
		if(_unObservableData_p)
			_unObservableData_p->setLforMissingData(_et,&_sp);
		return oldL; // keep the old tree, and old likelihood
	} 
	else 
        return currL;
}

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::allocatePlace() {
	
	_computeCountsV.resize(_et.getNodesNum());//initiateTablesOfCounts
	for (int node=0; node < _computeCountsV.size(); ++node) {
	{
		_computeCountsV[node].resize(_sp.alphabetSize()); 
		for (int letterAtRoot = 0; letterAtRoot < _computeCountsV[node].size(); ++letterAtRoot)
			_computeCountsV[node][letterAtRoot].countTableComponentAllocatePlace(_sp.alphabetSize(),_sp.categories()); //_computeCountsV[node][letterAtRoot][rate][alph][alph]
			//_computeCountsV[i][letterAtRoot].zero();	// removed, a BUG, done later	
		}
	}
	
	_cup.allocatePlace(_sc.seqLen(),_sp.categories(), _et.getNodesNum(), _sc.alphabetSize());	
	_cdown.resize(_sp.categories());
	for (int categor = 0; categor < _sp.categories(); ++categor)
	{
		// stay with the convention of fillComputeDownNonReversible where the first index is for rate cat and the second is for letterAtRoot
		_cdown[categor].allocatePlace(_sp.alphabetSize(), _et.getNodesNum(), _sc.alphabetSize()); //_cdown[categ][letter@root][nodeid][letter][prob]
	}
}

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::bblEM_it(const MDOUBLE tollForPairwiseDist){
	string costTable =  "costTableBBLEMit.txt";	//DEBUG 
	ofstream costTableStream(costTable.c_str());	//DEBUG
	//cout<<"before zero\n";
	for (int node=0; node < _computeCountsV.size(); ++node) {
		for (int letAtRoot=0; letAtRoot < _computeCountsV[node].size(); ++letAtRoot) {
			_computeCountsV[node][letAtRoot].zero();
			_computeCountsV[node][letAtRoot].printTable(costTableStream);	//DEBUG
		}
	}
	//cout<<"after zero\n";

	for (int i=0; i < _sc.seqLen(); ++i) {
		computeDown(i);
		addCounts(i); // computes the counts and adds to the table.
	}
	//cout<<"after add counts\n";
	for (int node=0; node < _computeCountsV.size(); ++node) {
		for (int letAtRoot=0; letAtRoot < _computeCountsV[node].size(); ++letAtRoot) {
			_computeCountsV[node][letAtRoot].printTable(costTableStream);	//DEBUG
		}
	}


	optimizeBranches(tollForPairwiseDist);
	if(_unObservableData_p)
		_unObservableData_p->setLforMissingData(_et,&_sp);
}

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::optimizeBranches(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistancefixRoot from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father(),_unObservableData_p);
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());

			if(false){	//DEBUG
				if(_unObservableData_p)
					_unObservableData_p->setLforMissingData(_et,&_sp);
				computeUp();
				MDOUBLE bL = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc,_sp,_cup,_posLike,_weights, _unObservableData_p);
				LOG(6,<<"  node "<<mynode->name()<<" L= "<<bL<<endl);
			}
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::computeUp(){
	_pij.fillPij(_et,_sp,0); // 0 is becaues we compute Pij(t) and not its derivations...
	computeUpAlg cupAlg;
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
        for (int categor = 0; categor < _sp.categories(); ++categor) {
			cupAlg.fillComputeUp(_et,_sc,pos,_pij[categor],_cup[pos][categor]);
		}
	}
 }

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::computeDown(const int pos){
	computeDownAlg cdownAlg;
		for (int categor = 0; categor < _sp.categories(); ++categor) {
			cdownAlg.fillComputeDownNonReversible(_et,_sc,pos,_pij[categor],_cdown[categor],_cup[pos][categor]);
		}
 }

/********************************************************************************************
*********************************************************************************************/
void bblEMfixRoot::addCounts(const int pos){
	//MDOUBLE posProb = 
	//	likelihoodComputation::getProbOfPosWhenUpIsFilledGam(pos,_et,_sc,_sp,_cup);
						
	MDOUBLE weig = (_weights ? (*_weights)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCountsFixedRoot(pos,mynode,_posLike[pos],weig);
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
// fill _computeCountsV: specific node, letterAtRoot and categor at a time
void bblEMfixRoot::addCountsFixedRoot(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig){

	computeCounts cc;
	for(int letterAtRoot = 0; letterAtRoot < _sp.alphabetSize(); letterAtRoot++)
	{
		for (int categor =0; categor< _sp.categories(); ++ categor) 
		{
				cc.computeCountsNodeFatherNodeSonHomPos(_sc,
											_pij[categor],
											_sp,
											_cup[pos][categor],
											_cdown[categor][letterAtRoot],
											weig,
											posProb,
											mynode,
											_computeCountsV[mynode->id()][letterAtRoot][categor],
											_sp.ratesProb(categor),
											letterAtRoot);	// letterInFather is used in freq? or already by _cdown?
		}
	}
}          

