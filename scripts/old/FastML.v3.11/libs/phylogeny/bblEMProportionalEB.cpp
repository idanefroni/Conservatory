// $Id: bblEMProprtional.cpp 962 2006-11-07 15:13:34Z privmane $
#include "bblEM.h"
#include "bblEMProportionalEB.h"
#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "computeUpAlg.h"
#include "computeDownAlg.h"
#include "computeCounts.h"
#include "treeIt.h"
#include "fromCountTableComponentToDistance.h"
#include <ctime>//#define VERBOS
#include "fromCountTableComponentToDistancePropEB.h"

bblEMProportionalEB::bblEMProportionalEB(tree& et,
									const vector<sequenceContainer>& sc,
									multipleStochasticProcess* msp,
									const gammaDistribution* pProportionDist,
									const bool optimizeSelectedBranches,
									const vector<Vdouble *> * weights,
									const int maxIterations,
									const MDOUBLE epsilon,
									const MDOUBLE tollForPairwiseDist,
									const MDOUBLE* likelihoodLast):

_et(et),_sc(sc),_msp(msp),_pProportionDist(pProportionDist),_weights (weights),_optimizeSelectedBranches(optimizeSelectedBranches) {
	_numberOfGenes = _sc.size();
	assert(_msp->getSPVecSize() == _sc.size());
	_treeLikelihoodVec = compute_bblEMPropEB(maxIterations,epsilon,tollForPairwiseDist,likelihoodLast);
}

Vdouble bblEMProportionalEB::compute_bblEMPropEB(
			const int maxIterations,
			const MDOUBLE epsilon,
			const MDOUBLE tollForPairwiseDist,
			const MDOUBLE* likelihoodLast){
	LOGnOUT(5,<<"Allocating place"<<endl);
	allocatePlacePropEB();
	LOGnOUT(5,<<"Done Allocating place"<<endl);
	Vdouble oldLvec(_numberOfGenes,VERYSMALL);
	Vdouble currLvec;
	currLvec.resize(_numberOfGenes);
	tree oldT = _et;
	//doubleRep epsilonDR(epsilon);//DR
	for (int i=0; i < maxIterations; ++i) {
		LOGnOUT(5,<<"Calling computeUpPropEB on iteration "<<i<<endl);
		computeUpPropEB();
		for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
			currLvec[geneN] = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc[geneN],*_msp->getSp(geneN),_cup[geneN],_pProportionDist,_posLike[geneN],(_weights?(*_weights)[geneN]:NULL));
		}
		LOGnOUT(5,<<"--- Iter="<<i<<" logL="<<sumVdouble(currLvec)<<endl);
		if(sumVdouble(oldLvec)<=sumVdouble(currLvec)){	// make sure not to use tree with lower likelihood then last computed likelihood (before BBL-EM)
			LOGnOUT(4,<<"Likelihood improved. oldL = "<<sumVdouble(oldLvec)<<" newL = "<<sumVdouble(currLvec)<<endl);
			if(likelihoodLast){
				if(*likelihoodLast<=sumVdouble(currLvec))
					oldT = _et;	 // L didn't go down
				else
					LOGnOUT(4,<<"Likelihood went down compared pre-BBL oldL="<<*likelihoodLast<<" newL="<<sumVdouble(currLvec)<<" Do not update tree"<<endl);
			}
			else{
				oldT = _et;	 // L didn't go down
				LOGnOUT(7,<<"Tree Updated"<<endl);
			}
		}
		else 
			LOGnOUT(4,<<"Likelihood did not improve. oldL="<<sumVdouble(oldLvec)<<" newL="<<sumVdouble(currLvec)<<" Do not update tree"<<endl);				

		if (sumVdouble(currLvec) < sumVdouble(oldLvec) + epsilon) { // need to break
			if (sumVdouble(currLvec)<sumVdouble(oldLvec)) {
				_et = oldT; //return to older tree
				LOGnOUT(4,<<"Finished bblEMPropEB. Likelihood ="<<sumVdouble(oldLvec)<<endl);
				return oldLvec; // keep the old tree, and old likelihood
			} else {
                //update the tree and likelihood and return
				LOGnOUT(4,<<"Finished bblEMPropEB. Likelihood ="<<sumVdouble(currLvec)<<endl);
				return currLvec;
			}
		}
		bblEM_itPropEB(tollForPairwiseDist);
		oldLvec = currLvec;
	}
	// in the case were we reached max_iter, we have to recompute the likelihood of the new tree...
	computeUpPropEB();
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		currLvec[geneN] = likelihoodComputation::getTreeLikelihoodFromUp2(_et,_sc[geneN],*_msp->getSp(geneN),_cup[geneN],_pProportionDist,_posLike[geneN],(_weights?(*_weights)[geneN]:NULL));
	}
	if (sumVdouble(currLvec)<sumVdouble(oldLvec)) {
		_et = oldT;
		LOGnOUT(4,<<"Finished bblEMPropEB max iter. Likelihood ="<<sumVdouble(oldLvec)<<endl);
		return oldLvec; // keep the old tree, and old likelihood
	} 
	else{
		LOGnOUT(4,<<"Finished bblEMPropEB max iter. Likelihood ="<<sumVdouble(currLvec)<<endl);
        return currLvec;
	}
}

void bblEMProportionalEB::allocatePlacePropEB() {
	_computeCountsV.resize(_numberOfGenes);
	_cup.resize(_numberOfGenes);
	_cdown.resize(_numberOfGenes);
	_pij.resize(_numberOfGenes);
	_posLike.resize(_numberOfGenes);
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		_posLike[geneN].resize(_sc[geneN].seqLen());
		for(int pos = 0;pos < _sc[geneN].seqLen();++pos){
			_posLike[geneN][pos].resize(_pProportionDist->categories(),0.0);
		}
		stochasticProcess * sp = _msp->getSp(geneN);
		_pij[geneN].resize(_pProportionDist->categories());
		_computeCountsV[geneN].resize(_et.getNodesNum()); //initiateTablesOfCounts
		for (int i=0; i < _computeCountsV[geneN].size(); ++i) {
			_computeCountsV[geneN][i].countTableComponentAllocatePlace(sp->alphabetSize(),_pProportionDist->categories(),sp->categories());
		}
		_cup[geneN].allocatePlace(_sc[geneN].seqLen(),_pProportionDist->categories(),sp->categories(),_et.getNodesNum(), _sc[geneN].alphabetSize());
		_cdown[geneN].allocatePlace(_pProportionDist->categories(),sp->categories(),_et.getNodesNum(), _sc[geneN].alphabetSize());
	}
}

void bblEMProportionalEB::computeUpPropEB(){
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		for(int globalRateCategor = 0;globalRateCategor < _pProportionDist->categories();++globalRateCategor){
			_msp->getSp(geneN)->setGlobalRate(_pProportionDist->rates(globalRateCategor));
			_pij[geneN][globalRateCategor].fillPij(_et,*_msp->getSp(geneN),0); // 0 is becaues we compute Pij(t) and not its derivations...
			computeUpAlg cupAlg;
			for (int pos=0; pos < _sc[geneN].seqLen(); ++pos) {
				for (int localRateCategor = 0; localRateCategor < _msp->getSp(geneN)->categories(); ++localRateCategor) {
					cupAlg.fillComputeUp(_et,_sc[geneN],pos,_pij[geneN][globalRateCategor][localRateCategor],_cup[geneN][pos][globalRateCategor][localRateCategor]);
				}
			}
		}
	}
}

void bblEMProportionalEB::bblEM_itPropEB(const MDOUBLE tollForPairwiseDist){
	for (int geneN=0; geneN < _numberOfGenes; ++geneN) {
		for (int treeNode=0; treeNode < _computeCountsV[geneN].size(); ++treeNode) {
            _computeCountsV[geneN][treeNode].zero();
		}
		for (int pos=0; pos < _sc[geneN].seqLen(); ++pos) {
			computeDownPropEB(geneN,pos);
			addCountsPropEB(geneN,pos); // computes the counts and adds to the table.
		}
	}
	optimizeBranchesPropEB(tollForPairwiseDist);
}

void bblEMProportionalEB::computeDownPropEB(const int gene, const int pos){
	computeDownAlg cdownAlg;
	stochasticProcess * sp = _msp->getSp(gene);
	for (int globalRateCategor = 0; globalRateCategor < _pProportionDist->categories(); ++globalRateCategor) {
		for (int localRateCategor = 0; localRateCategor < sp->categories(); ++localRateCategor) {
			//no need to set the global rate for each sp cause it was already performed in computeUpPropEB
			cdownAlg.fillComputeDown(_et,_sc[gene],pos,_pij[gene][globalRateCategor][localRateCategor],_cdown[gene][globalRateCategor][localRateCategor],_cup[gene][pos][globalRateCategor][localRateCategor]);
		}
	}
}

void bblEMProportionalEB::addCountsPropEB(const int gene, const int pos){
	vector<MDOUBLE> * weightsOfGene = (_weights?(*_weights)[gene]:NULL);					
	MDOUBLE weig = (weightsOfGene ? (*weightsOfGene)[pos] : 1.0);
	if (weig == 0) return;
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			addCountsPropEB(gene,pos,mynode,_posLike[gene][pos],weig);
		}
	}
}

void bblEMProportionalEB::addCountsPropEB(const int gene,const int pos, tree::nodeP mynode, const VdoubleRep posProb, const MDOUBLE weig){
	computeCounts cc;
	stochasticProcess * sp = _msp->getSp(gene);
	for (int globalRateCategor =0; globalRateCategor< _pProportionDist->categories(); ++globalRateCategor) {
		for (int localRateCategor =0; localRateCategor < sp->categories(); ++localRateCategor) {
			//cc.computeCountsNodeFatherNodeSonHomPosProportionalEB(_sc[gene],
			//							_pij[gene][globalRateCategor][localRateCategor],
			//							*sp,
			//							_cup[gene][pos][globalRateCategor][localRateCategor],
			//							_cdown[gene][globalRateCategor][localRateCategor],
			//							weig,
			//							posProb,
			//							mynode,
			//							_computeCountsV[gene][mynode->id()][globalRateCategor][localRateCategor],
			//							_pProportionDist->ratesProb(globalRateCategor)*sp->ratesProb(localRateCategor));
			cc.computeCountsNodeFatherNodeSonHomPosProportionalEB(_sc[gene],
										_pij[gene][globalRateCategor][localRateCategor],
										*sp,
										_cup[gene][pos][globalRateCategor][localRateCategor],
										_cdown[gene][globalRateCategor][localRateCategor],
										weig,
										posProb,
										mynode,
										_computeCountsV[gene][mynode->id()][globalRateCategor][localRateCategor]);
		}
	}
}

/*
//tal's old implementation, where i think there's a bug cause he sends _computeCountsV[mynode->id()] to the 
//fromCountTableComponentToDistanceProp constructor, but the first dimension of _computeCountsV is the genes and
//the tree nodes is only the second dimension
void bblEMProportionalEB::optimizeBranchesPropEB(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			fromCountTableComponentToDistanceProp from1(_computeCountsV[mynode->id()],_sp,tollForPairwiseDist,mynode->dis2father());
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
		}
	}
}
*/

void bblEMProportionalEB::optimizeBranchesPropEB(const MDOUBLE tollForPairwiseDist){
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (!tIt->isRoot()) {
			if((_optimizeSelectedBranches) && (tIt->getComment() != "1")) continue; //only selected branhes will be optimized
			fromCountTableComponentToDistancePropEB from1(_computeCountsV,mynode->id(),_msp,_pProportionDist,tollForPairwiseDist,mynode->dis2father());
			from1.computeDistance();
			mynode->setDisToFather(from1.getDistance());
		}
	}
}

