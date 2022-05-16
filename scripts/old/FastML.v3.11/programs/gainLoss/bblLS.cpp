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
#include "bblLS.h"
#include "numRec.h"
#include "likelihoodComputation.h"
#include "likelihoodComputationGL.h"
#include "gainLossOptions.h"
#include <cmath>


bblLS::bblLS()
{}


MDOUBLE bblLS::optimizeBranches(tree& tr, stochasticProcess* sp, const sequenceContainer &sc, Vdouble* weights, unObservableData* unObservableData_p,
                                const int outerIter,
								const MDOUBLE epsilonOptimizationBranch, const int numIterations,
								MDOUBLE curL)
{
	_weights = weights;
	MDOUBLE prevIterL = VERYSMALL;
	if (curL == NULL)
		_treeLikelihood = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,*sp,_weights,unObservableData_p);

	else
		_treeLikelihood  = curL;
	LOGnOUT(4,<<"============================="<<endl;);
	LOGnOUT(4,<<"ll before bbl = "<<_treeLikelihood<<endl;);
	vector<tree::nodeP> nodesV;
	tr.getAllNodes(nodesV,tr.getRoot());
	int numberOfBranchs = nodesV.size();
	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/1.5; // (was 2) for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationIter = epsilonOptimizationBranch*epsilonOptimizationIterFactor;	// for eBranch=0.2 next iteration only for 10 logL points
	LOGnOUT(4,<<"BBL starts with epsilon branch= "<<epsilonOptimizationBranch<<" and epsilon iter="<<epsilonOptimizationIter<<endl;);
	int iter;
	for (iter = 1; iter <= numIterations; ++iter) 
	{
		if (_treeLikelihood < prevIterL + epsilonOptimizationIter){
			LOGnOUT(3,<<" BBL optimization converged. Iter= "<<iter<<" Likelihood="<<_treeLikelihood<<endl);
			return _treeLikelihood; //likelihood converged
		}
		prevIterL = _treeLikelihood;
		LOG(4,<<"---- BBL iteration: "<<iter<<endl;);
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<nodesV.size(); i++)
		{
			if (nodesV[i]->isRoot()) 
				continue;
			oldBl = nodesV[i]->dis2father();
			if(gainLossOptions::_isBblForceFactorCorrection){
				newL = -brent((oldBl+gainLossOptions::_minBranchLength)/gainLossOptions::_BblFactorCorrection, 
					oldBl, 
					(oldBl+gainLossOptions::_minBranchLength)*gainLossOptions::_BblFactorCorrection,
					evalBranch(nodesV[i],&tr, sc, sp,_weights,unObservableData_p), epsilonOptimizationBranch, &paramFound); 
			}
			else{
				newL = -brent(gainLossOptions::_minBranchLength, oldBl, gainLossOptions::_maxBranchLength, evalBranch(nodesV[i],&tr, sc, sp,_weights,unObservableData_p), epsilonOptimizationBranch, &paramFound); 
			}
			if (newL >= _treeLikelihood) 
			{
				_treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
				if(unObservableData_p) unObservableData_p->setLforMissingData(tr,sp);
				LOGnOUT(4,<<"BL old... "<<oldBl<<" BL done... "<<nodesV[i]->dis2father()<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				unObservableData_p->setLforMissingData(tr,sp);
				LOGnOUT(4,<<"*** WARNING: L went down : "<<endl;);
				LOGnOUT(4,<<" BL Found... "<<paramFound<<"...LL="<<newL<<"...";);
				LOGnOUT(4,<<" BL old... "<<oldBl<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			}			
		}		
		string treeINodes = gainLossOptions::_outDir + "//" + "TheTree.INodes.iter" +int2string(outerIter)+ ".Inner"+ int2string(iter) + ".ph"; 
		printTree (tr, treeINodes);
		LOGnOUT(3,<<"BBL iter "<<iter<<"...LL="<<_treeLikelihood<<"..."<<endl;);

	}
	if (iter>numIterations)
		LOGnOUT(4,<<" Too many="<<iter-1<<" iterations in BBL. Last optimized tree is used."<<endl);
	return _treeLikelihood;
}
//////////////////////////////////////////////////////////////////////////
MDOUBLE bblLS::optimizeBranches(tree& tr, vector<vector<stochasticProcess*> >& spVVec,
						 const distribution * gainDist, const distribution * lossDist,
						 const sequenceContainer &sc, 
						 Vdouble* weights, unObservableData* unObservableData_p,
						 const int outerIter,
						 const MDOUBLE epsilonOptimizationBranch , const int numIterations ,
						 MDOUBLE curL) 	
{
	_weights = weights;
	MDOUBLE prevIterL = VERYSMALL;
	if (curL == NULL)
		_treeLikelihood = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(tr,sc,spVVec,gainDist,lossDist,weights,unObservableData_p);
	else
		_treeLikelihood  = curL;
	LOGnOUT(4,<<"============================="<<endl;);
	LOGnOUT(4,<<"ll before bbl = "<<_treeLikelihood<<endl;);
	vector<tree::nodeP> nodesV;
	tr.getAllNodes(nodesV,tr.getRoot());
	int numberOfBranchs = nodesV.size();
	MDOUBLE epsilonOptimizationIterFactor = numberOfBranchs/2.0; // for 100 branches (~50 species) the epsilon for the entire iter is 50 times the one for branch
	epsilonOptimizationIterFactor = max(5.0,epsilonOptimizationIterFactor);
	MDOUBLE epsilonOptimizationIter = epsilonOptimizationBranch*epsilonOptimizationIterFactor;	// for eBranch=0.2 next iteration only for 10 logL points
	LOGnOUT(4,<<"BBL starts with epsilon branch= "<<epsilonOptimizationBranch<<" and epsilon iter="<<epsilonOptimizationIter<<endl;);
	int iter;
	for (iter = 1; iter <= numIterations; ++iter) 
	{
		if (_treeLikelihood < prevIterL + epsilonOptimizationIter){
			LOGnOUT(3,<<" BBL optimization converged. Iter= "<<iter<<" Likelihood="<<_treeLikelihood<<endl);
			return _treeLikelihood; //likelihood converged
		}
		prevIterL = _treeLikelihood;
		LOG(4,<<"---- BBL iteration: "<<iter<<endl;);
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<numberOfBranchs; i++)
		{
			if (nodesV[i]->isRoot()) 
				continue;
			oldBl = nodesV[i]->dis2father();
			if(gainLossOptions::_isBblForceFactorCorrection){
				newL = -brent((oldBl+gainLossOptions::_minBranchLength)/gainLossOptions::_BblFactorCorrection, 
					oldBl, 
					(oldBl+gainLossOptions::_minBranchLength)*gainLossOptions::_BblFactorCorrection, evalBranchSPvv(nodesV[i],&tr, sc, spVVec,gainDist,lossDist,weights,unObservableData_p), epsilonOptimizationBranch, &paramFound); 
			}
			else{
				newL = -brent(gainLossOptions::_minBranchLength, oldBl, gainLossOptions::_maxBranchLength, evalBranchSPvv(nodesV[i],&tr, sc, spVVec,gainDist,lossDist,weights,unObservableData_p), epsilonOptimizationBranch, &paramFound); 
			}
			if (newL >= _treeLikelihood) 
			{
				_treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
				if(unObservableData_p)	unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"BL old... "<<oldBl<<" BL done... "<<nodesV[i]->dis2father()<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				if(unObservableData_p)	unObservableData_p->setLforMissingData(tr,spVVec,gainDist,lossDist);
				LOGnOUT(4,<<"*** WARNING: L went down: "<<endl;);
				LOGnOUT(4,<<" BL Found... "<<paramFound<<"...LL="<<newL<<"...";);
				LOGnOUT(4,<<" BL old... "<<oldBl<<"...LL="<<_treeLikelihood<<"..."<<endl;);
			}			
		}
		string treeINodes = gainLossOptions::_outDir + "//" + "TheTree.INodes.iter" +int2string(outerIter)+ ".Inner"+ int2string(iter) + ".ph";
		printTree (tr, treeINodes);
		LOGnOUT(3,<<"BBL iter "<<iter<<"...LL="<<_treeLikelihood<<"..."<<endl;);
	}
	if (iter>numIterations)
		LOGnOUT(4,<<" Too many="<<iter-1<<" iterations in BBL. Last optimized tree is used."<<endl);
	return _treeLikelihood;
}


//////////////////////////////////////////////////////////////////////////
MDOUBLE evalBranch::operator()(MDOUBLE x)
{
	_pNode->setDisToFather(x);
	if(_unObservableData_p)_unObservableData_p->setLforMissingData(*_tr,_sp);
	MDOUBLE LL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_tr,_sc,*_sp,_weights,_unObservableData_p);
	return -LL;
}

//////////////////////////////////////////////////////////////////////////
MDOUBLE evalBranchSPvv::operator()(MDOUBLE x)
{
	_pNode->setDisToFather(x);
	if(_unObservableData_p)	_unObservableData_p->setLforMissingData(*_tr,_spVVec,_gainDist,_lossDist);
	MDOUBLE LL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(*_tr,_sc,_spVVec,_gainDist,_lossDist,_weights,_unObservableData_p);
	return -LL;
}
//////////////////////////////////////////////////////////////////////////
MDOUBLE evalBranchProportionExponent::operator()(MDOUBLE x)
{
	
	MDOUBLE factorBL = pow(10,x);
	_tr->multipleAllBranchesByFactor(factorBL);
	if(_unObservableData_p)_unObservableData_p->setLforMissingData(*_tr,_sp);
	MDOUBLE LL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(*_tr,_sc,*_sp,_weights,_unObservableData_p);
	_tr->multipleAllBranchesByFactor(1/factorBL);
	LOG(5,<<"Branch factor val = "<<factorBL<<" logL = "<<LL<<endl);
	return -LL;
}

//////////////////////////////////////////////////////////////////////////
MDOUBLE evalBranchProportionExponentSPvv::operator()(MDOUBLE x)
{
	MDOUBLE factorBL = pow(10,x);
	_tr->multipleAllBranchesByFactor(factorBL);
	if(_unObservableData_p)	_unObservableData_p->setLforMissingData(*_tr,_spVVec,_gainDist,_lossDist);
	MDOUBLE LL = likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(*_tr,_sc,_spVVec,_gainDist,_lossDist,_weights,_unObservableData_p);
	LOG(5,<<"Branch factor val = "<<factorBL<<" logL = "<<LL<<endl);
	_tr->multipleAllBranchesByFactor(1/factorBL);
	return -LL;
}


