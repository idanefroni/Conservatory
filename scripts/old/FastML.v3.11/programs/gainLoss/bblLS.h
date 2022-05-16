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


#ifndef ___BBL_LS__
#define ___BBL_LS__

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "unObservableData.h"
#include "gainLossUtils.h"

using namespace std;

//#define MAX_BRANCH_LENGTH 50.0	//20.0

/*
This class optimize the branches using "naive" line search methodology.
go over each branch and optimize it using brent.
In one iteration it optimze seperatly all branches.
This procedure continues until convergence is reached or until the maximum number of iteration is reached.
*/
class bblLS {
public:
	
	explicit bblLS();
	~bblLS() {};
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}


	MDOUBLE optimizeBranches(tree& tr, stochasticProcess* sp, const sequenceContainer &sc, Vdouble* weights, unObservableData* unObservableData_p,
		const int outerIter,
		const MDOUBLE epsilonOptimization =0.1, const int numIterations =10,
		MDOUBLE curL =NULL);

	MDOUBLE optimizeBranches(tree& tr, vector<vector<stochasticProcess*> >& spVVec,
		const distribution * gainDist, const distribution * lossDist,
		const sequenceContainer &sc, 
		Vdouble* weights, unObservableData* unObservableData_p,
		const int outerIter,
		const MDOUBLE epsilonOptimization =0.1, const int numIterations =10,
		MDOUBLE curL =NULL); 	


private:
	Vdouble* _weights; 
	MDOUBLE _treeLikelihood;
};

//////////////////////////////////////////////////////////////////////////
class evalBranch{
public:
	explicit evalBranch(tree::nodeP pNode, tree* tr, const sequenceContainer &sc, stochasticProcess* sp, Vdouble* weights, unObservableData* unObservableData_p )
		:_pNode(pNode),_tr(tr), _sc(sc), _sp(sp),_weights(weights)
	{
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;

	};	
	virtual ~evalBranch(){
		if(_unObservableData_p)	delete _unObservableData_p;
	}	
	
	MDOUBLE operator() (MDOUBLE x);

private:
	tree::nodeP _pNode;
	tree* _tr;
	const sequenceContainer& _sc;
	const stochasticProcess* _sp;
	Vdouble* _weights;
	unObservableData* _unObservableData_p;
};

//////////////////////////////////////////////////////////////////////////
class evalBranchSPvv{
public:
	explicit evalBranchSPvv(tree::nodeP pNode, tree* tr, const sequenceContainer &sc, vector<vector<stochasticProcess*> >& spVVec,
		const distribution * gainDist, const distribution * lossDist,
		Vdouble* weights, unObservableData* unObservableData_p)
		:_pNode(pNode),_tr(tr),_sc(sc),_spVVec(spVVec), _gainDist(gainDist), _lossDist(lossDist),_unObservableData_p(unObservableData_p),_weights(weights)
	{
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;
	};
	virtual ~evalBranchSPvv(){
		if(_unObservableData_p)	delete _unObservableData_p;
	}
	MDOUBLE operator() (MDOUBLE x);

private:
	tree::nodeP _pNode;
	tree* _tr;
	const sequenceContainer& _sc;
	const vector<vector<stochasticProcess*> >& _spVVec;
	const distribution * _gainDist;
	const distribution * _lossDist;
	Vdouble* _weights;
	unObservableData* _unObservableData_p;
};


//////////////////////////////////////////////////////////////////////////
class evalBranchProportionExponent{
public:
	explicit evalBranchProportionExponent(tree* tr, const sequenceContainer &sc, stochasticProcess* sp, Vdouble* weights, unObservableData* unObservableData_p )
		:_tr(tr), _sc(sc), _sp(sp),_weights(weights)
	{
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;

	};	
	virtual ~evalBranchProportionExponent(){
		if(_unObservableData_p)	delete _unObservableData_p;
	}	

	MDOUBLE operator() (MDOUBLE x);

private:	
	tree* _tr;
	const sequenceContainer& _sc;
	const stochasticProcess* _sp;
	Vdouble* _weights;
	unObservableData* _unObservableData_p;
};

//////////////////////////////////////////////////////////////////////////
class evalBranchProportionExponentSPvv{
public:
	explicit evalBranchProportionExponentSPvv(tree* tr, const sequenceContainer &sc, vector<vector<stochasticProcess*> >& spVVec,
		const distribution * gainDist, const distribution * lossDist,
		Vdouble* weights, unObservableData* unObservableData_p)
		:_tr(tr),_sc(sc),_spVVec(spVVec), _gainDist(gainDist), _lossDist(lossDist),_unObservableData_p(unObservableData_p),_weights(weights)
	{
		if(unObservableData_p)
			_unObservableData_p = unObservableData_p->clone();
		else
			_unObservableData_p = NULL;
	};
	virtual ~evalBranchProportionExponentSPvv(){
		if(_unObservableData_p)	delete _unObservableData_p;
	}
	MDOUBLE operator() (MDOUBLE x);

private:
	tree* _tr;
	const sequenceContainer& _sc;
	const vector<vector<stochasticProcess*> >& _spVVec;
	const distribution * _gainDist;
	const distribution * _lossDist;
	Vdouble* _weights;
	unObservableData* _unObservableData_p;
};




#endif
