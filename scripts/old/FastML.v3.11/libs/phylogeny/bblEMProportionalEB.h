// $Id: bblEMProportional.h 962 2006-11-07 15:13:34Z privmane $
#ifndef ___BBL_EM_PROPORTIONALEB_H
#define ___BBL_EM_PROPORTIONALEB_H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "multipleStochasticProcess.h"

#include <vector>
using namespace std;


class bblEMProportionalEB {
public:
	explicit bblEMProportionalEB(tree& et,
									const vector<sequenceContainer>& sc,
									multipleStochasticProcess* msp,
									const gammaDistribution* pProportionDist,
									const bool optimizeSelectedBranches=false,
									const vector<Vdouble *> * weights = NULL,
									const int maxIterations=50,
									const MDOUBLE epsilon=0.05,
									const MDOUBLE tollForPairwiseDist=0.0001,
									const MDOUBLE* likelihoodLast=NULL);
	Vdouble getTreeLikelihood() const {return _treeLikelihoodVec;}

private:
	Vdouble compute_bblEMPropEB(const int maxIterations,const MDOUBLE epsilon,const MDOUBLE tollForPairwiseDist,const MDOUBLE* likelihoodLast=NULL);
	void allocatePlacePropEB();
	void computeUpPropEB();
	void bblEM_itPropEB(const MDOUBLE tollForPairwiseDist);
	void computeDownPropEB(const int gene, const int pos);
	void addCountsPropEB(const int gene, const int pos);
	void addCountsPropEB(const int gene,const int pos, tree::nodeP mynode, const VdoubleRep posProb, const MDOUBLE weig);
	void optimizeBranchesPropEB(const MDOUBLE tollForPairwiseDist);

	Vdouble _treeLikelihoodVec;
	tree& _et;
	const vector<sequenceContainer>& _sc;
	multipleStochasticProcess* _msp;
	const gammaDistribution* _pProportionDist;
	const vector<Vdouble *> * _weights;
	int _numberOfGenes;
	vector<	vector<countTableComponentGamProportional> > _computeCountsV; // for each gene, for each node - a table of globalRate*localRate*alph*alph - [globalRateCategory][localRateCategory][character]
	vector<suffStatGlobalGamProportional> _cup; //[gene][pos][globalRateCategory][localRateCategory][nodeID][character]
	vector<suffStatGlobalGamProportionalPos> _cdown; //[gene][globalRateCategory][localRateCategory][nodeID][character]
	vector< vector<computePijGam> > _pij;//[gene][globalRateCategory]
	VVVdoubleRep _posLike;//[gene][pos][globalRateCategory]
	const bool _optimizeSelectedBranches;

};

#endif
