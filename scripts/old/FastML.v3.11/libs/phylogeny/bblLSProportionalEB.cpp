#include "bblLSProportionalEB.h"
#include "numRec.h"
#include "logFile.h"
#include "errorMsg.h"


bblLSProportionalEB::bblLSProportionalEB(tree& et, const vector<sequenceContainer>& sc, multipleStochasticProcess* msp, const gammaDistribution* pProportionDist, Vdouble& treeLikelihoodVec, const bool optimizeSelectedBranches, int maxIter, MDOUBLE epsilon)
{
	_treeLikelihoodVec = optimizeBranches(et,sc,msp,pProportionDist,treeLikelihoodVec,optimizeSelectedBranches,maxIter,epsilon);
}


Vdouble bblLSProportionalEB::optimizeBranches(tree& et, const vector<sequenceContainer>& sc, multipleStochasticProcess* msp, const gammaDistribution* pProportionDist, Vdouble& inTreeLikelihoodVec, const bool optimizeSelectedBranches, int maxIter, MDOUBLE epsilon)
{
	Vdouble treeLikelihoodVec;
	if (inTreeLikelihoodVec.empty()){
        treeLikelihoodVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist);
	}
	else{
		treeLikelihoodVec  = inTreeLikelihoodVec;
	}
	MDOUBLE treeLikelihood = sumVdouble(treeLikelihoodVec);
	LOGnOUT(5,<<"ll before bblLSr4sp"<<" logL="<<treeLikelihood<<endl);
	vector<tree::nodeP> nodesV;
	et.getAllNodes(nodesV,et.getRoot());
	MDOUBLE prevIterL = VERYSMALL;
	for (int iter = 0; iter < maxIter; ++iter) 
	{
		if (treeLikelihood < prevIterL + epsilon){
			treeLikelihoodVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist);
			return treeLikelihoodVec; //likelihood converged
		}
		prevIterL = treeLikelihood;
		MDOUBLE paramFound;
		MDOUBLE oldBl;
		MDOUBLE newL;
		for (int i=0; i<nodesV.size(); i++)
		{
			if (nodesV[i]->isRoot()) continue;
			if((optimizeSelectedBranches) && (nodesV[i]->getComment() != "1")) continue; //only selected branhes will be optimized
			oldBl = nodesV[i]->dis2father();			
			newL = -brent(0.0,oldBl,MAX_BRANCH_LENGTH,evalR4SPBranch(nodesV[i],et,sc,msp,pProportionDist),epsilon,&paramFound);
			LOGnOUT(4,<<"oldL="<<treeLikelihood<<" newL="<<newL<<" BL="<<nodesV[i]->dis2father()<<endl);
			if (newL >= treeLikelihood) 
			{
				treeLikelihood = newL;
				nodesV[i]->setDisToFather(paramFound);
			} 
			else //likelihood went down!
			{
				nodesV[i]->setDisToFather(oldBl); //return to previous BL
				LOGnOUT(4,<<"Likelihood went down. oldL="<<treeLikelihood<<" newL="<<newL<<" Do not update tree"<<endl);
			}
		}
	}
	treeLikelihoodVec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist);
	return treeLikelihoodVec;
}

