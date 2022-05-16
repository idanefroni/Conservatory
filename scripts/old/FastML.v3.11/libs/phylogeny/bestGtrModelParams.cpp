// $Id: bestGtrModelparams.cpp 2008-29-04 10:57:00Z nimrod $

#include "bestGtrModelParams.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "bblEMProportionalEB.h"
#include "bblLSProportionalEB.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestGtrModel::bestGtrModel(tree& et, // find best Gtr Model Params
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights,
										const int maxTotalIterations,
										const MDOUBLE epsilonLikelihoodImprovment,
										const MDOUBLE epsilonLoglikelihoodForGTRParam,
										const MDOUBLE upperBoundGTRParam,
										const bool optimizeTree,
										const bool optimizeAlpha){
	LOG(5,<<"Starting bestGtrModel: find Best replacement matrix parameters"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	_bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(et,sc,sp,weights);

	MDOUBLE prev_a2c = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2c();
	MDOUBLE prev_a2g = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2g();
	MDOUBLE prev_a2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2t();
	MDOUBLE prev_c2g = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_c2g();
	MDOUBLE prev_c2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_c2t();
	MDOUBLE prev_g2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_g2t();

	MDOUBLE prevAlpha = epsilonLoglikeForBBL;

	for (int i=0; i < maxTotalIterations; ++i) {
		//optimize a2c
		newL = -brent(0.0, prev_a2c, upperBoundGTRParam,
					  C_evalGTRParam(a2c,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2c);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2c(_best_a2c);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2c(prev_a2c);
			LOG(5,<<"likelihood went down in optimizing a2c"<<endl<<"oldL = "<<_bestL);
		}

		//optimize a2t
		newL = -brent(0.0, prev_a2t, upperBoundGTRParam,
					  C_evalGTRParam(a2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2t(_best_a2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2t(prev_a2t);
			LOG(5,<<"likelihood went down in optimizing a2t"<<endl<<"oldL = "<<_bestL);
		}

		//optimize a2g
		newL = -brent(0.0, prev_a2g, upperBoundGTRParam,
					  C_evalGTRParam(a2g,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2g);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2g(_best_a2g);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2g(prev_a2g);
			LOG(5,<<"likelihood went down in optimizing a2g"<<endl<<"oldL = "<<_bestL);
		}

		//optimize c2g
		newL = -brent(0.0, prev_c2g, upperBoundGTRParam,
					  C_evalGTRParam(c2g,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_c2g);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2g(_best_c2g);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2g(prev_c2g);
			LOG(5,<<"likelihood went down in optimizing c2g"<<endl<<"oldL = "<<_bestL);
		}

		//optimize c2t
		newL = -brent(0.0, prev_c2t, upperBoundGTRParam,
					  C_evalGTRParam(c2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_c2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2t(_best_c2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2t(prev_c2t);
			LOG(5,<<"likelihood went down in optimizing c2t"<<endl<<"oldL = "<<_bestL);
		}

		//optimize g2t
		newL = -brent(0.0, prev_g2t, upperBoundGTRParam,
					  C_evalGTRParam(g2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_g2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_g2t(_best_g2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_g2t(prev_g2t);
			LOG(5,<<"likelihood went down in optimizing g2t"<<endl<<"oldL = "<<_bestL);
		}
		if(optimizeAlpha)
		{
			newL = -brent(0.0, prevAlpha, upperBoundForAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonLoglikeForAlphaOptimization,
					  &_bestAlpha);
			(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha); 

			if (newL >= _bestL) 
			{
				_bestL = newL;
				(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha); //safety
			} 
			else
			{//likelihood went down!
	            (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(prevAlpha);
				LOG(5,<<"likelihood went down in optimizing alpha"<<endl<<"oldL = "<<_bestL);
			}
		}

		if(optimizeTree)
		{
			bblEM bblEM1(et,sc,sp,weights,maxBBLIt,epsilonLoglikeForBBL);
			_bestL = bblEM1.getTreeLikelihood();
		}


		// check for improvement in the likelihood
		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
			prev_a2c = _best_a2c;
			prev_a2g = _best_a2g;
			prev_a2t = _best_a2t;
			prev_c2g = _best_c2g;
			prev_c2t = _best_c2t;
			prev_g2t = _best_g2t;
			prevAlpha = _bestAlpha;
		} else {
			break;
		}
	}
}

bestGtrModelProportional::bestGtrModelProportional(tree& et, // find best Gtr Model Params under a proportional model
										vector<sequenceContainer>& sc,
										multipleStochasticProcess* msp,
 										gammaDistribution* pProportionDist,
										Vdouble initLocalAlphas,
										Vdouble initLocala2cs,
										Vdouble initLocala2gs,
										Vdouble initLocala2ts,
										Vdouble initLocalc2gs,
										Vdouble initLocalc2ts,
										Vdouble initLocalg2ts,
										const MDOUBLE upperBoundOnLocalAlpha,
										const MDOUBLE initGlobalAlpha,
										const MDOUBLE upperBoundOnGlobalAlpha,
										const MDOUBLE upperBoundGTRParam,
										const int maxTotalIterations,
										const int maxBBLIterations,
										const bool optimizeSelectedBranches,
										const bool optimizeTree,
										const string branchLengthOptimizationMethod,
										const bool optimizeLocalParams,
										const bool optimizeGlobalAlpha,
										const Vdouble * weights,
										const MDOUBLE epsilonLikelihoodImprovment,
										const MDOUBLE epsilonLoglikelihoodForGTRParam,
										const MDOUBLE epsilonLoglikelihoodForLocalAlphaOptimization,
										const MDOUBLE epsilonLoglikelihoodForGlobalAlphaOptimization,
										const MDOUBLE epsilonLoglikelihoodForBBL){
	LOG(5,<<"Starting bestGtrModelProportional"<<endl);
	Vdouble current_a2cVec,current_a2gVec,current_a2tVec,current_c2gVec,current_c2tVec,current_g2tVec,currentLocalAlphaVec;
	MDOUBLE currentGlobalAlpha = initGlobalAlpha;
	currentLocalAlphaVec = initLocalAlphas;
	current_a2cVec = initLocala2cs;
	current_a2gVec = initLocala2gs;
	current_a2tVec = initLocala2ts;
	current_c2gVec = initLocalc2gs;
	current_c2tVec = initLocalc2ts;
	current_g2tVec = initLocalg2ts;

	Vdouble newLvec;
	//doubleRep epsilonLoglikelihoodForGlobalAlphaOptimizationDR(epsilonLoglikelihoodForGlobalAlphaOptimization);//DR
	newLvec.resize(msp->getSPVecSize());
	//doubleRep oldL(VERYSMALL);//DR
	//doubleRep newL;//DR
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL;
	_bestLvec.resize(msp->getSPVecSize(),0.0);
	_bestLocalAlphaVec = initLocalAlphas;
	_bestGlobalAlpha = initGlobalAlpha;
	int spIndex;
	_best_a2cVec = current_a2cVec;
	_best_a2gVec = current_a2gVec;
	_best_a2tVec = current_a2tVec;
	_best_c2gVec = current_c2gVec;
	_best_c2tVec = current_c2tVec;
	_best_g2tVec = current_g2tVec;
	pProportionDist->setAlpha(_bestGlobalAlpha);
	for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
		(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2c(_best_a2cVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2g(_best_a2gVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2t(_best_a2tVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_c2g(_best_c2gVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_c2t(_best_c2tVec[spIndex]);
		(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_g2t(_best_g2tVec[spIndex]);
	}
	//first compute the likelihood;
	_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);

	MDOUBLE ax_local = 0.0;
	MDOUBLE c_GTRParam_x = upperBoundGTRParam;
	MDOUBLE c_localAlpha_x = upperBoundOnLocalAlpha;

	for (int i=0; i < maxTotalIterations; ++i) {
		if(optimizeLocalParams){
			for(spIndex = 0;spIndex < msp->getSPVecSize();++spIndex){
				//optimize a2c
				MDOUBLE a2c_x = _best_a2cVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,a2c_x,c_GTRParam_x,
						  C_evalGTRParamProportional(a2c,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_a2cVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_a2cVec[spIndex] = current_a2cVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing a2c"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2c(_best_a2cVec[spIndex]);//safety
				
				//optimize a2t
				MDOUBLE a2t_x = _best_a2tVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,a2t_x,c_GTRParam_x,
						  C_evalGTRParamProportional(a2t,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_a2tVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_a2tVec[spIndex] = current_a2tVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing a2t"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2t(_best_a2tVec[spIndex]);//safety

				//optimize a2g
				MDOUBLE a2g_x = _best_a2gVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,a2g_x,c_GTRParam_x,
						  C_evalGTRParamProportional(a2g,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_a2gVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_a2gVec[spIndex] = current_a2gVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing a2g"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_a2g(_best_a2gVec[spIndex]);//safety

				//optimize c2g
				MDOUBLE c2g_x = _best_c2gVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,c2g_x,c_GTRParam_x,
						  C_evalGTRParamProportional(c2g,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_c2gVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_c2gVec[spIndex] = current_c2gVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing c2g"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_c2g(_best_c2gVec[spIndex]);//safety

				//optimize c2t
				MDOUBLE c2t_x = _best_c2tVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,c2t_x,c_GTRParam_x,
						  C_evalGTRParamProportional(c2t,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_c2tVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_c2tVec[spIndex] = current_c2tVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing c2t"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_c2t(_best_c2tVec[spIndex]);//safety

				//optimize g2t
				MDOUBLE g2t_x = _best_g2tVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,g2t_x,c_GTRParam_x,
						  C_evalGTRParamProportional(g2t,et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
						  epsilonLoglikelihoodForGTRParam,
						  &current_g2tVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_best_g2tVec[spIndex] = current_g2tVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing g2t"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gtrModel*>(msp->getSp(spIndex)->getPijAccelerator()->getReplacementModel()))->set_g2t(_best_g2tVec[spIndex]);//safety

				//optimize local alpha
				MDOUBLE localAlpha_x = _bestLocalAlphaVec[spIndex];
				newLvec[spIndex] = -brent(ax_local,localAlpha_x,c_localAlpha_x,
					  C_evalLocalAlpha(et,sc[spIndex],*msp->getSp(spIndex),pProportionDist,weights),
					  epsilonLoglikelihoodForLocalAlphaOptimization,
					  &currentLocalAlphaVec[spIndex]);
				if (newLvec[spIndex] >= _bestLvec[spIndex]) 
				{
					_bestLvec[spIndex] = newLvec[spIndex];
					_bestLocalAlphaVec[spIndex] = currentLocalAlphaVec[spIndex];
				} 
				else
				{//likelihood went down!
					LOG(2,<<"likelihood went down in optimizing local alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
				}
				(static_cast<gammaDistribution*>(msp->getSp(spIndex)->distr()))->setAlpha(_bestLocalAlphaVec[spIndex]); //safety
			}
			LOGnOUT(2,<<"Done with GTR local params optimization"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Local Params:"<<endl);
			LOGnOUT(2,<<"a2c:");
			for(spIndex = 0;spIndex < _best_a2cVec.size();++spIndex){
				LOGnOUT(2,<<_best_a2cVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"a2g:");
			for(spIndex = 0;spIndex < _best_a2gVec.size();++spIndex){
				LOGnOUT(2,<<_best_a2gVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"a2t:");
			for(spIndex = 0;spIndex < _best_a2tVec.size();++spIndex){
				LOGnOUT(2,<<_best_a2tVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"c2g:");
			for(spIndex = 0;spIndex < _best_c2gVec.size();++spIndex){
				LOGnOUT(2,<<_best_c2gVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"c2t:");
			for(spIndex = 0;spIndex < _best_c2tVec.size();++spIndex){
				LOGnOUT(2,<<_best_c2tVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"g2t:");
			for(spIndex = 0;spIndex < _best_g2tVec.size();++spIndex){
				LOGnOUT(2,<<_best_g2tVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
			LOGnOUT(2,<<"local alpha:");
			for(spIndex = 0;spIndex < _bestLocalAlphaVec.size();++spIndex){
				LOGnOUT(2,<<_bestLocalAlphaVec[spIndex]<<",";);
			}
			LOGnOUT(2,<<endl);
		}
		if(optimizeGlobalAlpha){
			//doubleRep ax_global(0.0);//DR
			//doubleRep c_globalAlpha_x(upperBoundOnGlobalAlpha);//DR
			//doubleRep minusOne(-1.0);//DR
			MDOUBLE ax_global = 0.0;
			MDOUBLE c_globalAlpha_x = upperBoundOnGlobalAlpha;

			//optimize global alpha
			//doubleRep globalAlpha_x(prevGlobalAlpha);//DR
			MDOUBLE globalAlpha_x = _bestGlobalAlpha;
			//newL = minusOne*brentDoubleRep(ax_global,globalAlpha_x,c_globalAlpha_x,
			//		C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
			//		epsilonLoglikelihoodForGlobalAlphaOptimizationDR,
			//		&_bestGlobalAlpha);
			newL = -brent(ax_global,globalAlpha_x,c_globalAlpha_x,
					C_evalGlobalAlpha(et,sc,msp,pProportionDist,weights),
					epsilonLoglikelihoodForGlobalAlphaOptimization,
					&currentGlobalAlpha);
			if (newL >= sumVdouble(_bestLvec))
			{
				_bestGlobalAlpha = currentGlobalAlpha;
			} 
			else
			{//likelihood went down!
				LOG(2,<<"likelihood went down in optimizing global alpha"<<endl<<"oldL = "<<sumVdouble(_bestLvec));
			}
			pProportionDist->setAlpha(_bestGlobalAlpha); //safety
			//whether or not likelihood has improved we need to update _bestLvec 
			_bestLvec = likelihoodComputation::getTreeLikelihoodProportionalAllPosAlphTheSame(et,sc,msp,pProportionDist,weights);
			LOGnOUT(2,<<"Done with global alpha optimization"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			LOGnOUT(2,<<"Global Alpha:"<<_bestGlobalAlpha<<endl);
		}

		if(optimizeTree)
		{
			if(branchLengthOptimizationMethod == "bblLS"){
				bblLSProportionalEB bblLSPEB1(et,sc,msp,pProportionDist,_bestLvec,optimizeSelectedBranches,maxBBLIterations,epsilonLoglikelihoodForBBL);
				_bestLvec = bblLSPEB1.getTreeLikelihoodVec();
				LOGnOUT(2,<<"Done with bblLS"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			else if(branchLengthOptimizationMethod == "bblEM"){
				bblEMProportionalEB bblEMPEB1(et,sc,msp,pProportionDist,optimizeSelectedBranches,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);
				_bestLvec = bblEMPEB1.getTreeLikelihood();
				LOGnOUT(2,<<"Done with bblEM"<<endl<<"LL:"<<sumVdouble(_bestLvec)<<endl);
			}
			LOGnOUT(2,<<et.stringTreeInPhylipTreeFormat()<<endl);
		}
		// check for improvement in the likelihood
		if (sumVdouble(_bestLvec) > oldL+epsilonLikelihoodImprovment) {
			//all params have already been updated
			oldL = sumVdouble(_bestLvec);
		} else {
			break;
		}
		LOGnOUT(4,<<"Done with optimization iteration "<<i<<". LL: "<<sumVdouble(_bestLvec)<<endl);
	}
}

