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
#include "gainLossUtils.h"
#include "gainLossOptions.h"
#include "gainLossModel.h"
#include "gammaDistributionPlusInvariant.h"
#include "Parameters.h"
#include <cmath>

/********************************************************************************************
*********************************************************************************************/
void printProgramInfo(){
	LOGnOUT(4,<<"+=============================================================+"<<endl);
	LOGnOUT(4,<<"+  The gainLoss project:										"<<endl);
	LOGnOUT(4,<<"+  Analysis of Phyletic Patterns in a Likelihood Framework		"<<endl);
	LOGnOUT(4,<<"+  "<<PROG_INFO<<"						"<<endl);
	LOGnOUT(4,<<"+  Ofir Cohen - ofircohe@tau.ac.il								"<<endl);
	LOGnOUT(4,<<"+=============================================================+"<<endl);
}
/********************************************************************************************
*********************************************************************************************/

void printTree (tree &tr, string treeFile){
	ofstream treeINodesStream(treeFile.c_str());
	printTree(tr,treeINodesStream);	//@ gainLoss - the branch lengths are lost	
}
void printTree (tree &tr){
	string treeINodes_st = gainLossOptions::_outDir + "//" + "TheTree.INodes.ph"; 
	ofstream treeINodesStream(treeINodes_st.c_str());
	printTree(tr,treeINodesStream);	//Note: @ gainLoss method - the branch lengths are lost	

	//string tree_st = gainLossOptions::_outDir + "//" + "TheTree.ph"; 
	//ofstream treeStream(tree_st.c_str());
	ofstream treeStream(gainLossOptions::_treeOutFile.c_str());
	tr.output(treeStream);
	treeStream.close();
}

void printTree (tree &tr,ostream &out){
	vector<tree::nodeP> vAllNodes;
	tr.getAllNodes(vAllNodes,tr.getRoot());
	Vstring Vnames(vAllNodes.size());
	for (int i = 0; i<vAllNodes.size();++i)
		Vnames[vAllNodes[i]->id()] = vAllNodes[i]->name();
	printTreeWithValuesAsBP(out,tr,Vnames);
	out<<endl;
}

/********************************************************************************************
*********************************************************************************************/
void printTreeWithValuesAsBP(ostream &out, tree &tr, Vstring values, VVVdouble *probs, bool printGains) {
	printTreeWithValuesAsBP(out,tr.getRoot(), values,probs);
	out<<"["<<values[tr.getRoot()->id()]<<"];";
}
void printTreeWithValuesAsBP(ostream &out, const tree::nodeP &myNode, Vstring values,  VVVdouble *probs, bool printGains)  {
	if (myNode->isLeaf()) {
		out<< myNode->name();
		if(probs)
			if (printGains)
				out<<"_P_"<<(*probs)[myNode->id()][0][1]; 
			else //print losses
				out<<"_P_"<<(*probs)[myNode->id()][1][0];
		out<< ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printTreeWithValuesAsBP(out, myNode->getSon(i), values,probs);
		}
		out <<")";
		if (myNode->isRoot()==false) {
			out<< myNode->name();
			if(probs)
				if (printGains)
					out<<"_P_"<<(*probs)[myNode->id()][0][1]; 
				else //print losses
					out<<"_P_"<<(*probs)[myNode->id()][1][0];
			out<< ":"<<myNode->dis2father();
//			out << "["<<values[myNode->id()]<<"]"; 
		}
	}
}

/********************************************************************************************
used For AncestralRec
*********************************************************************************************/
void printTreeStatesAsBPValues(ostream &out, Vint &states, tree &tr, 
VVVdouble *probs,bool printGains) {
	printTreeStatesAsBPValues(out,states, tr.getRoot(), probs);
	out<<"["<<(tr.getRoot())->name()<<"-"<<states[(tr.getRoot())->id()]<<"];";
}

void printTreeStatesAsBPValues(ostream &out, Vint &states, const tree::nodeP &myNode,	VVVdouble *probs,bool printGains)  {
	if (myNode->isLeaf()) {
		out << myNode->name()<<"-"<<states[myNode->id()]<< ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printTreeStatesAsBPValues(out,states,myNode->getSon(i),probs);
		}
		out <<")";
		if (myNode->isRoot()==false) {
			out.precision(3);
			if (probs){
				if (printGains)
					out<<(*probs)[myNode->id()][0][2]<<"//"<<(*probs)[myNode->id()][1][2]; 
				else //print losses
					out<<(*probs)[myNode->id()][2][0]<<"//"<<(*probs)[myNode->id()][2][1]; 
			}
			out << "["<<myNode->name()<<"-"<<states[myNode->id()]<<"]"; 
			out<<":"<<myNode->dis2father();
		}
	}
}

/********************************************************************************************
used For AncestralRec - Double (posterior size)
*********************************************************************************************/
void printTreeStatesAsBPValues(ostream &out, Vdouble &states, tree &tr,  
	 VVVdouble *probs, bool printGains) 
{
	printTreeStatesAsBPValues(out,states, tr.getRoot(), probs);
	out<<"["<<(tr.getRoot())->name()<<"-"<<states[(tr.getRoot())->id()]<<"];";
}

void printTreeStatesAsBPValues(ostream &out, Vdouble &states, const tree::nodeP &myNode, 
							   VVVdouble *probs,bool printGains)  
{
	if (myNode->isLeaf()) {
		out << myNode->name()<<"-"<<states[myNode->id()]<< ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			printTreeStatesAsBPValues(out,states,myNode->getSon(i),probs);
		}
		out <<")";
		if (myNode->isRoot()==false) {
			out.precision(3);
			if (probs){
				if (printGains)
					out<<(*probs)[myNode->id()][0][2]<<"//"<<(*probs)[myNode->id()][1][2]; 
				else //print losses
					out<<(*probs)[myNode->id()][2][0]<<"//"<<(*probs)[myNode->id()][2][1]; 
			}
			out << "["<<myNode->name()<<"-"<<states[myNode->id()]<<"]"; 
			out<<":"<<myNode->dis2father();
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE factorial (MDOUBLE num){
	if (num==1)
		return 1.0;
	return factorial(num-1)*num; 
}


/********************************************************************************************
*********************************************************************************************/
MDOUBLE getRateAlpha(distribution* dist)
{
	MDOUBLE res;
	//switch (gainLossOptions::_rateDistributionType)
	//{
	//case (gainLossOptions::GAMMA_PLUS_INV):
	//	res = static_cast<gammaDistributionPlusInvariant*>(dist)->getAlpha();
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):
	//	res = static_cast<generalGammaDistributionPlusInvariant*>(dist)->getAlpha();
	//	break;
	//case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
	//	res = static_cast<gammaDistributionFixedCategories*>(dist)->getAlpha();
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
	//	res = static_cast<generalGammaDistributionFixedCategories*>(dist)->getAlpha();
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA):
	//	res = static_cast<generalGammaDistribution*>(dist)->getAlpha();
	//	break;
	//case (gainLossOptions::GAMMA):
	//	res = static_cast<gammaDistribution*>(dist)->getAlpha();
	//	break;
	//default:
	//	errorMsg::reportError("unknown type in gainLossOptions::getDistributionType");
	//}
	if(dynamic_cast<gammaDistributionPlusInvariant*>(dist)){
		res = static_cast<gammaDistributionPlusInvariant*>(dist)->getAlpha();
	}
	else if(dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)){
		res = static_cast<generalGammaDistributionPlusInvariant*>(dist)->getAlpha();
	}
	else if (dynamic_cast<gammaDistributionFixedCategories*>(dist)){
		res = static_cast<gammaDistributionFixedCategories*>(dist)->getAlpha();
	}
	else if (dynamic_cast<generalGammaDistributionFixedCategories*>(dist)){
		res = static_cast<generalGammaDistributionFixedCategories*>(dist)->getAlpha();
	}
	else if (dynamic_cast<generalGammaDistribution*>(dist)){
		res = static_cast<generalGammaDistribution*>(dist)->getAlpha();
	}
	else if (dynamic_cast<gammaDistribution*>(dist)){
		res = static_cast<gammaDistribution*>(dist)->getAlpha();
	}
	else{
		LOGnOUT(4,<<"unknown type in gainLossOptions::getDistributionType, zero is filled for Alpha\n");
		res = 0;
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
void setRateAlpha(distribution* dist, MDOUBLE paramAlpha)
{	
	//switch (gainLossOptions::_rateDistributionType)
	//{
	//case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):
	//	static_cast<generalGammaDistributionPlusInvariant*>(dist)->setAlpha(paramAlpha);
	//	break;
	//case (gainLossOptions::GAMMA_PLUS_INV):
	//	static_cast<gammaDistributionPlusInvariant*>(dist)->setAlpha(paramAlpha);
	//	break;
	//case (gainLossOptions::GAMMA_FIXED_CATEGORIES):
	//	static_cast<gammaDistributionFixedCategories*>(dist)->setAlpha(paramAlpha);
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
	//	static_cast<generalGammaDistributionFixedCategories*>(dist)->setAlpha(paramAlpha);
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA):
	//	static_cast<generalGammaDistribution*>(dist)->setAlpha(paramAlpha);
	//	break;
	//case (gainLossOptions::GAMMA):
	//	static_cast<gammaDistribution*>(dist)->setAlpha(paramAlpha);
	//	break;
	//default:
	//	errorMsg::reportError("unknown type in distributionType");
	//}
	if (dynamic_cast<gammaDistributionPlusInvariant*>(dist)){
		static_cast<gammaDistributionPlusInvariant*>(dist)->setAlpha(paramAlpha);
	}
	else if (dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)){
		static_cast<generalGammaDistributionPlusInvariant*>(dist)->setAlpha(paramAlpha);
	}
	else if (dynamic_cast<generalGammaDistributionFixedCategories*>(dist)){
		static_cast<generalGammaDistributionFixedCategories*>(dist)->setAlpha(paramAlpha);
	}
	else if (dynamic_cast<gammaDistributionFixedCategories*>(dist)){
		static_cast<gammaDistributionFixedCategories*>(dist)->setAlpha(paramAlpha);
	}
	else if (dynamic_cast<gammaDistribution*>(dist)){
		static_cast<gammaDistribution*>(dist)->setAlpha(paramAlpha);
	}
	else if (dynamic_cast<generalGammaDistribution*>(dist)){
		static_cast<generalGammaDistribution*>(dist)->setAlpha(paramAlpha);
	}
	else{
		errorMsg::reportError("unknown type in distributionType");
	}	

}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE getRateBeta(distribution* dist)
{
	MDOUBLE res;
	//switch (gainLossOptions::_rateDistributionType)
	//{
	//case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):
	//	res = static_cast<generalGammaDistributionPlusInvariant*>(dist)->getBeta();
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
	//	res = static_cast<generalGammaDistributionFixedCategories*>(dist)->getBeta();
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA):
	//	res = static_cast<generalGammaDistribution*>(dist)->getBeta();
	//	break;
	//default:
	//	errorMsg::reportError("unknown type in gainLossOptions::getDistributionType");
	//}
	if(dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)){
		res = static_cast<generalGammaDistributionPlusInvariant*>(dist)->getBeta();
	}
	else if(dynamic_cast<generalGammaDistributionFixedCategories*>(dist)){
		res = static_cast<generalGammaDistributionFixedCategories*>(dist)->getBeta();
	}
	else if (dynamic_cast<generalGammaDistribution*>(dist)){
		res = static_cast<generalGammaDistribution*>(dist)->getBeta();
	}
	else{
		errorMsg::reportError("unknown type in gainLossOptions::getDistributionType");
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
void setRateBeta(distribution* dist, MDOUBLE paramBeta)
{	
	//switch (gainLossOptions::_rateDistributionType)
	//{
	//case (gainLossOptions::GENERAL_GAMMA_PLUS_INV):
	//	static_cast<generalGammaDistributionPlusInvariant*>(dist)->setBeta(paramBeta);
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA_FIXED_CATEGORIES):
	//	static_cast<generalGammaDistributionFixedCategories*>(dist)->setBeta(paramBeta);
	//	break;
	//case (gainLossOptions::GENERAL_GAMMA):
	//	static_cast<generalGammaDistribution*>(dist)->setBeta(paramBeta);
	//	break;
	//default:
	//	errorMsg::reportError("unknown type in distributionType");
	//}

	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)){
		static_cast<generalGammaDistributionPlusInvariant*>(dist)->setBeta(paramBeta);
	}
	else if (dynamic_cast<generalGammaDistributionFixedCategories*>(dist)){
		static_cast<generalGammaDistributionFixedCategories*>(dist)->setBeta(paramBeta);
	}
	else if (dynamic_cast<generalGammaDistribution*>(dist)){
		static_cast<generalGammaDistribution*>(dist)->setBeta(paramBeta);
	}
	else{
		errorMsg::reportError("unknown type in distributionType");
	}

}
/********************************************************************************************
*********************************************************************************************/
bool isAlphaOptimization(distribution* dist)
{
	if ((dynamic_cast<gammaDistribution*>(dist)) ||
		(dynamic_cast<generalGammaDistribution*>(dist)) || 
		(dynamic_cast<gammaDistributionFixedCategories*>(dist)) || 
		(dynamic_cast<gammaDistributionPlusInvariant*>(dist)) ||
		(dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)) ||
		(dynamic_cast<generalGammaDistributionFixedCategories*>(dist)) )
		return true;
	else
		return false;
}
/********************************************************************************************
*********************************************************************************************/
bool isBetaOptimization(distribution* dist)
{
	if( ((dynamic_cast<generalGammaDistribution*>(dist)) || 
		 (dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)) ||
		 (dynamic_cast<generalGammaDistributionFixedCategories*>(dist)) ) 
		 &&
		 !( (dynamic_cast<gammaDistributionFixedCategories*>(dist)) ||
		    (dynamic_cast<gammaDistribution*>(dist)) ||
		    (dynamic_cast<gammaDistributionPlusInvariant*>(dist))	 ) 	)
		return true;
	else
		return false;
}
/********************************************************************************************
*********************************************************************************************/
bool isMixOptimization(distribution* dist)
{
	if (dynamic_cast<mixtureDistribution*>(dist) )
		return true;
	else
		return false;
}
/********************************************************************************************
*********************************************************************************************/
bool isInvariantOptimization(distribution* dist, bool onlyForPrintVal)
{
	bool isInvariantDist = false;
	if (! Parameters::getInt("_isOptimizeInvariantCategoryProb") && !onlyForPrintVal )
		return false;

	if ( (dynamic_cast<generalGammaDistributionPlusInvariant*>(dist)) ||
		 (dynamic_cast<gammaDistributionPlusInvariant*>(dist))  )
		 isInvariantDist =true;		
	
	return isInvariantDist;
}
/********************************************************************************************
*********************************************************************************************/
bool isThetaOptimization()
{
	if (gainLossOptions::_characterFreqEval==gainLossOptions::optimizeOverTree && !gainLossOptions::_isRootFreqEQstationary)
		return true;
	else
		return false;
}

/********************************************************************************************
*********************************************************************************************/
void printHelp(){
	cout <<"+-------------------------------------------+"<<endl;
	cout<<"+  The gainLoss project:										"<<endl;
	cout<<"+  Analysis of Phyletic Patterns in a Likelihood Framework		"<<endl;
	cout<<"+  "<<PROG_INFO<<"						"<<endl;
	cout<<"+  Ofir Cohen - ofircohe@tau.ac.il								"<<endl;
	cout <<"use a parameter file with these options:    "<<endl;
	cout <<"+-------------------------------------------+"<<endl;
	cout <<"_seqFile									"<<endl;
	cout <<"|-------------------------------------------|"<<endl;
	cout <<"_treeFile									"<<endl;
	cout <<"|------------------------------------------|"<<endl;
	cout <<"_rootAt										"<<endl;
	cout <<"_logFile									"<<endl;
	cout <<"_logValue									"<<endl;
	cout <<"_referenceSeq								"<<endl;
	cout <<"_outDir										"<<endl;
	cout <<"_outFile									"<<endl;
	cout <<"_treeOutFile								"<<endl;
	cout <<"_numberOfGainCategories						"<<endl;
	cout <<"_numberOfLossCategories						"<<endl;
	cout <<"_numberOfRateCategories						"<<endl;
	cout <<"_maxNumOfIterationsModel					"<<endl;
	cout <<"_epsilonOptimizationModel					"<<endl;
	cout <<"_maxNumOfIterationsBBL						"<<endl;
	cout <<"_epsilonOptimizationBBL						"<<endl;
	cout <<"_epsilonOptimization						"<<endl;
	cout <<"|------------------------------------------|"<<endl;
	cout <<"_gainLossDist								"<<endl;
	cout <<"_calculateRate4site							"<<endl;
	cout <<"_calculeGainLoss4site						"<<endl;
	cout <<"_printTree									"<<endl;
	cout <<"_printPij_t									"<<endl;
	cout <<"_printLofPos								"<<endl;
	cout <<"_performOptimizations						"<<endl;
	cout <<"_isHGT_normal_Pij							"<<endl;
	cout <<"_isReversible								"<<endl;
	cout <<"...(a very partial list)					"<<endl;
	cout <<"+------------------------------------------+"<<endl;

}

/********************************************************************************************
*********************************************************************************************/
void updateGainAlpha(MDOUBLE param,
					 vector<vector<stochasticProcess*> >& spVVec,
					distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{
	bool isReversible = spVVec[0][0]->isReversible();
	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(gainDist))
		static_cast<generalGammaDistributionPlusInvariant*>(gainDist)->setAlpha(param);
	else
		static_cast<generalGammaDistribution*>(gainDist)->setAlpha(param);	

	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
		static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel())->setMu1(gainDist->rates(gainIndex),isReversible);
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}
/********************************************************************************************
*********************************************************************************************/
void updateGainBeta(MDOUBLE param, 
					vector<vector<stochasticProcess*> >& spVVec,
					distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{
	bool isReversible = spVVec[0][0]->isReversible();
	MDOUBLE normFactor;
	
	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(gainDist))
		static_cast<generalGammaDistributionPlusInvariant*>(gainDist)->setBeta(param);
	else
		static_cast<generalGammaDistribution*>(gainDist)->setBeta(param);

	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());		
		static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel())->setMu1(gainDist->rates(gainIndex),isReversible);
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normFactor = normalizeQ(spVVec, gainDist, lossDist);

}
/********************************************************************************************
*********************************************************************************************/
void updateGainProbInvariant(MDOUBLE param, distribution* gainDist)
{
	static_cast<generalGammaDistributionPlusInvariant*>(gainDist)->setInvProb(param);
}


/********************************************************************************************
*********************************************************************************************/
void updateLossAlpha(MDOUBLE param, 
									vector<vector<stochasticProcess*> >& spVVec,
									distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{

	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(lossDist))
		static_cast<generalGammaDistributionPlusInvariant*>(lossDist)->setAlpha(param);
	else
		static_cast<generalGammaDistribution*>(lossDist)->setAlpha(param);
	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());		
		static_cast<gainLossModelNonReversible*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel())->setMu2(lossDist->rates(lossIndex));
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}
/********************************************************************************************
*********************************************************************************************/
void updateLossBeta(MDOUBLE param, 
								   vector<vector<stochasticProcess*> >& spVVec,
								   distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{

	if (dynamic_cast<generalGammaDistributionPlusInvariant*>(gainDist))
		static_cast<generalGammaDistributionPlusInvariant*>(lossDist)->setBeta(param);
	else
		static_cast<generalGammaDistribution*>(lossDist)->setBeta(param);	
	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());		
		static_cast<gainLossModelNonReversible*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel())->setMu2(lossDist->rates(lossIndex));
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}
/********************************************************************************************
*********************************************************************************************/
void updateLossProbInvariant(MDOUBLE param, distribution* lossDist)
{
	static_cast<generalGammaDistributionPlusInvariant*>(lossDist)->setInvProb(param);
}


/********************************************************************************************
*********************************************************************************************/
void updateRateAlpha(MDOUBLE param, 
					vector<vector<stochasticProcess*> >& spVVec,
					distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{
	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
		setRateAlpha(spVVec[gainIndex][lossIndex]->distr(), param);
		//static_cast<gammaDistribution*>(spVVec[gainIndex][lossIndex]->distr())->setAlpha(param);
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}
/********************************************************************************************
*********************************************************************************************/
void updateRateProbInvariant(MDOUBLE param, 
					 vector<vector<stochasticProcess*> >& spVVec,
					 distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{
	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
		static_cast<generalGammaDistributionPlusInvariant*>(spVVec[gainIndex][lossIndex]->distr())->setInvProb(param);
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}
/********************************************************************************************
*********************************************************************************************/
void updateTheta(MDOUBLE param, 
				vector<vector<stochasticProcess*> >& spVVec,
				distribution * gainDist, distribution * lossDist, bool isNormalizeQ)
{
	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());		
		(static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel()))->setTheta(param);
	}
	if(gainLossOptions::_isNormalizeQinSpVVec && isNormalizeQ)
		normalizeQ(spVVec, gainDist, lossDist);
}

/********************************************************************************************
*********************************************************************************************/
void cloneSpVVec(vector<vector<stochasticProcess*> >& spVVec, vector<vector<stochasticProcess*> >& neWspVVec){
	
	neWspVVec.resize(spVVec.size());
	for (int gainCategor=0; gainCategor<spVVec.size(); gainCategor++){
		neWspVVec[gainCategor].resize(spVVec[0].size());
		for (int lossCategor=0; lossCategor<spVVec[0].size(); lossCategor++){		
			neWspVVec[gainCategor][lossCategor] = spVVec[gainCategor][lossCategor]->clone();
		}
	}
}

/********************************************************************************************
*********************************************************************************************/
void deleteSpVVec(vector<vector<stochasticProcess*> >* spVVec_p){

	if(spVVec_p){
		for (int gainCategor=0; gainCategor<spVVec_p->size(); gainCategor++){
			for (int lossCategor=0; lossCategor<(*spVVec_p)[0].size(); lossCategor++){		
				delete (*spVVec_p)[gainCategor][lossCategor];
			}
		}
	}
}


/********************************************************************************************
*********************************************************************************************/
void clearVVVV(VVVVdouble& vetor){	
	for (int i=0;i<vetor.size();++i){
		if(vetor.size()==0)
			break;
		for (int j=0;j<vetor[i].size();++j){
			if(vetor[i].size()==0)
				break;
			for (int k=0;j<vetor[i][j].size();++k){
				if(vetor[i][j].size()==0)
					break;
				if(vetor[i][j][k].size()==0)
					break;
				vetor[i][j][k].clear();
			}
			vetor[i][j].clear();
		}
		vetor[i].clear();
	}
	vetor.clear();
}

/********************************************************************************************
*********************************************************************************************/
void clearVVV(VVVdouble& vetor){	
	for (int i=0;i<vetor.size();++i){
		if(vetor.size()==0)
			break;
		for (int j=0;j<vetor[i].size();++j){
			if(vetor[i].size()==0)
				break;
			vetor[i][j].clear();
		}
		vetor[i].clear();
	}
	vetor.clear();
}


/********************************************************************************************
*********************************************************************************************/
void resizeVVVV(int dim1, int dim2, int dim3, int dim4,  VVVVdouble& vetor){
	
	vetor.resize(dim1);
	for (int posNum=0;posNum<vetor.size();++posNum){
		vetor[posNum].resize(dim2);
		for (int n=0;n<vetor[posNum].size();++n){
			resizeMatrix(vetor[posNum][n],dim3,dim4);
		}
	}
}
/********************************************************************************************
*********************************************************************************************/
void resizeVVV(int dim1, int dim2, int dim3, VVVdouble& vetor){	
	vetor.resize(dim1);
	for (int n=0;n<vetor.size();++n){
		resizeMatrix(vetor[n],dim2,dim3);
	}
}

///********************************************************************************************
//*********************************************************************************************/
//MDOUBLE getDistance2ROOT(const tree::nodeP &myNode){	
//	if(myNode->isRoot())
//		return 0.0;
//	else
//		return ( myNode->dis2father() + getDistance2ROOT(myNode->father()) );
//}
///********************************************************************************************
//getMinimalDistance2OTU()
//This implementation is only for binary trees.
//Can easily be generalized to arbitrary number of sons.
//*********************************************************************************************/
//MDOUBLE getMinimalDistance2OTU(const tree::nodeP &myNode){	
//	if(myNode->isLeaf())
//		return 0.0;
//	else{
//		if(myNode->getNumberOfSons()>2)
//			LOGnOUT(3,<<" ERROR: getMinimalDistance2OTU is only for binary trees, and this node "
//			<<myNode->name()<<" is with "<<myNode->getNumberOfSons()<<"sons.\n The return value is only for first 2 sons\n");
//
//		return ( min(
//			myNode->getSon(0)->dis2father() + getMinimalDistance2OTU(myNode->getSon(0)),
//			myNode->getSon(1)->dis2father() + getMinimalDistance2OTU(myNode->getSon(1))
//			) );
//
//	}
//}			



/********************************************************************************************
*********************************************************************************************/
void fillVnames(Vstring& Vnames,const tree& tr){
	vector<tree::nodeP> vAllNodes;
	tr.getAllNodes(vAllNodes,tr.getRoot());
	Vnames.resize(vAllNodes.size());
	for (int i = 0; i<vAllNodes.size();++i)
		Vnames[vAllNodes[i]->id()] = vAllNodes[i]->name();
}			
/********************************************************************************************
*********************************************************************************************/
void P11forgain(ostream& out)   {
	string P11forgain = gainLossOptions::_outDir + "//" + "P11forgain.txt"; 
	ofstream P11forgainStream(P11forgain.c_str());
	P11forgainStream.precision(PRECISION);
	
	MDOUBLE loss = 0.0;
	MDOUBLE dist = 0.3;
	MDOUBLE increment = 0.1;
	P11forgainStream <<"gain"<<"\t"<<"loss"<<"\t"<<"dist"<<"\t"<<"P11"<<endl;
	for (int ind = 1; ind<10000; ind++){
		MDOUBLE gain = ind*increment;		
		MDOUBLE eigenvalue =  -(gain + loss);
		MDOUBLE P11 = gain/(-eigenvalue) + exp(eigenvalue*dist)*(1 - gain/(-eigenvalue));
		P11forgainStream <<gain<<"\t"<<loss<<"\t"<<dist<<"\t"<<P11<<endl;
	}
}


/********************************************************************************************
normalize the Q matrix so average rate of substitution = 1
*********************************************************************************************/
MDOUBLE normalizeQ(vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist){
	MDOUBLE sumPijQij=0.0;
	MDOUBLE scale;

	//int numOfSPs = gainDist->categories()*lossDist->categories();
	//for (int i=0; i < numOfSPs; ++i) {
	//	int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
	//	int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
	//	sumPijQij+=gainDist->ratesProb(gainIndex)*lossDist->ratesProb(lossIndex)
	//		*(static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel()))->sumPijQij();	
	//}	
	//if (sumPijQij ==0){
	//	errorMsg::reportError("Error in normalizeMatrices - sumPijQij=0");
	//}
	sumPijQij = sumPijQijVec(spVVec,  gainDist,  lossDist);
	scale = (1.0 / sumPijQij);
	normVec(scale, spVVec,  gainDist,  lossDist);
	//for (int i=0; i < numOfSPs; ++i) {
	//	int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
	//	int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
	//	(static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel()))->norm(scale);	
	//}
	////MDOUBLE AlphaGainLossRatio = getRateAlpha(gainDist)/getRateAlpha(lossDist);
	//MDOUBLE newGainBeta = getRateBeta(gainDist)/scale;
	//updateGainBeta(newGainBeta,spVVec,gainDist,lossDist,false);	// BUG fixed. If only Q matrices are corrected -> problem
	//MDOUBLE newLossBeta = getRateBeta(lossDist)/scale;
	//updateLossBeta(newLossBeta,spVVec,gainDist,lossDist,false);
	return sumPijQij;
}


/********************************************************************************************
normalize the Q matrix so average rate of substitution = 1
*********************************************************************************************/
MDOUBLE sumPijQijVec(vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist){
	MDOUBLE sumPijQij=0.0;
	MDOUBLE scale;

	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
		sumPijQij+=gainDist->ratesProb(gainIndex)*lossDist->ratesProb(lossIndex)
			*(static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel()))->sumPijQij();	
	}	
	if (sumPijQij ==0){
		errorMsg::reportError("Error in normalizeMatrices - sumPijQij=0");
	}
	return sumPijQij;
}


/********************************************************************************************
normalize the Q matrix so average rate of substitution = 1
*********************************************************************************************/
void normVec(const MDOUBLE scale, vector<vector<stochasticProcess*> >& spVVec, distribution * gainDist, distribution * lossDist){

	int numOfSPs = gainDist->categories()*lossDist->categories();
	for (int i=0; i < numOfSPs; ++i) {
		int gainIndex =fromIndex2gainIndex(i,gainDist->categories(),lossDist->categories());
		int lossIndex =fromIndex2lossIndex(i,gainDist->categories(),lossDist->categories());
		(static_cast<gainLossModel*>(spVVec[gainIndex][lossIndex]->getPijAccelerator()->getReplacementModel()))->norm(scale);	
	}
	
	MDOUBLE newGainBeta = getRateBeta(gainDist)/scale;
	updateGainBeta(newGainBeta,spVVec,gainDist,lossDist,false);	// BUG fixed. If only Q matrices are corrected -> problem
	MDOUBLE newLossBeta = getRateBeta(lossDist)/scale;
	updateLossBeta(newLossBeta,spVVec,gainDist,lossDist,false);	
}

/********************************************************************************************/
MDOUBLE normalizeQ(stochasticProcess*  sp){	
	MDOUBLE sumPijQij=(static_cast<gainLossModel*>(sp->getPijAccelerator()->getReplacementModel()))->sumPijQij();	
	(static_cast<gainLossModel*>(sp->getPijAccelerator()->getReplacementModel()))->norm( 1.0/sumPijQij );	
	return sumPijQij;
}


/********************************************************************************************
*********************************************************************************************/
MDOUBLE computeExpectationOfStationaryFrequency(distribution* gainDist, distribution* lossDist){
	MDOUBLE estimatedStationaryFreq=0;
	//if(gainDist->categories() == lossDist->categories()){
		for(int i=0; i<gainDist->categories(); ++i){
			for(int j=0; j<lossDist->categories(); ++j){
			//if(gainDist->ratesProb(i) == lossDist->ratesProb(i)){
				estimatedStationaryFreq += (gainDist->rates(i)/(gainDist->rates(i)+lossDist->rates(j)))*  gainDist->ratesProb(i)*lossDist->ratesProb(j);
			//}
			//else{
			//	LOGnOUT(4,<<" WARN: computeExpectationOfStationaryFrequency did not compute Theta" <<endl);
			//}
			}
		}
	//}
	//else{
	//	LOGnOUT(4,<<" WARN: computeExpectationOfStationaryFrequency did not compute Theta" <<endl);
	//}
	if(estimatedStationaryFreq<0 || estimatedStationaryFreq>1){
		LOGnOUT(4,<<" ERROR: computeExpectationOfStationaryFrequency <0 or >1" <<estimatedStationaryFreq<<endl);
		return 0;
	}
	return estimatedStationaryFreq;
}


/********************************************************************************************
exp(gain/loss) (not exp(gain)/exp(loss) )
Each gain/loss ratio is weighted by the rateCategories probability 
(their duplication is the matrix probability)
*********************************************************************************************/
MDOUBLE computeExpectationOfGainLossRatio(distribution* gainDist, distribution* lossDist){

	MDOUBLE compGainLossRatio=0;

	for(int i=0; i<gainDist->categories(); ++i){
		for(int j=0; j<lossDist->categories(); ++j){				
			compGainLossRatio += gainDist->rates(i)/lossDist->rates(j) *gainDist->ratesProb(i)*lossDist->ratesProb(j);
		}
	}
	if(compGainLossRatio<0 ){
		LOGnOUT(4,<<" ERROR: compGainLossRatio <0 " <<compGainLossRatio<<endl);
		return 0;
	}
	return compGainLossRatio;
}


/********************************************************************************************
exp(gain)/exp(loss) (not exp(gain/loss)  )
Each gain/loss ratio is weighted by the rateCategories probability 
(their duplication is the matrix probability)
*********************************************************************************************/
MDOUBLE computeExpOfGainByExpOfLossRatio(distribution* gainDist, distribution* lossDist){
	
	MDOUBLE compGainLossRatio = 1;
	MDOUBLE ExpGain=0;
	MDOUBLE ExpLoss=0;
	
	//for(int i=0; i<gainDist->categories(); ++i){
	//	ExpGain += gainDist->rates(i) *gainDist->ratesProb(i);
	//}
	ExpGain = rateExpectation(gainDist);
	//for(int j=0; j<lossDist->categories(); ++j){				
	//	ExpLoss += lossDist->rates(j) *lossDist->ratesProb(j);		
	//}
	ExpLoss = rateExpectation(lossDist);
	compGainLossRatio = ExpGain/ExpLoss;
	if(compGainLossRatio<0 ){
		LOGnOUT(4,<<" ERROR: compGainLossRatio <0 " <<compGainLossRatio<<endl);
		return 0;
	}
	return compGainLossRatio;
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE rateExpectation(distribution* dist){
	MDOUBLE ExpRate=0;
	bool isWithInvariant = isInvariantOptimization(dist);
	if(isWithInvariant){
		for(int i=0; i<dist->categories(); ++i){
			ExpRate += dist->rates(i) *dist->ratesProb(i);
		}
	}else{
		ExpRate = getRateAlpha(dist)/getRateBeta(dist);
	}
	return ExpRate;
}



/********************************************************************************************
Mixture
*********************************************************************************************/
void printMixtureParams(stochasticProcess*  sp) 
{	
	mixtureDistribution * pMixture = static_cast<mixtureDistribution*>(sp->distr());
	for (int k = 0; k < pMixture->getComponentsNum(); ++k)
	{
		LOGnOUT(4, << "comp="<<k<<" Alp/Beta= "<<pMixture->getAlpha(k)/pMixture->getBeta(k)<<" alpha= "<<pMixture->getAlpha(k) << " beta= " <<pMixture->getBeta(k)<<" Prob= "<<pMixture->getComponentProb(k)<<endl);  
	}
}

/********************************************************************************************
*********************************************************************************************/
stochasticProcess* startStochasticProcessSimpleGamma(MDOUBLE init_gain, MDOUBLE init_loss, Vdouble& freq, int numberOfRateCategories)
{ 
	LOGnOUT(4,<<"startStochasticProcess SimpleGamma of "<<numberOfRateCategories<<" categories... \nwith gain="
		<<init_gain<<" loss="<<init_loss<<" root(1)="<<freq[1]<<endl);
	stochasticProcess *spSimple;
	replacementModel* glm;
	if(!gainLossOptions::_isReversible){
		glm = new gainLossModelNonReversible(init_gain,init_loss,freq,gainLossOptions::_isRootFreqEQstationary,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}
	else{
		glm = new gainLossModel(init_gain, freq,gainLossOptions::_isRootFreqEQstationary, true,gainLossOptions::_isHGT_normal_Pij,gainLossOptions::_isHGT_with_Q);
	}	
	trivialAccelerator* pijAcc = new trivialAccelerator(glm);
	MDOUBLE initAlphaRate = gainLossOptions::_userAlphaRate;	
	distribution* rateDist = new gammaDistribution(initAlphaRate,numberOfRateCategories);
	spSimple = new stochasticProcess(rateDist,pijAcc,gainLossOptions::_isReversible);	
	if (rateDist) delete rateDist;	//at r4s after the sp object is created all other objects dynamically constructed are deleted
	if (pijAcc) delete pijAcc;
	if (glm) delete glm;
	if(gainLossOptions::_isNormalizeQ)
		normalizeQ(spSimple);
	return spSimple;
}

/********************************************************************************************
Assume the first site number =1 in selectedSites file
*********************************************************************************************/
void readIntegersFromFileIntoVector(Vint& intVector, const int maxAllowed, const int minAllowed, string* inFile,Vint* evolvingSites){
	if(inFile){
		ifstream myReadFile;
		myReadFile.open(inFile->c_str());
		
		bool isSiteLegal = true;
		vector<int>::iterator vec_iter;
		if (myReadFile.is_open()) {
			while (!myReadFile.eof()) {
				int site = -1; // thus, only if read a new site
				myReadFile>>site;
				--site; // since count starts from 1 denoted sites
				isSiteLegal = true;
				if(evolvingSites){
					vec_iter = evolvingSites->begin();
					while(vec_iter < evolvingSites->end() && *vec_iter != site)
						vec_iter++;
					if (vec_iter==evolvingSites->end() )
						isSiteLegal = false;
				}
				if(isSiteLegal && site <= maxAllowed && site>=minAllowed){
					intVector.push_back(site);
				}
				else
					LOGnOUT(4,<<" WARN selectedSitesForCorrelation - "<<site<<" is not in seq length or not legal (not among evolvingSites). Thus not included"<<maxAllowed<<endl);
			}
		}
		else
			LOGnOUT(4,<<" Error selectedSitesForCorrelation file "<<inFile<<" can't be opened"<<endl);
		myReadFile.close();
	}
	else{
		for(int site=minAllowed; site<maxAllowed; ++site){
			intVector.push_back(site);
		}
	}
}


/********************************************************************************************
Flat Tree BeforeOpt
*********************************************************************************************/
void FlatTree(tree& trForSM , MDOUBLE defaultBranchLength){
	LOGnOUT(4,<<"\nNote: FlatTreeBeforeOpt..  with defaultBranchLength="<<defaultBranchLength<<endl);
	treeIterDownTopConst tIt(trForSM);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		mynode->setDisToFather(defaultBranchLength);
	}
}


/********************************************************************************************
fill map_PosXY
*********************************************************************************************/
void computeRateValPerPos(VVVVdouble& expChanges_PosNodeXY, VVVdouble& map_PosXY){
	int numOfPositions = expChanges_PosNodeXY.size();
	int numOfBranches = expChanges_PosNodeXY[0].size();
	int AlphSize = expChanges_PosNodeXY[0][0].size(); // =2
	resizeVVV(numOfPositions,AlphSize,AlphSize,map_PosXY);

	for (int pos = 0; pos <numOfPositions; ++pos){
		for(int b=0;b<numOfBranches;++b){
			for(int j=0;j<AlphSize;++j){
				for(int k=0;k<AlphSize;++k){
					if(gainLossOptions::_isNminBasedOnCountBranchesOverCutOff && expChanges_PosNodeXY[pos][b][j][k]>gainLossOptions::_probCutOffCounts)
						map_PosXY[pos][j][k] += 1;
					else if(!gainLossOptions::_isNminBasedOnCountBranchesOverCutOff)
						map_PosXY[pos][j][k] += expChanges_PosNodeXY[pos][b][j][k];
				}
			}
		}		
			
	}

}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE computeNminRforCorrelWithGainAndLoss(MDOUBLE gainVal, MDOUBLE lossVal){ 
	return (gainVal+lossVal)/2.0; 
}
