// $Id: suffStatComponent.h 9253 2011-01-31 01:37:21Z rubi $

#ifndef ___SUFF_STAT_COMPONENT
#define ___SUFF_STAT_COMPONENT

#include "definitions.h"
#include <vector>
using namespace std;

// spec = for a specific node. global = for all the nodes
// hom = no rate variation. gam = with rate variation. gamProportional = with gobal and local rate variation
// pos = for one position
//-------------------------------------------------------------
class suffStatSpecHomPos{ // this is for a specific node.
	public:
		void set(const int letter,const doubleRep& val) {
			_V[letter]=val;
		}
		
		doubleRep get(const int letter) const	{
		    doubleRep tmp=_V[letter];
//		    cout << "tmp =";
//		    tmp.outputn(cout);
		    
		    return tmp;
		}

		void allocatePlace(const int alphabetSize) {
			_V.resize(alphabetSize);
		}
		bool isEmpty (){return (_V.empty());};
		int size() const {return _V.size();}

	private:
		vector<doubleRep> _V;//size = letter
};
//-------------------------------------------------------------
/*
class suffStatSpecGamPos{// this is for a specific node with rates
	public:
		void set(const int rateCategor,
			const int letter,const MDOUBLE val) {
			_V[rateCategor].set(letter,val);
		}
		
		MDOUBLE get(const int rateCategor,
			const int letter) const	{
			return _V[rateCategor].get(letter);
		}
		void allocatePlace(const int numberOfrateCategories,const int alphabetSize) {
			_V.resize(numberOfrateCategories);
			for (int i=0; i < numberOfrateCategories; ++i) {
				_V[i].allocatePlace(alphabetSize);
			}
		}
		bool isEmpty (){return (_V.empty());};
	private:
		vector<suffStatSpecHomPos> _V;//rateCategor,letter
};
*/
//-------------------------------------------------------------
/*
class suffStatSpecGam{// this is for a specific node with rates
	public:
		void set(const int pos,const int rateCategor,
			const int letter,const MDOUBLE val) {
			_V[pos].set(rateCategor,letter,val);
		}
		
		MDOUBLE get(const int pos,const int rateCategor,
			const int letter) const	{
			return _V[pos].get(rateCategor,letter);
		}

		void allocatePlace(const int pos,const int numberOfrateCategories,const int alphabetSize) {
			_V.resize(pos);
			for (int i=0;i<pos;++i) _V[i].allocatePlace(numberOfrateCategories,alphabetSize);
		}
		bool isEmpty (){return (_V.empty());};
		suffStatSpecGamPos& operator[] (int index) {return _V[index];}
		const suffStatSpecGamPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatSpecGamPos> _V;//pos,rateCategor,letter
};
*/
//-------------------------------------------------------------
/*
class suffStatGlobalGam {
public:
	MDOUBLE get(const int nodeId, const int pos,const int rateCategor,
			const int letter) const	{
		return _V[nodeId].get(pos,rateCategor,letter);
	}
	void allocatePlace(const int numOfNodes,
						const int pos,
						const int numberOfrateCategories,
						const int alphabetSize) {
		_V.resize(numOfNodes);
		for (int i=0;i<numOfNodes;++i) _V[i].allocatePlace(pos,numberOfrateCategories,alphabetSize);
	}
	int size() const {return _V.size();}
	suffStatSpecGam& operator[] (int index) {return _V[index];}
	const suffStatSpecGam& operator[] (int index) const {return _V[index];}

private:
	vector<suffStatSpecGam> _V;
};
*/
//-------------------------------------------------------------
class suffStatGlobalHomPos{ // this is for all nodes
	public:
		void set(const int nodeId,const int letter,const doubleRep val) {
			_V[nodeId].set(letter,val);
		}
		
		doubleRep get(const int nodeId,const int letter) const	{
		    doubleRep tmp(_V[nodeId].get(letter));
//		    tmp;
		    
//		    cout << "tmp2=";
//		    tmp.outputn(cout);
		    return tmp;
		}

		void allocatePlace(const int numOnNodes,const int alphabetSize) {
			_V.resize(numOnNodes);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());};
		int size() const {return _V.size();}
	private:
		vector<suffStatSpecHomPos> _V;//size = number of nodes.
};
//-------------------------------------------------------------
class suffStatGlobalGamPos{ // this is for all nodes
	public:
		void set(const int categor,const int nodeId,const int letter,const doubleRep val) {
			_V[categor].set(nodeId,letter,val);
		}
		
		doubleRep get(const int categor,const int nodeId,const int letter) const	{
			return _V[categor].get(nodeId,letter);
		}

		void allocatePlace(const int categor,const int numOnNodes,const int alphabetSize) {
			_V.resize(categor);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(numOnNodes,alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());}
		int size() const {return _V.size();}

	suffStatGlobalHomPos& operator[] (int index) {return _V[index];}
	const suffStatGlobalHomPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatGlobalHomPos> _V;//size = number of categories
};
//-------------------------------------------------------------
class suffStatGlobalGamProportionalPos{ // this is for all nodes
	public:
		void set(const int globalRateCategor,const int localRateCategor,const int nodeId,const int letter,const doubleRep val) {
			_V[globalRateCategor].set(localRateCategor,nodeId,letter,val);
		}
		
		doubleRep get(const int globalRateCategor,const int localRateCategor,const int nodeId,const int letter) const	{
			return _V[globalRateCategor].get(localRateCategor,nodeId,letter);
		}

		void allocatePlace(const int globalRateCategor,const int localRateCategor,const int numOnNodes,const int alphabetSize) {
			_V.resize(globalRateCategor);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(localRateCategor,numOnNodes,alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());}
		int size() const {return _V.size();}

	suffStatGlobalGamPos& operator[] (int index) {return _V[index];}
	const suffStatGlobalGamPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatGlobalGamPos> _V;//size = number of global rate categories
};
//-------------------------------------------------------------
class suffStatGlobalGam{ // this is for all positions (and for all nodes).
	public:
		void set(const int pos,const int categor,const int nodeId,const int letter,const doubleRep val) {
			_V[pos].set(categor,nodeId,letter,val);
		}
		
		doubleRep get(const int pos,const int categor,const int nodeId,const int letter) const	{
			return _V[pos].get(categor,nodeId,letter);
		}

		void allocatePlace(const int pos,const int categor,const int numOnNodes,const int alphabetSize) {
			_V.resize(pos);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(categor,numOnNodes,alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());}
		int size() const {return _V.size();}
	suffStatGlobalGamPos& operator[] (int index) {return _V[index];}
	const suffStatGlobalGamPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatGlobalGamPos> _V;
};

//-------------------------------------------------------------
class suffStatGlobalGamProportional{ // this is for all positions (and for all nodes).
	public:
		void set(const int pos,const int globalRateCategor,const int localRateCategor,const int nodeId,const int letter,const doubleRep val) {
			_V[pos].set(globalRateCategor,localRateCategor,nodeId,letter,val);
		}
		
		doubleRep get(const int pos,const int globalRateCategor,const int localRateCategor,const int nodeId,const int letter) const	{
			return _V[pos].get(globalRateCategor,localRateCategor,nodeId,letter);
		}

		void allocatePlace(const int pos,const int globalRateCategor,const int localRateCategor,const int numOnNodes,const int alphabetSize) {
			_V.resize(pos);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(globalRateCategor,localRateCategor,numOnNodes,alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());}
		int size() const {return _V.size();}
	suffStatGlobalGamProportionalPos& operator[] (int index) {return _V[index];}
	const suffStatGlobalGamProportionalPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatGlobalGamProportionalPos> _V;
};

// from ItayM not to use with the EM algorithm.
class suffStatGlobalHom{ // this is for all positions (and for all nodes).
	public:
		void set(const int pos, const int nodeId, const int letter,const doubleRep val) {
			_V[pos].set(nodeId, letter, val);
		}
		
		doubleRep get(const int pos, const int nodeId, const int letter) const {
			return _V[pos].get(nodeId, letter);
		}

		void allocatePlace(const int pos, const int numOnNodes, const int alphabetSize) {
			_V.resize(pos);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(numOnNodes, alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());};
	suffStatGlobalHomPos& operator[] (int index) {return _V[index];}
	const suffStatGlobalHomPos& operator[] (int index) const {return _V[index];}
	private:
		vector<suffStatGlobalHomPos> _V;
};


#endif

