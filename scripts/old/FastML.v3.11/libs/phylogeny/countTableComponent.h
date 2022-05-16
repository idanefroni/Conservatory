// $Id: countTableComponent.h 9595 2011-06-30 18:56:40Z rubi $

#ifndef ___COUNT_TABLE_COMPONENT
#define ___COUNT_TABLE_COMPONENT

#include "definitions.h"
#include <iostream>
#include <cassert>

class countTableComponentHom{
public:  

	void setCount( const int letter1,
					const int letter2,
					const MDOUBLE val) {
		_countValues[letter1][letter2]=val;
	}
	int alphabetSize() const {return _countValues.size();}
	void zero();
	MDOUBLE getCounts(	const int letter1,
						const int letter2) const	{
		return _countValues[letter1][letter2];
	}
	void addToCounts(const int let1,const int let2,const MDOUBLE val) {
		_countValues[let1][let2]+=val;
	}
	int getSize() const {return _countValues.size();}
	bool isEmpty (){return (_countValues.empty());};
	void countTableComponentAllocatePlace(const int alphabetSize);
	void printTable(ostream & out) const;
	const Vdouble& operator[] (int i) const {return _countValues[i];}
private:					
	VVdouble _countValues;//letter1,letter2

};

class countTableComponentGam{
public:  

	void setCount( const int letter1,
					const int letter2,
					const int rateCategor,
					const MDOUBLE val) {
		_countValues[rateCategor].setCount(letter1,letter2,val);
	}
	
	int alphabetSize() const {return _countValues.empty()?0:_countValues[0].alphabetSize();}
	void zero(){
		for (int rateCat=0; rateCat < _countValues.size(); ++rateCat) _countValues[rateCat].zero();
	}


	MDOUBLE getCounts(	const int letter1,
						const int letter2,
						const int rateCategor) const {
		assert(_countValues[rateCategor].getCounts(letter1,letter2)>=0);
		return _countValues[rateCategor].getCounts(letter1,letter2);
	}

	void addToCounts(const int let1,const int let2,
		const int rate,const MDOUBLE val) {
		_countValues[rate].addToCounts(let1,let2,val);
	}

	bool isEmpty (){return (_countValues.empty());};

	void countTableComponentAllocatePlace(const int alphabetSize,
		const int numberOfrateCategories) {
		_countValues.resize(numberOfrateCategories);
		for (int rateCat=0; rateCat < _countValues.size(); ++rateCat){
			_countValues[rateCat].countTableComponentAllocatePlace(alphabetSize);
		}
	}
	void printTable(ostream & out) const {
		for (int rateCat=0; rateCat < _countValues.size(); ++rateCat) {
			_countValues[rateCat].printTable(out);
		}
	}
	int getSize() const {return _countValues.size();}
	countTableComponentHom& operator[] (int i) {return _countValues[i];}
	const countTableComponentHom& operator[] (int i) const {return _countValues[i];}
private:
	vector<countTableComponentHom> _countValues;//letter1,letter2,rateCategor 

};

class countTableComponentGamProportional{
public:  

	void setCount( const int letter1,
					const int letter2,
					const int globalRateCategor,
					const int localRateCategor,
					const MDOUBLE val) {
		_countValues[globalRateCategor].setCount(letter1,letter2,localRateCategor,val);
	}
	
	int alphabetSize() const {return _countValues.empty()?0:_countValues[0].alphabetSize();}
	void zero(){
		for (int globalRateCat=0; globalRateCat < _countValues.size(); ++globalRateCat) _countValues[globalRateCat].zero();
	}


	MDOUBLE getCounts(	const int letter1,
						const int letter2,
						const int globalRateCategor,
						const int localRateCategor) const {
		assert(_countValues[globalRateCategor].getCounts(letter1,letter2,localRateCategor)>=0);
		return _countValues[globalRateCategor].getCounts(letter1,letter2,localRateCategor);
	}

	void addToCounts(const int let1,const int let2,
		const int globalRate,const int localRate,const MDOUBLE val) {
		_countValues[globalRate].addToCounts(let1,let2,localRate,val);
	}

	bool isEmpty (){return (_countValues.empty());}

	void countTableComponentAllocatePlace(const int alphabetSize,
		const int numberOfGlobalRateCategories,const int numberOfLocalRateCategories) {
		_countValues.resize(numberOfGlobalRateCategories);
		for(int globalRateCat = 0;globalRateCat < _countValues.size(); ++globalRateCat){
			_countValues[globalRateCat].countTableComponentAllocatePlace(alphabetSize,numberOfLocalRateCategories);
		}
	}
	void printTable(ostream & out) const {
		for (int globalRateCat=0; globalRateCat < _countValues.size(); ++globalRateCat) {
			_countValues[globalRateCat].printTable(out);
		}
	}
	int getSize() const {return _countValues.size();}
	countTableComponentGam& operator[] (int i) {return _countValues[i];}
	const countTableComponentGam& operator[] (int i) const {return _countValues[i];}
private:
	vector<countTableComponentGam> _countValues;//letter1,letter2,globalRateCategor,localRateCategor

};

#endif

