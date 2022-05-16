#ifndef ___GAP__
#define ___GAP__


#include "definitions.h"

#include <iostream>
using namespace std;



class gaps {
public:

	explicit gaps()	{};
	~gaps() {
		for(int i=0; i<_gaps.size();++i){
			gap* gap2delete = _gaps[i];
			delete gap2delete;
		}
	};
	////////////////////////////////////////////////////////////////////////// inner class
	class gap {
	public:

		explicit gap(int coord_5p, int coord_3p, int seqID, int coord_5Abs):
			_coord_5p(coord_5p), _coord_3p(coord_3p), _seqID(seqID),_coord_5Abs(coord_5Abs)	{};
		~gap() {};
		int getCoord5() const {return _coord_5p;};
		int getCoord3()const {return _coord_3p;};
		int getSeqID() const {return _seqID;};
		int getCoord5Abs() const {return _coord_5Abs;};
		int getLength() const {return _coord_3p-_coord_5p+1;};

	private:
		int _coord_5p; 
		int _coord_3p;
		int _seqID;
		int _coord_5Abs; 
	};
	////////////////////////////////////////////////////////////////////////// end	
	
	
	gap* operator[](const int i) {return  _gaps[i];} // get the ID of the gap. Return the gap itself.
	
	int numOfGaps(){return _gaps.size();}
	
	/********************************************************************************************
	insertNewGap
	// Sort the vector containing all indels by  I =(i1,i2), K =(k1,k2), I<K iff i1<k1 or i1=k1 and i2<k2
	*********************************************************************************************/
	void insertNewGap(int coord_5p, int coord_3p, int seqID, int coord_5Abs)
	{
		gap* gap_p = new gap(coord_5p, coord_3p,seqID, coord_5Abs);
		//_gaps.push_back(gap_p);

		vector<gap*>::iterator iter;
		int position = 0;
		iter = _gaps.begin();
		while( iter!=_gaps.end() &&
			  ( (*iter)->getCoord5() < coord_5p || 
			    ((*iter)->getCoord5() <= coord_5p &&  (*iter)->getCoord3() < coord_3p) ) )
		{
			iter++;
			position++;
		}
		_gaps.insert(iter, gap_p);
	};

	//////////////////////////////////////////////////////////////////////////
	void insertNewGap(gap* gap_p){
		vector<gap*>::iterator iter;
		int position = 0;
		iter = _gaps.begin();
		while( iter!=_gaps.end() &&
			( (*iter)->getCoord5() < gap_p->getCoord5() || 
			((*iter)->getCoord5() <= gap_p->getCoord5() &&  (*iter)->getCoord3() < gap_p->getCoord3()) ) )
		{
			iter++;
			position++;
		}
		_gaps.insert(iter, gap_p);
	};
	
	void insertNewGapNotSorted(gap* gap_p){
		_gaps.push_back(gap_p);
	};

	//////////////////////////////////////////////////////////////////////////
	void printGaps(){
		vector<gap*>::iterator iter;
		iter = _gaps.begin();
		while( iter!=_gaps.end())
		{
			cout<<"Gap  "<<(*iter)->getCoord5()<<" "<<(*iter)->getCoord3()<<endl;
			iter++;
		}
	};

private:
	vector<gap*> _gaps;


};


#endif
