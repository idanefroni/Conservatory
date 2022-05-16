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


#ifndef ___GAIN_LOSS_ALPH
#define ___GAIN_LOSS_ALPH

#include "alphabet.h"
#include "errorMsg.h"

class gainLossAlphabet : public alphabet {
public:
	explicit gainLossAlphabet(); 
	virtual ~gainLossAlphabet() {}
	virtual alphabet* clone() const { return new gainLossAlphabet(*this); }
	int unknown() const  {return -2;}
	int gap() const  {errorMsg::reportError("The method indel::gap() is used"); return -1;} // What is it for ? I don't need this !!!
	int size() const {return 2;} // presence or absence only
	int stringSize() const {return 1;} // one letter code.
	int relations(const int charInSeq, const int charToCheck) const;
	int fromChar(const string& str, const int pos) const;
	int fromChar(const char s) const;
	string fromInt(const int in_id) const;
	vector<int> fromString(const string& str) const;
	bool isSpecific(const int id) const {return (id>=0 && id < size());}

};

#endif
