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
#include "gainLossAlphabet.h"

gainLossAlphabet::gainLossAlphabet() {}

int gainLossAlphabet::fromChar(const char s) const{
	switch (s) {
		case '0': return 0; break;	
		case '1': return 1; break;
		case '-' : case'_' : return -1; break;

		default:
			vector<string> err;
			err.push_back(" The gainLoss sequences contained the character: ");
			err[0]+=s;
			err.push_back(" gainLoss was not one of the following: ");
			err.push_back(" 0, 1");
			errorMsg::reportError(err);
	}// end of switch
	return -99; // never suppose to be here.	
}// end of function

vector<int> gainLossAlphabet::fromString(const string &str) const {
	vector<int> vec;
	for (int i=0;i<str.size();i++)
		vec.push_back(fromChar(str[i]));
	return vec;
}

string gainLossAlphabet::fromInt(const int in_id) const{
	char res = 0;
	switch (in_id) {
		case 0 : res = '0'  ; break;
		case 1 : res = '1'  ; break;
		case -2 : res = '-'; break;
		default:
			vector<string> err;
			err.push_back("unable to print gainLoss_id. gainLossl_id was not one of the following: ");
			err.push_back("0,1,2");
			errorMsg::reportError(err);
	}//end of switch
	string vRes;
	vRes.append(1,res);
	return vRes;
}// end of function

// There are no relations here.
int gainLossAlphabet::relations(const int charInSeq, const int charToCheck) const{
	if (charInSeq == charToCheck)
		return 1;
	return 0;
}

int gainLossAlphabet::fromChar(const string& str, const int pos) const{
	return fromChar(str[pos]);
}




