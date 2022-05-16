// $Id: readTreeWithComments.cpp 753 2006-06-29 14:41:56Z ninio $


#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
using namespace std;

#include "logFile.h"
#include "tree.h"
//#include "readTree.h"

int main(int argc,char*argv[]) {
  if(argc<2) exit(1);
  if(argc>2) myLog::setLog("-",atoi(argv[2]));
  string treeName(argv[1]);
  tree t(treeName);
  t.output(cout);

  vector<tree::nodeP> nv;
  t.getAllNodes(nv, t.getRoot());
  cout <<"got "<<nv.size()<<" noded"<<endl;
  for (vector<tree::nodeP>::iterator i=nv.begin();i!=nv.end();++i)
	cout << (*i)->getComment()<<endl;
  exit(0);
}
