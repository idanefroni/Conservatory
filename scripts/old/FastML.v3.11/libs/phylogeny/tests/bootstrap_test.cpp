#include "bootstrap.h"
#include "treeUtil.h"
#include "someUtil.h"
using namespace std;

int main()
{
  cout << "creating a bootstrap object from a file"<<endl;
  
  string filename("bootstrap_test.txt");

  vector<tree> tv(getStartingTreeVecFromFile(filename));

  // first constractor
  cout << " first constractor"<<endl;
  bootstrap b1(filename);
  b1.print();
  cout <<endl;

  // secound constractor
  cout << " secound constractor" <<endl;
  bootstrap b2(tv);
  b2.print();
  cout <<endl;
  
  cout << "getting weights from a tree" << endl;
  map<int,MDOUBLE> v1(b1.getWeightsForTree(tv[0])) ;
  for (map<int,MDOUBLE>::iterator i = v1.begin();i!=v1.end();++i)
    cout << " "<<i->second;
  cout << endl;

  cout << "print the support of a tree" <<endl;
  b1.printTreeWithBPvalues(cout, tv[0], v1);
  cout <<endl <<endl;
  

  cout<< "remove the first tree from the list, and use is as bases for additional computation"<<endl;
  tree t(tv[0]);
  cout<< "use the secound tree twice"<<endl;
  tv[0]=tv[1];
  
  // secound constractor
  bootstrap b3(tv);
  b3.print_names(cout);
  b3.print();
  map<int,MDOUBLE> v3(b3.getWeightsForTree(t)) ;
  //  for (map<int,MDOUBLE>::iterator i = v3.begin();i!=v3.end();++i)
  //  cout << " "<<i->second;
  //cout << endl;
  
  cout << "print the support of the removed tree"<<endl;
  b3.printTreeWithBPvalues(cout, t, v3);
  cout <<endl;

  cout <<endl<<endl<<endl<<"compatability"<<endl;
  tree t2(b3.consensusTree());
  //  cout << t2<<endl;
  map<int,MDOUBLE> support(b3.getWeightsForTree(t2));


  b3.printTreeWithBPvalues(cout, t2, support);
  cout <<endl;

  //  for (map<int,MDOUBLE>::const_iterator ii= support.begin(); ii != support.end();++ii)
  // cout << ii->second <<" ";
  // cout << endl;

  cout <<"compatability 0.0"<<endl;
  t=b3.consensusTree(0.);
  support=b3.getWeightsForTree(t);
  b3.printTreeWithBPvalues(cout, t, support);
  cout <<endl;

//   for (map<int,MDOUBLE>::iterator i=support.begin();i!=support.end();++i)
//     {
//       cout << "<"<<i->first<<","<<i->second <<">:"<<support[i->first]<<endl;
//     }

  double c=0.8;
  cout <<"compatability "<<c<<endl;
  t=b3.consensusTree(c);
  support=b3.getWeightsForTree(t);
  b3.printTreeWithBPvalues(cout, t, support);
  cout <<endl;
//   for (map<int,MDOUBLE>::iterator i=support.begin();i!=support.end();++i)
//     {
//       cout << "<"<<i->first<<","<<i->second <<">:"<<support[i->first]<<endl;
//     }

//   for (map<int,MDOUBLE>::const_iterator i=support.begin();i!=support.end();++i)
//     {
//       cout << "<"<<i->first<<","<<">:"<<support[i->first]<<endl;
//     }

  
  return (0);
}
