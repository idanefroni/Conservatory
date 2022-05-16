using namespace std;
//#include "bootstrap.h"
#include "splitMap.h"

int main()
{
    
  // create a split one way
  split s1(5);
  s1.reverseMembership(0);
  s1.reverseMembership(1);
  s1.reverseMembership(4);


  // and an other split the other way

  vector<int> v(3,0);
  v[0]=2;  v[1]=0;   v[2]=3;
  vector<int>::const_iterator vbeg = v.begin();
  vector<int>::const_iterator vend = v.end();
  split s2(vbeg,vend,5);
  
  cout << endl << "Test the splitMap" << endl;

  splitMap sm1;

  cout <<"s1: ";
  s1.print();
  cout <<"s2: ";
  s2.print();
  cout << endl;

  cout <<"add s1"<<endl;
  sm1.add(s1);
  sm1.print();

  cout <<"add s2"<<endl;
  sm1.add(s2);
  sm1.print();

  cout <<"add s1"<<endl;
  sm1.add(s1);
  sm1.print();

  cout <<"add s1"<<endl;
  sm1.add(s1);
  sm1.print();

  cout <<"add s1"<<endl;
  sm1.add(s1);
  sm1.print();
  cout << endl;

  // print test
  cout << "print test"<<endl;
  cout << sm1;
  cout << endl;
  
  // reverse

  cout << "reverse the map"<<endl;

  vector<pair<split,int> >  rmap =  sm1.sortSplits();
  for (vector<pair<split,int> >::const_iterator i=rmap.begin();i!=rmap.end();++i)
    cout <<i->second<<" "<<i->first<<endl;


  return (0);
}
