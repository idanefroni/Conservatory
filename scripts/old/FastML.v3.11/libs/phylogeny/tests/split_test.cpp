using namespace std;
#include "split.h"

int main()
{
  cout << " testing the \"split\" class" <<endl;
    
  cout << "make set - max=5" <<endl;
  split s1(5);
  s1.print();

  cout << "toggle(4)" <<endl;
   {
    s1.reverseMembership(4);
  }
  s1.print();

  cout << "toggle(4)" <<endl;
   {
    s1.reverseMembership(4);
  }
  s1.print();

  cout << "toggle(4)" <<endl;
   {
    s1.reverseMembership(4);
  }
  s1.print();

  cout << "toggle(3)" <<endl;
   {
    s1.reverseMembership(3);
  }
  s1.print();

  cout << "toggle(3)" <<endl;
   {
    s1.reverseMembership(3);
  }
  s1.print();

  cout << "toggle(0)" <<endl;
   {
    s1.reverseMembership(0);
  }
  s1.print();

  cout << "toggle(1);" <<endl;
   {
    s1.reverseMembership(1);
  }
  s1.print();

  cout << "toggle(1);" <<endl;
   {
    s1.reverseMembership(1);
  }
  s1.print();

  cout << "toggle(1);" <<endl;
   {
    s1.reverseMembership(1);
  }
  s1.print();

  cout << "toggle(0)" <<endl;
 {
    s1.reverseMembership(0);
  }
  s1.print();

  cout << "toggle(0)" <<endl;
 {
    s1.reverseMembership(0);
  }
  s1.print();



  // part II - from iterator

  cout <<endl << "test split constractor from iterator"<<endl;
  vector<int> v(3,0);
  v[0]=2;  v[1]=0;   v[2]=4;
  vector<int>::const_iterator vbeg = v.begin();
  vector<int>::const_iterator vend = v.end();
  split s2(vbeg,vend,5);
  s2.print();

  v[0]=2;  v[1]=3;   v[2]=4;
  vbeg = v.begin();
  vend = v.end();
  split s3(vbeg,vend,5);

  
  cout << s3 <<endl;

  cout <<endl<<"Testing competability"<<endl;

  cout << s1<<" and "<<s2<<"\t:"<<s1.compatible(s2)<<endl;
  cout << s1<<" and "<<s3<<"\t:"<<s1.compatible(s3)<<endl;
  
  return (0);
}
