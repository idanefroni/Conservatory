// we want to make sure that we are using the doubleRep class with DOUBLEREP enabled.
#define DOUBLEREP t
#include "../doubleRep.cpp"


int main()
{
  double   d=5.352e-30;
  doubleRep k;
  k= d;
  k.output(cout);cout<<endl;
  cout << k.mantissa() <<" "<< k.expon() <<endl<<endl;

  cout <<"as double    "<<d<<endl;
  cout <<"as doubleRep "<<convert(k)<<endl;
  return(0);
}
