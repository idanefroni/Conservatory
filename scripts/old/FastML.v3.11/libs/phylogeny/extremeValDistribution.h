#ifndef _EXTREME_VAL_DISTRIBUTION
#define _EXTREME_VAL_DISTRIBUTION
#include "definitions.h"
#define EULER_CONSTANT 0.5772156649015328606065120900824024310421593359399235
/*
The extreme value distribution is used to model the largest 
value from a large collection of random observations from the 
same distribution.
1. The distribution has two parameters: 
a location parameter alpha and a scale parameter beta. 
2. The cumulative distribution function (CDF) is:
exp(-exp(-(x-alpha)/beta))
3. Mean: E(x) = alpha + beta*EULER_CONSTANT
   STD(x) = beta * pi / sqrt(6)  
*/
class extremeValDistribution
{
public:
	extremeValDistribution();
	extremeValDistribution(const extremeValDistribution& other);
	extremeValDistribution& operator=(const extremeValDistribution& other);
	virtual ~extremeValDistribution();
	//fit the alpha and beta parameters from a population mean and std
	void fitParametersFromMoments(MDOUBLE exp, MDOUBLE std);

	MDOUBLE getCDF(MDOUBLE score) const;
	MDOUBLE getInverseCDF(MDOUBLE pVal) const;
private:
	MDOUBLE _alpha;
	MDOUBLE _beta;

};
#endif //_EXTREME_VAL_DISTRIBUTION
