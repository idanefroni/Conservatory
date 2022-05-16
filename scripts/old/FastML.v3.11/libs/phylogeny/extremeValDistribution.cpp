#include "extremeValDistribution.h"
#include <cmath>
using namespace std;



extremeValDistribution::extremeValDistribution()
: _alpha(0), _beta(0)
{
}

extremeValDistribution::extremeValDistribution(const extremeValDistribution& other)
{
	_alpha = other._alpha;
	_beta = other._beta;
}

extremeValDistribution& extremeValDistribution::operator=(const extremeValDistribution& other)
{
	_alpha = other._alpha;
	_beta = other._beta;
	return *this;
}

extremeValDistribution::~extremeValDistribution()
{
}

/*fits the _alpha and _beta parameters based on a population mean and std.
Based on the following arguments:
1. If variable Z has a cumulative distribution F(Z) = exp(-exp((-z))
     Then E(Z) = EULER_CONSTANT
	 Var(Z) = pi^2/6
2. We assign Z = (X-_alpha) / _beta --> X = _beta*Z + _alpha
   and we get:
   E(X) = _beta*E(Z)+_alpha = _beta*EULER_CONSTANT+_alpha
   Var(X) = _beta^2*pi^2/6
3. We can now find _alpha and _beta based on the method of moments:
   mean = _beta*EULER_CONSTANT+_alpha
   s = _beta * pi / sqrt(6)
4. And solve:
   _beta = s * qsrt(6) / pi
   _alpha = mean - _beta*EULER_CONSTANT
*/
void extremeValDistribution::fitParametersFromMoments(MDOUBLE mean, MDOUBLE s)
{
	_beta = s * sqrt(6.0) / PI;
	_alpha = mean - (_beta * EULER_CONSTANT);
}

MDOUBLE extremeValDistribution::getCDF(MDOUBLE score) const
{
	MDOUBLE res = exp(-exp(-(score-_alpha) / _beta));
	return res;
}

//get y such that pVal = CDF(y):
//	pVal = exp(-exp(-(y-alpha)/beta))
//	ln(-ln(pVal)) = -(y-alpha)/beta
//	y = alpha - beta*ln(-ln(pVal))
MDOUBLE extremeValDistribution::getInverseCDF(MDOUBLE pVal) const
{
	MDOUBLE res = _alpha - _beta * log(-log(pVal));
	return res;
}

