#include "KH_calculation.h"

namespace KH_calculation {

double get_phi (double z)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of z
    int sign = 1;
    if (z < 0)
	{
		sign = -1;
	}
    z = fabs(z)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*z);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-z*z);

    return 0.5*(1.0 + sign*y);
}


double calc_p_value_kh (const Vdouble & LogLikePerPositionA, const Vdouble & LogLikePerPositionB)
{
	//calc esteemated variance of delta of KH (Assessing the Uncertainty in Phylogenetic Inference, Nielsen, pg 484)
	//delta(X) = LL(A) - LL(B)
	//H0: E(delta(X)) <= 0   ---> tree B is either better or equal to tree A
	//H1: E(delta(X)) > 0    ---> tree A is better than tree B
	int num_pos = LogLikePerPositionA.size();
	double varDeltaX = 0;
	double sum_diffs = 0;
	double avg_diff = 0;
	for (int i=0; i < num_pos; ++i)
	{
		sum_diffs += (LogLikePerPositionA[i] - LogLikePerPositionB[i]);
	}
	avg_diff = sum_diffs / num_pos;

	double sum_squares = 0;
	double sqr_diff = 0;
	for (int i=0; i < num_pos; ++i)
	{
		sqr_diff = pow (LogLikePerPositionA[i] - LogLikePerPositionB[i] - avg_diff, 2);
		sum_squares += sqr_diff;
	}
	varDeltaX = (num_pos / (num_pos - 1)) * sum_squares;
	//end calc esteemated variance of delta of KH (Assessing the Uncertainty in Phylogenetic Inference, Nielsen, pg 484)

	//obtain the standard test statistic, z:
	double stdDeltaX = sqrt (varDeltaX);
	double z = sum_diffs / stdDeltaX; //let's hope stdDeltaX is not a zero
	
	double phi_of_z = get_phi (z);
	double p_value = 1 - phi_of_z; //one-sided test to see if A is better than B

	return p_value;
}

};