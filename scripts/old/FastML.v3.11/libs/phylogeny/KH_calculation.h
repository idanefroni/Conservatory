// Kishino-Hasegawa Test 2013 02 27 Eli Levy Karin

#ifndef ___KH_CALCULATION
#define ___KH_CALCULATION

#include "math.h"
#include <cmath>
#include "definitions.h"

namespace KH_calculation {
	double calc_p_value_kh (const Vdouble & LogLikePerPositionA, const Vdouble & LogLikePerPositionB);
	
	double get_phi (double z);
};

#endif
