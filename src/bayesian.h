#ifndef __U_BAYESIAN__
#define __U_BAYESIAN__

#include "gaussian.h"

class bayesian {
public:
	static gaussian posterior_simple0(const gaussian& g, vec x, vec speed_elastic, vec speed_transition);
	static gaussian posterior_outlier0(const gaussian& g, vec x, vec speed_elastic, vec speed_transition, double outlier_percentage);

};

#endif
