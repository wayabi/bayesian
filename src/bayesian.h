#ifndef __U_BAYESIAN__
#define __U_BAYESIAN__

#include "gaussian.h"

class bayesian {
public:
	static gaussian posterior_simple0(const gaussian& g, const vec& p, double effectiveness);
	static gaussian posterior_outlier0(const gaussian& g, const vec& p, double effectiveness, double outlier_percentile);

};

#endif
