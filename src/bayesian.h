#ifndef __U_BAYESIAN__
#define __U_BAYESIAN__

#include <deque>
#include <vector>
#include <memory>
#include "gaussian.h"

class bayesian {
public:
	static gaussian posterior_simple0(const gaussian& g, const vec& p, double effectiveness);
	static gaussian posterior_outlier0(const gaussian& g, const vec& p, double effectiveness, double outlier_sigma);

public:
	bayesian(int size_data_queue);
	void clear_data();
	gaussian posterior_outlier1(const gaussian& g, const vec* data, int num_data, double outlier_percentile);

private:
	int size_data_queue_;
	std::deque<vec> data_;
};

#endif
