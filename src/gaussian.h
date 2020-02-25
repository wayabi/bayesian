#ifndef __U_GAUSSIAN__
#define __U_GAUSSIAN__

#include <vector>

typedef std::vector<double> vec;
#define SQRT2 = 2.50662827463;

struct gaussian {
public:
	double get_y(vec x);
	vec get_random_x();

public:
	vec mu;
	vec sigma;
};

#endif
