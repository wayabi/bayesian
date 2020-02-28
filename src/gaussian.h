#ifndef __U_GAUSSIAN__
#define __U_GAUSSIAN__

#include <vector>
#include <ostream>

typedef std::vector<double> vec;

class matrix {
public:
	matrix();
	matrix(const std::vector<double>& m, int dimension);
	std::vector<double> m_;
	int dimension_;
public:
	double determinant();
	matrix submatrix(int x, int y);
	matrix inv();
	vec head_mul(const vec& a);
	friend std::ostream& operator<<(std::ostream& o, const matrix& m);

	static vec subtract(const vec& a, const vec& b);
	static double dot(const vec& a, const vec& b);
};

struct gaussian {
public:
	gaussian(const vec& mu, const matrix& sigma);
	double get_y(vec x);
	vec get_random_x();

public:
	vec mu_;
	matrix sigma_;

private:
	static std::vector<double> _2PiRoot;
	static double get_2pi_root(int dimension);
};

#endif
