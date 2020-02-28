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
	double determinant() const;
	matrix submatrix(int x, int y) const;
	matrix inv() const;
	vec head_mul(const vec& a);
	friend std::ostream& operator<<(std::ostream& o, const matrix& m);

	static vec subtract(const vec& a, const vec& b);
	static double dot(const vec& a, const vec& b);
};

struct gaussian {
public:
	gaussian();
	gaussian(const vec& mu, const matrix& sigma);
	double get_y(vec p);
	vec get_random_x();
	void calc_variance_matrix(int dimension, std::vector<vec> data);
	double calc_percentile(const vec& p) const;

	friend std::ostream& operator<<(std::ostream& o, const gaussian& m);

public:
	vec mu_;
	matrix sigma_;

private:
	static std::vector<double> _2PiRoot;
	static double get_2pi_root(int dimension);
	static double _Percentile[];
};

#endif
