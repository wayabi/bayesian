#ifndef __U_GAUSSIAN__
#define __U_GAUSSIAN__

#include <vector>
#include <ostream>
#include <random>
#include <memory>

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
	std::shared_ptr<matrix> inv() const;
	vec head_mul(const vec& a);
	vec tail_mul(const vec& a);
	matrix cholesky() const;

	friend std::ostream& operator<<(std::ostream& o, const matrix& m);
	friend matrix operator+(const matrix& a, const matrix& b);
	friend matrix operator-(const matrix& a, const matrix& b);

	static vec subtract(const vec& a, const vec& b);
	static double dot(const vec& a, const vec& b);
};

struct gaussian {
public:
	gaussian();
	gaussian(const vec& mu, const matrix& sigma);
	double get_y(const vec& p);
	vec get_random_p();
	double calc_percentile(const vec& p) const;
	double calc_percentile_x(const vec& p) const;

	friend std::ostream& operator<<(std::ostream& o, const gaussian& m);
	static gaussian calc_variance_matrix(int dimension, const vec* data, int num_data);

public:
	vec mu_;
	matrix sigma_;

private:
	static std::vector<double> _2PiRoot;
	static double get_2pi_root(int dimension);
	static double _Percentile[];
	static std::random_device _RandomDevice;
	static std::mt19937 _RandomGen;
	static std::normal_distribution<> _NormalDistribution;
};

#endif
