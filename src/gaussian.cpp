#include "gaussian.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

std::vector<double> gaussian::_2PiRoot;

matrix::matrix()
{
	dimension_ = 0;
}

matrix::matrix(const std::vector<double>& m, int dimension)
{
	assert(m.size() == dimension * dimension);
	m_ = m;
	dimension_ = dimension;
}

double matrix::determinant()
{
	if(dimension_ <= 0) return 0;
	if(dimension_ == 1) return m_[0];
	if(dimension_ == 2){
		return m_[0]*m_[3]-m_[1]*m_[2];
	}
	if(dimension_ == 3){
		return
			m_[0]*m_[4]*m_[8] +
			m_[1]*m_[5]*m_[6] +
			m_[2]*m_[3]*m_[7] -

			m_[2]*m_[4]*m_[6] -
			m_[1]*m_[3]*m_[8] -
			m_[0]*m_[5]*m_[7];
	}
	double sum = 0;
	for(int i=0;i<dimension_;++i){
		sum += pow(-1, i*dimension_+1) * m_[dimension_*i+0] * submatrix(i, 0).determinant();
	}
	return sum;
}

matrix matrix::submatrix(int x, int y)
{
	vec v;
	for(int xx = 0;xx<dimension_;++xx){
		if(xx == x) continue;
		for(int yy = 0;yy<dimension_;++yy){
			if(yy == y) continue;
			v.push_back(m_[xx*dimension_+yy]);
		}
	}
	return matrix(v, dimension_-1);
}

matrix matrix::inv()
{
	vec v;
	double det = determinant();
	for(int x=0;x<dimension_;++x){
		for(int y=0;y<dimension_;++y){
			v.push_back(pow(-1, x+y+2)*submatrix(x, y).determinant());
		}
	}
	return matrix(v, dimension_);
}

vec matrix::subtract(const vec& a, const vec& b)
{
	assert(a.size() == b.size());
	int size = a.size();
	vec ret;
	for(int i=0;i<size;++i){
		ret.push_back(a[i] - b[i]);
	}
	return ret;
}

double matrix::dot(const vec& a, const vec& b)
{
	assert(a.size() == b.size());
	int size = a.size();
	double sum = 0;
	for(int i=0;i<size;++i){
		sum += a[i]*b[i];
	}
	return sum;
}

//transpose vec. then mul
vec matrix::head_mul(const vec& a)
{
	assert(a.size() == dimension_);
	vec ret;
	for(int y=0;y<dimension_;++y){
		double sum = 0;
		for(int x=0;x<dimension_;++x){
			sum += a[x] * m_[x*dimension_+y];
		}
		ret.push_back(sum);
	}
	return ret;
}

std::ostream& operator<< (std::ostream& o, const matrix& m)
{
	o << "[" << endl;
	for(int y=0;y<m.dimension_;++y){
		for(int x=0;x<m.dimension_;++x){
			o << m.m_[y*m.dimension_+x] << ", ";
		}
		o << endl;
	}
	o << "]" << endl;
	return o;
}

gaussian::gaussian(const vec& mu, const matrix& sigma)
{
	assert(mu.size() == sigma.dimension_);
	mu_ = mu;
	sigma_ = sigma;
}

double gaussian::get_y(vec x)
{
	assert(x.size() == mu_.size());

	double a = 1.0 / (get_2pi_root(x.size()) * sqrt(sigma_.determinant()));
	cout << "a:" << a << endl;
	cout << "2pi_root:" << get_2pi_root(x.size()) << endl;
	cout << "sigma_determinant:" << sigma_.determinant() << endl;
	matrix sigma_inv = sigma_.inv();
	cout << sigma_inv << endl;
	vec sub = matrix::subtract(x, mu_);
	double y = a * exp(-0.5 * matrix::dot(sigma_inv.head_mul(sub), sub));
	return y;
}

vec gaussian::get_random_x()
{
	return vec();
}

double gaussian::get_2pi_root(int dimension)
{
	if(_2PiRoot.size() == 0){
		for(int i=0;i<20;++i){
			_2PiRoot.push_back(pow(2*M_PI, i / 2.0));
		}
	}
	assert(dimension >= 0 && dimension < _2PiRoot.size());
	return _2PiRoot[dimension];
}
