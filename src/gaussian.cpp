#include "gaussian.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <random>

using namespace std;

std::vector<double> gaussian::_2PiRoot;
std::random_device gaussian::_RandomDevice;
std::mt19937 gaussian::_RandomGen{_RandomDevice()};
std::normal_distribution<> gaussian::_NormalDistribution{0, 1};
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

double matrix::determinant() const
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

matrix matrix::submatrix(int x, int y) const
{
	vec v;
	for(int yy = 0;yy<dimension_;++yy){
		if(yy == y) continue;
		for(int xx = 0;xx<dimension_;++xx){
			if(xx == x) continue;
			v.push_back(m_[yy*dimension_+xx]);
		}
	}
	return matrix(v, dimension_-1);
}

matrix matrix::inv() const
{
	vec v;
	double det = determinant();
	assert(det > 0);
	for(int y=0;y<dimension_;++y){
		for(int x=0;x<dimension_;++x){
			v.push_back(pow(-1, x+y+2)*submatrix(x, y).determinant() / det);
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

vec matrix::tail_mul(const vec& a)
{
	assert(a.size() == dimension_);
	vec ret;
	for(int y=0;y<dimension_;++y){
		double sum = 0;
		for(int x=0;x<dimension_;++x){
			sum += a[x] * m_[y*dimension_+x];
		}
		ret.push_back(sum);
	}
	return ret;
}

matrix matrix::cholesky() const
{
	double det = determinant();
	vec ret;
	for(int y=0;y<dimension_;++y){
		for(int x=0;x<dimension_;++x){
			ret.push_back(0);
		}
	}

	assert(det > 0);
	for(int y=0;y<dimension_;++y){
		for(int x=0;x<=y;++x){
			if(x==y){
				double sum = 0;
				for(int k=0;k<x;++k){
					sum = ret[y*dimension_+k] * ret[y*dimension_+k];
				}
				ret[y*dimension_+x] = sqrt(m_[y*dimension_+x] - sum);
			}else{
				double sum = 0;
				for(int k=0;k<x;++k){
					sum += ret[y*dimension_+k] * ret[x*dimension_+k];
				}
				ret[y*dimension_+x] = (m_[y*dimension_+x] - sum) / ret[x*dimension_+x];
			}
		}
	}
	return matrix(ret, dimension_);
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

matrix operator+(const matrix& a, const matrix& b)
{
	assert(a.dimension_ == b.dimension_);
	matrix m;
	m.dimension_ = a.dimension_;
	for(int y=0;y<m.dimension_;++y){
		for(int x=0;x<m.dimension_;++x){
			m.m_.push_back(a.m_[y*m.dimension_+x] + b.m_[y*m.dimension_+x]);
		}
	}
	return m;
}
matrix operator-(const matrix& a, const matrix& b)
{
	assert(a.dimension_ == b.dimension_);
	matrix m;
	m.dimension_ = a.dimension_;
	for(int y=0;y<m.dimension_;++y){
		for(int x=0;x<m.dimension_;++x){
			m.m_.push_back(a.m_[y*m.dimension_+x] - b.m_[y*m.dimension_+x]);
		}
	}
	return m;
}

gaussian::gaussian()
{
}

gaussian::gaussian(const vec& mu, const matrix& sigma)
{
	assert(mu.size() == sigma.dimension_);
	mu_ = mu;
	sigma_ = sigma;
}

double gaussian::get_y(const vec& p)
{
	assert(p.size() == mu_.size());
	double det = sigma_.determinant();
	assert(det > 0);
	double a = 1.0 / (get_2pi_root(p.size()) * sqrt(det));
	//cout << "a:" << a << endl;
	//cout << "2pi_root:" << get_2pi_root(x.size()) << endl;
	//cout << "sigma_determinant:" << sigma_.determinant() << endl;
	//cout << sigma_inv << endl;
	vec sub = matrix::subtract(p, mu_);
	double y = a * exp(-0.5 * matrix::dot(sigma_.inv().head_mul(sub), sub));
	return y;
}

vec gaussian::get_random_p()
{
	matrix chol = sigma_.cholesky();
	int num = mu_.size();
	vec p;
	for(int i=0;i<num;++i){
		double d = _NormalDistribution(_RandomGen);
		p.push_back(d);
	}
	p = chol.tail_mul(p);
	for(int i=0;i<num;++i){
		p[i] = p[i] + mu_[i];
	}
	return p;
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

gaussian gaussian::calc_variance_matrix(int dimension, const vec* data, int num_data)
{
	assert(dimension > 0 && num_data > 0);
	vec mu;
	for(int i=0;i<dimension;++i){
		mu.push_back(0);
	}
	for(int i=0;i<num_data;++i){
		for(int j=0;j<dimension;++j){
			mu[j] += (*(data+i))[j];
		}
	}
	for(int i=0;i<dimension;++i){
		mu[i] /= num_data;
	}

	vector<vec> center;
	for(int i=0;i<num_data;++i){
		vec c = *(data+i);
		for(int j=0;j<dimension;++j){
			c[j] -= mu[j];
		}
		center.push_back(c);
	}

	vec m;
	for(int y=0;y<dimension;++y){
		for(int x=0;x<dimension;++x){
			m.push_back(0);
		}
	}
	for(int i=0;i<num_data;++i){
		for(int y=0;y<dimension;++y){
			for(int x=0;x<dimension;++x){
				m[y*dimension+x] += center[i][x] * center[i][y];
			}
		}
	}
	for(int y=0;y<dimension;++y){
		for(int x=0;x<dimension;++x){
			m[y*dimension+x] /= num_data;
		}
	}
	return gaussian(mu, matrix(m, dimension));
}

std::ostream& operator<<(std::ostream& o, const gaussian& m)
{
	o << "mu:[";
	for(auto ite = m.mu_.begin();ite != m.mu_.end();++ite){
		o << *ite << ", ";
	}
	o << "]" << endl;
	o << m.sigma_ << endl;
	return o;
}

double gaussian::calc_percentile(const vec& p) const
{
	vec sub = matrix::subtract(p, mu_);
	double x = -0.5 * matrix::dot(sigma_.inv().head_mul(sub), sub);
	cout << "x:" << x << endl;
	int idx = (int)(-100 * x);
	double per = 0;
	if(idx >= 500){
		per = 0.5;
	}else{
		per = _Percentile[idx];
	}
	return per;
}

double gaussian::calc_percentile_x(const vec& p) const
{
	vec sub = matrix::subtract(p, mu_);
	double x = 0.5 * matrix::dot(sigma_.inv().head_mul(sub), sub);
	return x;
}
double gaussian::_Percentile[]{.0000,.0040,.0080,.0120,.0160,.0199,.0239,.0279,.0319,.0359,.0398,.0438,.0478,.0517,.0557,.0596,.0636,.0675,.0714,.0753,.0793,.0832,.0871,.0910,.0948,.0987,.1026,.1064,.1103,.1141,.1179,.1217,.1255,.1293,.1331,.1368,.1406,.1443,.1480,.1517,.1554,.1591,.1628,.1664,.1700,.1736,.1772,.1808,.1844,.1879,.1915,.1950,.1985,.2019,.2054,.2088,.2123,.2157,.2190,.2224,.2257,.2291,.2324,.2357,.2389,.2422,.2454,.2486,.2517,.2549,.2580,.2611,.2642,.2673,.2704,.2734,.2764,.2794,.2823,.2852,.2881,.2910,.2939,.2967,.2995,.3023,.3051,.3078,.3106,.3133,.3159,.3186,.3212,.3238,.3264,.3289,.3315,.3340,.3365,.3389,.3413,.3438,.3461,.3485,.3508,.3531,.3554,.3577,.3599,.3621,.3643,.3665,.3686,.3708,.3729,.3749,.3770,.3790,.3810,.3830,.3849,.3869,.3888,.3907,.3925,.3944,.3962,.3980,.3997,.4015,.4032,.4049,.4066,.4082,.4099,.4115,.4131,.4147,.4162,.4177,.4192,.4207,.4222,.4236,.4251,.4265,.4279,.4292,.4306,.4319,.4332,.4345,.4357,.4370,.4382,.4394,.4406,.4418,.4429,.4441,.4452,.4463,.4474,.4484,.4495,.4505,.4515,.4525,.4535,.4545,.4554,.4564,.4573,.4582,.4591,.4599,.4608,.4616,.4625,.4633,.4641,.4649,.4656,.4664,.4671,.4678,.4686,.4693,.4699,.4706,.4713,.4719,.4726,.4732,.4738,.4744,.4750,.4756,.4761,.4767,.4772,.4778,.4783,.4788,.4793,.4798,.4803,.4808,.4812,.4817,.4821,.4826,.4830,.4834,.4838,.4842,.4846,.4850,.4854,.4857,.4861,.4864,.4868,.4871,.4875,.4878,.4881,.4884,.4887,.4890,.4893,.4896,.4898,.4901,.4904,.4906,.4909,.4911,.4913,.4916,.4918,.4920,.4922,.4925,.4927,.4929,.4931,.4932,.4934,.4936,.4938,.4940,.4941,.4943,.4945,.4946,.4948,.4949,.4951,.4952,.4953,.4955,.4956,.4957,.4959,.4960,.4961,.4962,.4963,.4964,.4965,.4966,.4967,.4968,.4969,.4970,.4971,.4972,.4973,.4974,.4974,.4975,.4976,.4977,.4977,.4978,.4979,.4979,.4980,.4981,.4981,.4982,.4982,.4983,.4984,.4984,.4985,.4985,.4986,.4986,.4987,.4987,.4987,.4988,.4988,.4989,.4989,.4989,.4990,.4990,.4990,.4991,.4991,.4991,.4992,.4992,.4992,.4992,.4993,.4993,.4993,.4993,.4994,.4994,.4994,.4994,.4994,.4995,.4995,.4995,.4995,.4995,.4995,.4996,.4996,.4996,.4996,.4996,.4996,.4997,.4997,.4997,.4997,.4997,.4997,.4997,.4997,.4997,.4997,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4998,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.4999,.49991,.49992,.49992,.49992,.49992,.49992,.49993,.49993,.49993,.49994,.49994,.49994,.49994,.49995,.49995,.49995,.49995,.49995,.49996,.49996,.49996,.49996,.49996,.49996,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49997,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49998,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.49999,.499995,.499995,.499995,.499995,.499995,.499995,.499995,.499995,.499995,.499995};
