#include "bayesian.h"
#include <iostream>
#include <assert.h>

using namespace std;

gaussian bayesian::posterior_simple0(const gaussian& g, const vec& p, double effectiveness)
{
	double per = g.calc_percentile(p);

	matrix sig0 = g.sigma_;
	vec c = matrix::subtract(p, g.mu_);
	int d = g.sigma_.dimension_;
	for(int y=0;y<d;++y){
		for(int x=0;x<d;++x){
			sig0.m_[y*d+x] += effectiveness * c[x]*c[y];
		}
	}
	vec mu0 = g.mu_;
	for(int i=0;i<d;++i){
		mu0[i] = (1-effectiveness)*mu0[i] + effectiveness*p[i];
	}
	return gaussian(mu0, sig0);
}

gaussian bayesian::posterior_outlier0(const gaussian& g, const vec& p, double effectiveness, double outlier_percentile)
{
	double per = g.calc_percentile(p);
	if(per > outlier_percentile) return g;

	matrix sig0 = g.sigma_;
	int d = g.sigma_.dimension_;

	vec c = matrix::subtract(p, g.mu_);
	for(int y=0;y<d;++y){
		for(int x=0;x<d;++x){
			//sig0.m_[y*d+x] = (1-effectiveness) * sig0.m_[y*d+x] + effectiveness * c[x]*c[y];
			sig0.m_[y*d+x] = ((99.0) * sig0.m_[y*d+x] + c[x]*c[y])/100.0;
			//sig0.m_[y*d+x] += effectiveness * c[x]*c[y];
		}
	}
	vec mu0 = g.mu_;
	for(int i=0;i<d;++i){
		mu0[i] = (1-effectiveness)*mu0[i] + effectiveness*p[i];
	}
	return gaussian(mu0, sig0);
}

bayesian::bayesian(int size_data_queue)
{
	size_data_queue_ = size_data_queue;
}

void bayesian::clear_data()
{
	data_.clear();
}

gaussian bayesian::posterior_outlier1(const gaussian& g, const vec* data, int num_data, double outlier_percentile)
{
	assert(num_data > 0);
	for(int i=0;i<num_data;++i){
		vec v = *(data+i);
		double perx = g.calc_percentile_x(v);
		if(perx < 4){
			//cout << "per:" << per << endl;
			data_.push_back(v);
			if(data_.size() > size_data_queue_){
				data_.pop_front();
			}
		}
	}
	vector<vec> dd;
	for(auto ite=data_.begin();ite!=data_.end();++ite){
		dd.push_back(*ite);
	}
	cout << "size dd:" << dd.size() << endl;
	if(dd.size() < 3) return g;
	return gaussian::calc_variance_matrix((data+0)->size(), &dd[0], dd.size());
	
}
