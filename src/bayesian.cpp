#include "bayesian.h"
#include <iostream>

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
	if(per < 0.45) return g;

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

