#include <sstream>
#include <iostream>
#include <random>

#include "gaussian.h"
#include "bayesian.h"

using namespace std;

int main(int argc, char** argv){
	random_device rnd;
	mt19937 mt(rnd());
	//uniform_real_distribution<> r(-1, 1);
	normal_distribution<> r(0.0, 1.0);

	vec mu{1, 2, 3};
	vec sig{
		1, 0.8, 0.4,
		0.8, 0.8, 0.5,
		0.4, 0.5, 1.8};

	matrix sigma(sig, 3);
	gaussian g(mu, sigma);

	vec p{1, 2, 3};
	double y = g.get_y(p);
	cout << "y:" << y << endl;

	for(int i=0;i<100;++i){
		vec pp{1+r(mt), 2+r(mt), 3+r(mt)};
		double yy = g.get_y(pp);
		cout << "yy:" << yy << endl;
	}
	//////
	vector<vec> data;
	int dimension = 3;
	for(int i=0;i<100;++i){
		vec pp{1+r(mt), 2+r(mt)*3, 3+r(mt)};
		data.push_back(pp);
	}
	g.calc_variance_matrix(3, data);
	cout <<  "g:" << g << endl;
	for(int i=0;i<100;++i){
		vec pp{1+r(mt), 2+r(mt)*3, 3+r(mt)};
		//cout << "percentile:" << g.calc_percentile(pp) << " (" << pp[0] << ", " << pp[1] << ", " << pp[2] << ")" << endl;
	}

	cout << "cholesky:" << g.sigma_.cholesky() << endl;

	for(int i=0;i<100;++i){
		vec p = g.get_random_p();
		cout << p[0] << "," << p[1] << "," << p[2] << endl;
	}
	return 0;
}
