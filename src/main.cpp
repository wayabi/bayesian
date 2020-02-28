#include "gaussian.h"
#include <sstream>
#include <iostream>
#include <random>

using namespace std;

int main(int argc, char** argv){
	random_device rnd;
	mt19937 mt(rnd());
	uniform_real_distribution<> r(-1, 1);

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
	return 0;
}
