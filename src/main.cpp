#include <stdio.h>

#include <sstream>
#include <iostream>
#include <random>
#include <vector>

#include "gaussian.h"
#include "bayesian.h"
#include "Util.h"

using namespace std;

void t0(){
	auto ss = Util::load_csv("./out/a");
	vector<vec> dd;
	int num = ss.size();
	for(int i=0;i<num;++i){
		vec d;
		d.push_back(atof(ss[i][0].c_str()));
		d.push_back(atof(ss[i][1].c_str()));
		d.push_back(atof(ss[i][2].c_str()));
		dd.push_back(d);
	}

	cout << "num data:" << num << endl;
	bayesian bayes(1000);
	gaussian g(vec{0, 0, 0}, matrix(vec{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3));
	g = g.calc_variance_matrix(3, &dd[0], 30);
	cout << "g0:" << g << endl;
	for(int i=0;i<100;++i){
		g = bayes.posterior_outlier1(g, &dd[i*10], 10, 0.45);
		cout << "g_next:" << g << endl;
	}
}

int main(int argc, char** argv){
	random_device rnd;
	mt19937 mt(rnd());
	//uniform_real_distribution<> r(-1, 1);
	normal_distribution<> r(0.0, 1.0);

	t0();
	return 0;
	
	vec mu{1, 2, 3};
	vec sig{
		1, 0.8, 0.4,
		0.8, 0.8, 0.5,
		0.4, 0.5, 1.8};

	matrix sigma(sig, 3);
	gaussian g(mu, sigma);

	vec p{1, 2, 3};
	double y = g.get_y(p);
	//cout << "y:" << y << endl;

	for(int i=0;i<100;++i){
		vec pp{1+r(mt), 2+r(mt), 3+r(mt)};
		double yy = g.get_y(pp);
		//cout << "yy:" << yy << endl;
	}
	//////
	/*
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
	*/

	//cout << "cholesky:" << g.sigma_.cholesky() << endl;

	for(int i=0;i<1000;++i){
		vec p = g.get_random_p();
		double per = g.calc_percentile(p);
		cout << "per:" << per << endl;
		//cout << p[0] << "," << p[1] << "," << p[2] << endl;
	}




	return 0;
}
