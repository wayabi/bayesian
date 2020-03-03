#include <stdio.h>

#include <sstream>
#include <iostream>
#include <random>
#include <vector>

#include <opencv2/core.hpp>
#include <opencv2/aruco.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "gaussian.h"
#include "bayesian.h"
#include "Util.h"
#include "thread_pool.h"

using namespace std;
using namespace cv;
using namespace boost;

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
	g = g.calc_variance_matrix(3, &dd[0], 10);
	cout << "g0:" << g << endl;
	for(int i=0;i<100;++i){
		g = bayes.posterior_outlier1(g, &dd[i*10], 10, 0.45);
		cout << "g_next:" << g << endl;
	}
}

void t1(const string& path)
{
	Mat m = imread(path.c_str());
	resize(m, m, Size(),0.5, 0.5);
	//cout << "channel:" << m.channels() << endl;
	int sub = 4;
	int w = m.cols;
	int h = m.rows;
	int ww = w/sub;
	int hh = h/sub;
	vector<vec> dd;
	for(int yy=0;yy<1;++yy){
		for(int xx=0;xx<1;++xx){
			for(int y=0;y<hh;++y){
				for(int x=0;x<ww;++x){
					int idx = (y*sub+yy) * w * 3 + (x*sub+xx) * 3;
					double r = *(m.data+idx+0);
					double g = *(m.data+idx+1);
					double b = *(m.data+idx+2);
					r += rand()%3-1;
					g += rand()%3-1;
					b += rand()%3-1;
					dd.push_back(vec{r, g, b});
				}
			}
		}
	}

	cout << "num data:" << dd.size() << endl;
	int num_step = 50;
	int num_unit = dd.size() / num_step;

	bayesian bayes(100000);
	gaussian g(vec{0, 0, 0}, matrix(vec{1, 0, 0, 0, 1, 0, 0, 0, 1}, 3));
	g = g.calc_variance_matrix(3, &dd[0], num_unit*10);
	//cout << "g0:" << g << endl;
	for(int i=10;i<num_step;++i){
		g = bayes.posterior_outlier1(g, &dd[i*num_unit], num_unit, 5);
		//cout << "g_next:" << g << endl;
	}
	cout << "g_infelenced:" << g << endl;
	if(g.sigma_.inv() == nullptr) cout << "inv() null" << endl;
	int divide_lut = 1;
	int scale_lut = 256 / divide_lut;
	int size_lut = scale_lut * scale_lut * scale_lut;
	vector<int> lut(size_lut, -1);
	lut.resize(size_lut);
	for(int i=0;i<size_lut;++i){
		lut[i] = -1;
	}
	for(int y=0;y<h;++y){
		for(int x=0;x<w;++x){
			int idx = y*w*3+x*3;
			double rr = *(m.data+idx+0);
			double gg = *(m.data+idx+1);
			double bb = *(m.data+idx+2);
			int idx_lut = scale_lut * scale_lut * ((int)rr / divide_lut) +
				scale_lut * ((int)gg / divide_lut) +
				(int)bb / divide_lut;
			if(lut[idx_lut] < 0){
				double len = g.calc_percentile_x(vec{rr, gg, bb});
				if(len < 3){
					lut[idx_lut] = 1;
				}else{
					lut[idx_lut] = 0;
				}
			}
			if(lut[idx_lut] == 1){
				*(m.data+idx+0) = 0;
				*(m.data+idx+1) = 255;
				*(m.data+idx+2) = 0;
			}
		}
	}

	auto pp = Util::split(path, '/');
	string path_out = "./out/" + pp[pp.size()-1];
	imwrite(path_out.c_str(), m);
	cout << path_out.c_str() << endl;
}

void t1_thread(thread_pool& pool)
{
	auto ff = Util::getListFile("./data");
	for(auto ite = ff.begin();ite != ff.end();++ite){
		pool.post(std::bind(t1, string("./data/") + *ite));
	}
}

int main(int argc, char** argv){
	random_device rnd;
	mt19937 mt(rnd());
	//uniform_real_distribution<> r(-1, 1);
	normal_distribution<> r(0.0, 1.0);

	asio::io_service io_service_;
	thread_pool pool(io_service_, 4);

	t1_thread(pool);
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
