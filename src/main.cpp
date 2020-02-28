#include "gaussian.h"
#include <sstream>
#include <iostream>

using namespace std;
int main(int argc, char** argv){
	vec mu{1, 2, 3};
	vec sig{
		1, 0.8, 0.4,
		0.8, 0.1, -0.5,
		0.4, -0.5, 1.8};

	matrix sigma(sig, 3);
	gaussian g(mu, sigma);

	vec p{1, 2, 3};
	double y = g.get_y(p);
	cout << "y:" << y << endl;

	return 0;
}
