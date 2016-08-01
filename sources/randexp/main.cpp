
#include "Random.h"

int main()	{

	cerr << "this is a uniform random number\n";
	double u = Random::Uniform();
	cerr << u << '\n';
	cerr << '\n';

	cerr << "this is a standard (unit mean) exponential random number\n";
	double e = Random::sExpo();
	cerr << e << '\n';
	cerr << '\n';

}


