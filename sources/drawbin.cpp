
#include "Random.h"

int main(int argc, char* argv[])	{

	int N = atoi(argv[1]);
	int n = atoi(argv[2]);
	double p = atof(argv[3]);
	double lim = atof(argv[4]);
	double f = 0;
	for (int i=0; i<N; i++)	{
		double x = 0;
		for (int j=0; j<n; j++)	{
			if (Random::Uniform() < p)	{
				x++;
			}
		}
		x /= n;
		x *= 100;
		cout << x << '\n';
		if (x <= lim)	{
			f++;
		}
	}
	cerr << '\n';
	cerr << f/N << '\n';
}
