
#include "Random.h"

int main(int argc, char* argv[])	{

	double s1 = atof(argv[1]);
	double s2 = atof(argv[2]);
	double alpha = atof(argv[3]);
	int N = atoi(argv[4]);
	double T = atof(argv[5]);

	double x = 0;
	double y = 0;
	for (int i=0; i<N; i++)	{
		
		x += s2 * Random::sNormal();
		if (Random::Uniform() < alpha)	{
			y += s1 * Random::sNormal();
		}
		/*
		double dx = s1 * Random::sNormal();
		x += dx;
		y += alpha * dx + s2 * Random::sNormal();
		// y = x + alpha * (y-x) + (1 - alpha) * s2 * Random::sNormal();
		cout << x << '\t' << y << '\n';
		*/
		cout << T * i / N << '\t' << x << '\t' << y << '\t' << x + y << '\n';
	}
}

