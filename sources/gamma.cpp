
#include "Random.h"
#include "IncompleteGamma.h"


int main(int argc, char* argv[])	{

	double alpha = atof(argv[1]);
	int discgam = atoi(argv[2]);

	double* xx = new double[discgam];
	double* yy = new double[discgam];
	double* rr = new double[discgam];

	double lg = Random::logGamma(alpha+1.0);
	for (int i=0; i<discgam; i++)	{
		xx[i] = PointGamma((i+1.0)/discgam,alpha,alpha);
	}
	for (int i=0; i<discgam; i++)	{
		yy[i] = IncompleteGamma(alpha*xx[i],alpha+1,lg);
	}
	yy[discgam-1] = 1.0;
	rr[0] = discgam * yy[0];
	for (int i=1; i<discgam; i++)	{
		rr[i] = discgam * (yy[i] - yy[i-1]);
	}
	for (int i=0; i<discgam; i++)	{
		cout << rr[i] << '\n';
	}
}

