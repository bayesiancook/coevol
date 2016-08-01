
#include "Random.h"
#include "stabledistribution.h"

double logNormalDensity(double scale, double x)	{
	return -0.5 * log(2 * Pi) - log(scale) - 0.5 * x * x / scale / scale;
}

double logCauchyDensity(double scale, double x)	{
	return -log (Pi * scale * (1 + x*x/scale/scale));
}

int main(int argc, char* argv[])	{

	double T = atof(argv[1]);
	double dt = atof(argv[2]);
	double alpha = atof(argv[3]);

	StableDistribution as(alpha,1.0);

	int NT = (int) (T/dt);
	for (int i=0; i<=NT; i++)	{
		double t = i*dt;
		double tmp1 = as.logPDF(t);
		double tmp2 = 0;
		if (alpha == 1.0)	{
			tmp2 = logCauchyDensity(1.0,t);
		}
		else if (alpha == 2.0)	{
			tmp2 = logNormalDensity(sqrt(2.0),t);
		}
		double error = fabs(tmp1 - tmp2);
		cerr << t << '\t' << error << '\t' << tmp1 << '\t' << tmp2 << '\n';
	}
	
}

