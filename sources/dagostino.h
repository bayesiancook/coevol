
using namespace std;

#include <vector>
#include <cmath>
#include <iostream>

double DAgostinosKandZ(vector<double>& ic)	{

	double n = ic.size();

	double m1 = 0;
	for (int i=0; i<n; i++)	{
		m1 += ic[i];
	}
	m1 /= n;

	double m2 = 0;
	double m3 = 0;
	double m4 = 0;

	for (int i=0; i<n; i++)	{
		double tmp = ic[i] - m1;
		m2 += tmp * tmp;
		m3 += tmp * tmp * tmp;
		m4 += tmp * tmp * tmp * tmp;
	}
	m2 /= n;
	m3 /= n;
	m4 /= n;

	// cerr << "mean : " << m1 << '\n';
	// cerr << "var  : " << m2 << '\n';

	double g1 = m3 / sqrt(m2 * m2 * m2);
	double g2 = m4 / m2 / m2 - 3;

	double mu2 = 6.0 * (n-2) / (n+1) / (n+3);
	double gamma2 = 36.0 * (n-7) * (n*n + 2*n -5) / (n-2) / (n+5) / (n+7) / (n+9);
	double W2 = sqrt(2*gamma2 +4)- 1;
	double W = sqrt(W2);
	double delta = 1.0 / sqrt(log(W));
	double alpha = sqrt(2 / (W2 - 1));
	double t = g1 / alpha / sqrt(mu2);
	double Z1 = delta * log(t + sqrt(t*t+1));

	double mmu1 = -6.0 / (n+1);
	double mmu2 = 24.0 * n * (n-2) * (n-3) / (n+1) / (n+1) / (n+3) / (n+5);
	double ggamma1 = 6.0 * (n*n -5*n + 2) / (n+7) / (n+9) * sqrt(6.0 * (n+3) * (n+5) / n / (n-2) / (n-3));
	double A = 6.0 + 8.0 / ggamma1* (2.0 / ggamma1 + sqrt(1 + 4.0 / ggamma1 / ggamma1));
	double Z2 = sqrt(9.0 * A / 2) * (1 - 2.0 / 9 / A - exp( 1.0 / 3 * log( (1 - 2.0 / A) / (1.0 + (g2 - mmu1) / sqrt(mmu2) * sqrt(2.0 / (A - 4))))));

	double ZZ1 = Z1 * Z1;
	double ZZ2 = Z2 * Z2;
	double KK = ZZ1 + ZZ2;
	// cout << KK << '\t' << ZZ1 << '\t' << ZZ2 << '\n';
	return KK;

}
