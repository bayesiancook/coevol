
# include <cstdlib>
# include <cstdio>
# include <cmath>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

double f(double alpha, double gamma, double bgc, double x)	{

	return (gamma * x + bgc) * exp(-gamma * x) / (exp(bgc) - exp(-gamma * x));
	// return (gamma * x + bgc) * exp(-gamma * x) / (exp(bgc) - exp(-gamma * x)) * exp((alpha - 1) * log(x));

}

int main(int argc, char* argv[])	{

	string quadfile = argv[1];
	int maxorder = atoi(argv[2]);
	double max = atoi(argv[3]);
	double alpha = atof(argv[4]);
	double gamma = atof(argv[5]);
	double bgc = atof(argv[6]);

	double k1 = 0;

	for (int order = maxorder/10; order<= maxorder; order+=maxorder / 10)	{

		double* x = new double[order];
		
		double total = 0;
		for (int i=1; i<order; i++)	{
			double x = ((double) max * i) / order;
			total += f(alpha,gamma,bgc,x) * exp(-x);
		}
		total *= max / order;
		cout.precision(20);
		cout << order << '\t' << total  << '\n';
		k1 = total;

		delete[] x;
	}

	int quadorder;
	double k2 = 0;
	ifstream is(quadfile.c_str());
	is >> quadorder;
	double* x = new double[quadorder];
	double* w = new double[quadorder];
	for (int i=0; i<quadorder; i++)	{
		is >> w[i] >> x[i];
	}
	
	double total = 0;
	for (int i=0; i<quadorder; i++)	{
		total += w[i] * f(alpha,gamma,bgc,x[i]);
	}
	cout.precision(20);
	cout << quadorder << '\t' << total << '\n';
	k2 = total;

	delete[] x;
	delete[] w;
	cout << '\n';
	cout << fabs(k2-k1)/ k2<< '\n';
	cout << '\n';
	return 0;

}

