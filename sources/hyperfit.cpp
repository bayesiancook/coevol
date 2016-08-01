
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

double* v;
double* t;
int N;
double v0;
double v1;

double score()	{

	double total = 0;
	for (int i=0; i<N; i++)	{
		double tmp = v[i] - v0 - v1 / t[i];
		// double tmp = log(v[i]) - log(v0 + v1 / t[i]);
		total += t[i] * tmp * tmp;
	}
	return total;

}

double d0()	{

	double total = 0;
	for (int i=0; i<N; i++)	{
		// double tmp = log(v[i]) - log(v0 + v1 / t[i]);
		// total -= tmp / (v0 + v1/t[i]);
		double tmp = v[i] - v0 - v1 / t[i];
		total -= t[i] * tmp;
	}
	return total;

}

double d1()	{

	double total = 0;
	for (int i=0; i<N; i++)	{
		// double tmp = log(v[i]) - log(v0 + v1 / t[i]);
		// total -= tmp / (v0 + v1/t[i]) / t[i];
		double tmp = v[i] - v0 - v1 / t[i];
		total -= t[i] * tmp / t[i];
	}
	return total;
}

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	double delta = atof(argv[2]);
	double epsilon = atof(argv[3]);
	is >> N;
	v = new double[N];
	t = new double[N];
	for (int i=0; i<N; i++)	{
		is >> t[i] >> v[i];
	}

	v0 = 1;
	v1 = 1;
	double s = score();
	int count = 0;
	while ((count < 10000) && (s > epsilon))	{
		count++;
		v0 -= delta * d0();
		s = score();
		cerr << "0 : " << s << '\t' << v0 << '\t' << v1 << '\n';
		v1 -= delta * d1();
		s = score();
		cerr << "1 : " << s << '\t' << v0 << '\t' << v1 << '\n';
	}
}


