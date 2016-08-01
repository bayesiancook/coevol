
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

using namespace std;


int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);

	double mean1[N];
	double stdev1[N];
	double mean2[N];
	double stdev2[N];
	double min1[N];
	double max1[N];
	double min2[N];
	double max2[N];
	double* zs = new double[N];
	string left1[N];
	string right1[N];
	string left2[N];
	string right2[N];

	for (int i=0; i<N; i++)	{
		is >> left1[i] >> right1[i] >> mean1[i] >> stdev1[i] >> min1[i] >> max1[i] >> left2[i] >> right2[i] >> mean2[i] >> stdev2[i] >> min2[i] >> max2[i];
		if (left1[i] != left2[i])	{
			cerr << "error : non matching names: " << left1[i] << '\t' << left2[i] << '\n';
		}
		if (right1[i] != right2[i])	{
			cerr << "error : non matching names: " << right1[i] << '\t' << right2[i] << '\n';
		}
		if (stdev1[i] + stdev2[i] > 1e-6)	{
			zs[i] = 2 * (mean1[i] - mean2[i]) / (stdev1[i] + stdev2[i]);
		}
		else	{
			zs[i] = 0;
		}
	}

	for (int i=0; i<N; i++)	{
		double max = 0;
		int imax = 0;
		for (int j=0; j<N; j++)	{
			if (max<fabs(zs[j]))	{
				max = fabs(zs[j]);
				imax = j;
			}
		}
		if (max > 1e-6)	{
			cout << zs[imax] << '\t' << left1[imax] << '\t' << right1[imax] << '\t' << mean1[imax] << '\t' << stdev1[imax] << '\t' << exp(min1[imax]) << '\t' << exp(max1[imax]) << '\t' << mean2[imax] << '\t' << stdev2[imax] << '\t' << exp(min2[imax]) << '\t' << exp(max2[imax]) << '\n';
		}
		zs[imax] = 0;
	}

}

