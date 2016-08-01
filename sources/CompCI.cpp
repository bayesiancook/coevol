
#include <iostream>
#include <cstdlib>
#include <fstream>

using namespace std;


int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int N = atoi(argv[2]);
	int P = atoi(argv[3]);

	double min[N][P];
	double max[N][P];

	for (int i=0; i<N; i++)	{
		for (int j=0; j<P; j++)	{
			is >> min[i][j] >> max[i][j];
		}
	}

	for (int k=0; k<P; k++)	{
		for (int l=k+1; l<P; l++)	{
			double mean = 0;
			double maxx = 0;
			double max1 = 0;
			double max2 = 0;
			for (int i=0; i<N; i++)	{
				double mininf = (min[i][k] > min[i][l]) ? min[i][l] : min[i][k];
				double minsup = (min[i][k] > min[i][l]) ? min[i][k] : min[i][l];
				double maxinf = (max[i][k] > max[i][l]) ? max[i][l] : max[i][k];
				double maxsup = (max[i][k] > max[i][l]) ? max[i][k] : max[i][l];
				double tmp = 2 * (maxinf - minsup) / (max[i][k] - min[i][k] + max[i][l] - min[i][l]);
				mean += tmp;
				if ((!i) || (maxx > tmp))	{
					maxx = tmp;
					max1 = min[i][k];
					max2 = max[i][k];
				}
			}
			mean /= N;
			cout << k << '\t' << l << '\t' << mean << '\t' << maxx << '\t' << max1 << '\t' << max2 << '\n';
		}
	}

}

