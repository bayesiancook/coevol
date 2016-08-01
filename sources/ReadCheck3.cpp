
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int n = atoi(argv[2]);
	int k1 = atoi(argv[3]);
	int k2 = atoi(argv[4]);
	int k = atoi(argv[5]);

	double q = atof(argv[6]);

	double tmp;
	int tot1 = 0;
	int tot2 = 0;
	for (int i=0; i<n; i++)	{
		for (int j=0; j<k1; j++)	{
			is >> tmp;
		}
		for (int j=k1; j<k2; j++)	{
			is >> tmp;
			if (tmp < q/2)	{
				tot1++;
			}
			if (tmp > (1 - q/2))	{
				tot2++;
			}
		}
		for (int j=k2; j<k; j++)	{
			is >> tmp;
		}
	}
	cout << "obs : " << tot1 + tot2 << '\n';
	cout << "inf : " << tot1 << '\n';
	cout << "sup : " << tot2 << '\n';
	cout << "exp : " << n * (k2 - k1) * q << '\n';
}




