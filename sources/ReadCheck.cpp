
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

void Cred95(int b, double p, int& infb, int& supb)	{

	double threshold = 0.05;

	double logmult[b+1];
	double logp = log(p);
	double logq = log(1 - p);

	double tot = 0;
	double max = 0;
	for (int i=0; i<=b; i++)	{
		logmult[i] = i * logp + (b-i) * logq;
		for (int j=0; j<i; j++)	{
			logmult[i] += log(b-j) - log(j+1);
		}

		if ((i == 0) || (max < logmult[i]))	{
			max = logmult[i];
		}
	}
	for (int i=0; i<b; i++)	{
		tot += exp(logmult[i] - max);
	}
	double mult[b+1];
	for (int i=0; i<=b; i++)	{
		mult[i] = exp(logmult[i] - max) / tot;
	}

	double total = 0;
	infb = 0;
	while (total < threshold)	{
		total += mult[infb];
		infb++;
	}

	total = 0;
	supb = b;
	while (total < threshold)	{
		total += mult[supb];
		supb--;
	}
}

int main(int argc, char* argv[])	{

	ifstream is(argv[1]);
	int n = atoi(argv[2]);
	int k = atoi(argv[3]);
	int ncat = atoi(argv[4]);

	double p[k][ncat];
	for (int cat=0; cat<ncat; cat++)	{
		for (int j=0; j<k; j++)	{
			p[j][cat] = 0;
		}
	}

	for (int i=0; i<n; i++)	{
		for (int j=0; j<k; j++)	{
			double tmp;
			is >> tmp;
			for (int cat = 0; cat<ncat; cat++)	{
				// is it in the x%-CI where x = 10*cat + 5 ?
				// is the value > x/2 and < 1-x/2 ?
				double cutoff = ((double) cat + 0.5) / ncat / 2;
				if ((tmp >= cutoff) && (tmp <= 1-cutoff))	{
					p[j][cat] ++;
				}
			}
		}
	}

	cout << '\n';
	for (int cat=0; cat<ncat; cat++)	{
		double x = (1 - ((double) cat + 0.5) / ncat);
		cout << 100*x << '\t';
		for (int j=0; j<k; j++)	{

			// expected
			int n1,n2;
		 	Cred95(n,x,n1,n2);

			// observed
			if ((p[j][cat] < n1) || (p[j][cat]>n2))	{
				cout << "*";
			}
			else	{
				cout << " ";
			}
			cout << p[j][cat] << '\t';

			cout << "[" << n1 << ":" << n2 << "]";
			cout << '\t';

		}
		cout << '\n';
	}

	cout << '\n' << '\n';

	for (int cat=0; cat<ncat; cat++)	{
		double x = (1 - ((double) cat + 0.5) / ncat);
		cout << 100*x << '\t';
		for (int j=0; j<k; j++)	{
			// expected
			int n1,n2;
		 	Cred95(n,x,n1,n2);

			// observed
			if ((p[j][cat] < n1) || (p[j][cat]>n2))	{
				cout << "*";
			}
			else	{
				cout << " ";
			}
			cout << ((int) (100 * p[j][cat] / n)) << '\t';
			cout << '\t';
		}
		cout << '\n';
	}
}




