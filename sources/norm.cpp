
#include "Random.h"
#include "dagostino.h"

int main(int argc, char* argv[])	{

	int N = atoi(argv[1]);
	int P = atoi(argv[2]);
	double var = atof(argv[3]);
	double obs = atof(argv[4]);
	double p = 0;
	for (int j=0; j<P; j++)	{
		vector<double> ic;
		for (int i=0; i<N; i++)	{
			ic.push_back(sqrt(var) * Random::sNormal());
		}
		if (DAgostinosKandZ(ic) > obs) 	{ // 6.34)	{
			p++;
		}
	}
	cout << "number of rejections : " << 100 * p / P  << '\n';


	/*
	for (int i=0; i<N; i++)	{
		for (int j=0; j<P; j++)	{
			cout << sqrt(var) * Random::sNormal() << '\t';
		}
		cout << '\n';
	}
	*/
}

