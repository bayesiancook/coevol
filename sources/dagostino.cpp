#include "dagostino.h"

#include <iostream>
#include <fstream>

int main(int argc, char* argv[])	{


	ifstream is(argv[1]);
	int N;
	is >> N;
	vector<double> ic;
	for (int i=0; i<N; i++)	{
		double tmp;
		is >> tmp;
		ic.push_back(tmp);
	}
	cout << "vector size: " << ic.size() << '\n';
	cout << "dago\n\n";
	cout << DAgostinosKandZ(ic) << '\n';
}

