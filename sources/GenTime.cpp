
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[])	{

	int N;
	int K;
	ifstream is(argv[1]);
	ofstream os(argv[2]);

	is >> N >> K;
	os << N << '\t' << 5 << '\n';
	if (K != 4)	{
		cerr << "error \n";
		exit(1);
	}
	for (int i=0; i<N; i++)	{
		string name;
		double mat, mass, longe, met;
		is >> name >> mat >> mass >> longe >> met;
		double gen;
		if ((mat != -1)	&& (longe != -1))	{
			gen = 0.5 * (mat + 365 * longe);
		}
		else	{
			gen = -1;
		}
		os << name << '\t' << mat << '\t' << longe << '\t' << gen << '\t' << mass << '\t' << met << '\n';
	}
	/*
	for (int i=0; i<N; i++)	{
		string name;
		double mat, mass, longe, met;
		is >> name >> mat >> mass >> longe >> met;
		if (met != -1)	{
			if (mass == -1)	{
				met = -1;
			}
			else	{
				met /= mass;
			}
		}
		os << name << '\t' << mat << '\t' << mass << '\t' << longe << '\t' << met << '\n';
	}
	*/
}

