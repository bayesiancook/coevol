
#include "BipoissonChain.h"

int main(int argc, char* argv[])	{


	cerr << "seed was : " << Random::GetSeed() << '\n';
	if (argc == 2)	{
		string filename = argv[1];
		BipoissonChain chain(filename);
	}
	else if (argc == 6)	{
		string datafile = argv[1];
		string filename = argv[2];
		int every = atoi(argv[3]);
		int until = atoi(argv[4]);
		int conjugate = atoi(argv[5]);
		BipoissonChain chain(datafile,filename,until,every,conjugate,1);
	}
	else if (argc == 4)	{
		int N = atoi(argv[1]);
		int P = atoi(argv[2]);
		string name = argv[3];
		BipoissonModel model(N,P);
		model.Simulation(name);
	}
	else	{
		cerr << '\n';
		cerr << "error in command:\n";
		cerr << "to run a new chain : mult <datafile>  <chainname>\n";
		cerr << "to run an already existing chain: mult <chainname>\n";
		cerr << '\n';
	}
}

