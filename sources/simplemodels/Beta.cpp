
#include "BetaChain.h"

int main(int argc, char* argv[])	{


	cerr << "seed was : " << Random::GetSeed() << '\n';
	if (argc == 2)	{
		string filename = argv[1];
		BetaChain chain(filename);
	}
	else if (argc == 5)	{
		string datafile = argv[1];
		string filename = argv[2];
		int burnin = atoi(argv[3]);
		int until = atoi(argv[4]);
		BetaChain chain(datafile,filename,burnin,until,1);
	}
	else	{
		cerr << '\n';
		cerr << "error in command:\n";
		cerr << "to run a new chain : mult <datafile>  <chainname>\n";
		cerr << "to run an already existing chain: mult <chainname>\n";
		cerr << '\n';
	}
}

