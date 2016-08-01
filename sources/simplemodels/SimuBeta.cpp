
#include "BetaChain.h"

int main(int argc, char* argv[])	{


	cerr << "seed was : " << Random::GetSeed() << '\n';
	if (argc == 5)	{
		int N = atoi(argv[1]);
		int n = atoi(argv[2]);
		double hyper = atof(argv[3]);
		string filename = argv[4];
		BetaModel mod(N,n,hyper);
		ofstream os((filename+".simu").c_str());
		mod.PrintSimu(os);
		ofstream dos((filename+".data").c_str());
		mod.PrintData(dos);
		os << '\n';
		os.close();
	}
	else	{
		cerr << '\n';
		cerr << "error in command:\n";
		cerr << '\n';
	}
}

