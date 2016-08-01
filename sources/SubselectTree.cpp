
#include "Tree.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	string treefile = argv[1];
	ifstream is(argv[2]);
	ofstream os(argv[3]);

	int Ntaxa;
	is >> Ntaxa;
	string* speciesnames = new string[N];
	for (int i=0; i<N; i++)	{
		is >> speciesnames[i];
	}

	Tree* tree = new Tree(treefile);
	
	tree->SelectTaxa(N,speciesnames);

}
