
#include "Tree.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	string treefile = argv[1];
	ofstream os(argv[2]);

	Tree* tree = new Tree(treefile);
	
	double min = tree->GetMinHeight(tree->GetRoot());
	double max = tree->GetMaxHeight(tree->GetRoot());

	cerr << min << '\t' << max << '\n';

	tree->ToStreamRenorm(tree->GetRoot(),os,1.0 / max);

}
