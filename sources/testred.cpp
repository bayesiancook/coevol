
#include "Tree.h"
#include "ReduceTree.h"

int main(int argc, char* argv[])	{

	string filename = argv[1];
	Tree tree(filename);

	string* names = new string[2];
	names[0] = "B";
	names[1] = "C";
	TaxonSet* taxset = new TaxonSet(names,2);

	ReduceTree red(&tree,taxset);

	red.ToStream(cerr);

}

