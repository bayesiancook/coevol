
#include "Tree.h"
#include <fstream>
#include <cstdlib>
using namespace std;


int main(int argc, char* argv[])	{

	string treefile = argv[1];
	Tree* tree = new Tree(treefile);
		cerr << "before : " << tree->GetFullSize(tree->GetRoot()) << '\n';
	tree->Subdivide(tree->GetRoot(),2);
		cerr << "after  : " << tree->GetFullSize(tree->GetRoot()) << '\n';
	tree->ToStream(cerr);


}


