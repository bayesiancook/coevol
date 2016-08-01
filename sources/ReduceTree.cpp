
#include "Tree.h"

#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char* argv[])	{

	string infile = argv[1];
	ofstream os(argv[2]);
	Tree* tree = new Tree(infile);

	tree->Reduce();
	tree->PrintReduced(os);

}

