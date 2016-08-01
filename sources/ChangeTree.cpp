
#include <fstream>
using namespace std;
#include "Tree.h"
#include "SequenceAlignment.h"
#include "ContinuousData.h"

int main(int argc, char* argv[])	{

	/*
	Tree tree(argv[1]);
	ofstream os(argv[2]);
	tree.ToStream(os);
	*/

	/*
	FileSequenceAlignment data(argv[1]);
	ofstream os(argv[2]);
	data.ToStream(os);
	*/

	FileContinuousData data(argv[1]);
	ofstream os(argv[2]);
	data.ToStream(os);

}
