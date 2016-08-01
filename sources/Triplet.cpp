

#include "SequenceAlignment.h"
#include "ContinuousData.h"
#include "Tree.h"
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;


int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);

	// string base = argv[2];

	cout << data->GetNsite() << '\t' << data->GetNonMissingTriplet() << '\n';

	/*
	ofstream os(argv[2]);
	data->ToStreamTriplet(os);
	*/

}

