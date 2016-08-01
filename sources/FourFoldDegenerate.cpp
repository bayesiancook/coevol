

#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"
#include "ContinuousData.h"
#include "Tree.h"
#include <fstream>
#include <cstdlib>
#include <sstream>

using namespace std;


int main(int argc, char* argv[])	{

	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);

	string base = argv[2];
	GeneticCodeType type = Universal;
	if (argc == 4)	{
		string gtype = argv[3];
		if (gtype == "mtmam")	{
			type = MtMam;
		}
	}

	CodonSequenceAlignment* codondata = new CodonSequenceAlignment(data,true,type);
	ostringstream s0;
	s0 << base << 3 << ".ali";
	ofstream os0(s0.str().c_str());
	codondata->ToStream(os0,2);

	ostringstream s;
	s << base << 3 << "_4fold.ali";
	ofstream os(s.str().c_str());
	codondata->ToStreamFourFoldThird(os);

	ostringstream s2;
	s2 << base << 3 << "_4foldwocpg.ali";
	ofstream os2(s2.str().c_str());
	codondata->ToStreamFourFoldThirdwoCpG(os2);

	/*
	ostringstream s3;
	s3 << base << 3 << "_4foldtriplet.ali";
	ofstream os3(s3.str().c_str());
	codondata->ToStreamFourFoldTriplet(os3);
	*/
}

