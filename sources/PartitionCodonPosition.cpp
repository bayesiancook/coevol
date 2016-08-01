

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
	CodonSequenceAlignment* codondata = new CodonSequenceAlignment(data,true,MtMam);

	string base = argv[2];

	for (int pos=0; pos<3; pos++)	{
		ostringstream s;
		s << base << pos+1 << ".ali";
		ofstream os(s.str().c_str());
		codondata->ToStream(os,pos);
	}
}

