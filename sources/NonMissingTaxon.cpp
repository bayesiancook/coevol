#include "SequenceAlignment.h"
#include "CodonSequenceAlignment.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* data = new FileSequenceAlignment(argv[1]);
    string taxon = argv[2];

    GeneticCodeType type = Universal;
    CodonSequenceAlignment* codondata = new CodonSequenceAlignment(data,true,type);
    SequenceAlignment* codondata2 = new SequenceAlignment(codondata,taxon);

	cerr << "number of sites remaining: " << codondata2->GetNsite() << '\n';
	ofstream os(argv[3]);
	codondata2->ToStream(os);

}

