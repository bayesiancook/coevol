

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

	if (argc == 1)	{
		cerr << '\n';
		cerr << "comp3 <datafile> <outfile> [mtmam]\n";
		cerr << "computes empirical composition of dataset\n";
		cerr << "if mtmam specified, uses the vertebrate mitochondrial code, otherwise, uses universal code\n";
		cerr << '\n';
		exit(1);
	}
	string datafile = argv[1];

	SequenceAlignment* data = new FileSequenceAlignment(datafile);

	string base = argv[2];
	GeneticCodeType type = Universal;
	if (argc == 4)	{
		string gtype = argv[3];
		if (gtype == "mtmam")	{
			type = MtMam;
			cerr << "mito\n";
		}
	}

	CodonSequenceAlignment* codondata = new CodonSequenceAlignment(data,true,type);
	ofstream thirdcomp4os((base + "_pos3.compACGT").c_str());
	ofstream thirdcomp3os((base + "_pos3.compGCskews").c_str());
	codondata->NucleotideCompositionalHeterogeneity(&thirdcomp4os,2,0,&thirdcomp3os);

	ofstream comp4os((base + ".compACGT").c_str());
	ofstream comp3os((base + ".compGCskews").c_str());
	codondata->NucleotideCompositionalHeterogeneity(&comp4os,-1,0,&comp3os);

	cerr << '\n';
	cerr << "in " << base << ".compACGT\n";
	cerr << "entries are : pi_A pi_C pi_G pi_T\n";
	cerr << '\n';
	cerr << "in " << base << ".compGCskews\n";
	cerr << "entries are : GC/(1-GC) exp(GCskew) exp(ATskew)\n";
	cerr << "      GC      = (pi_G + pi_C)\n";
	cerr << "      GC skew = (pi_G - pi_C) / (pi_G + pi_C)\n";
	cerr << "      AT skew = (pi_A - pi_T) / (pi_A + pi_T)\n";

	cerr << "\n";
	cerr << "in " << base << "_pos3 files: compositions at third coding positions\n";
	cerr << '\n';
}

