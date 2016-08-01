#include "CodonSequenceAlignment.h"
#include "ProteinSequenceAlignment.h"

int main(int argc, char* argv[])	{

	cerr << "opening file : " << argv[1] << '\n';
	FileSequenceAlignment* nuc = new FileSequenceAlignment(argv[1]);
	cerr << "using universal code\n";
	cerr << "making codon aligment\n";
	CodonSequenceAlignment* cod = new CodonSequenceAlignment(nuc,true,Universal);
	cerr << "protein alignment\n";
	ProteinSequenceAlignment* prot = new ProteinSequenceAlignment(cod);
	cerr << "output\n";
	ofstream os(argv[2]);
	prot->ToStream(os);
	cerr << "ok\n";
}

