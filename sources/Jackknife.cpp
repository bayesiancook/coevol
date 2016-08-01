#include "CodonSequenceAlignment.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* nuc = new FileSequenceAlignment(argv[1]);
	CodonSequenceAlignment* cod = new CodonSequenceAlignment(nuc,true);
	ofstream os(argv[2]);
	double p = atof(argv[3]);
	cod->ToStreamRandomJackknife(os,p);
}

