

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
	ifstream is(argv[2]);
	int N;
	is >> N;
	int begin[N];
	int end[N];
	int length[N];
	for (int i=0; i<N; i++)	{
		is >> begin[i] >> end[i];
		begin[i] --;
;		length[i] = end[i] - begin[i];
	}
	// int N = atoi(argv[2]);
	string base = argv[3];
	int Nsite = data->GetNsite();
	// int nsite = Nsite / N;

	int tot = 0;

	for (int i=0; i<N; i++)	{
		ostringstream s;
		s << base << i  << ".ali";
		ofstream os(s.str().c_str());
		// cerr << s.str() << '\t' << tot << '\t' <<  nsite << '\n';
		// SequenceAlignment* part = new SequenceAlignment(data,tot,nsite);
		SequenceAlignment* part = new SequenceAlignment(data,begin[i],length[i]);
		cerr << "ok\n";
		part->ToStream(os);

		/*
		CodonSequenceAlignment* cod = new CodonSequenceAlignment(part);
		ProteinSequenceAlignment* prot = new ProteinSequenceAlignment(cod);
		ofstream os2((s.str() + ".prot").c_str());
		prot->ToStream(os2);
		*/

		// tot += nsite;
		delete part;
	}
}

