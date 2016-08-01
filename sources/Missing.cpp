#include "SequenceAlignment.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* data = new FileSequenceAlignment(argv[1]);
	ifstream is(argv[2]);
	int N;
	int ngroup;
	is >> N >> ngroup;
	if (N != data->GetNtaxa())	{
		cerr << "error when reading groups\n";
		exit(1);
	}
	map<string,int> group;
	for (int i=0; i<N; i++)	{
		string temp;
		int tmp;
		is >> temp >> tmp;
		group[temp] = tmp;
		cerr << temp << '\t' << tmp << '\n';
	}

	SequenceAlignment* data2 = new SequenceAlignment(data,group,ngroup);
	cerr << "number of sites remaining: " << data2->GetNsite() << '\n';
	ofstream os(argv[3]);
	data2->ToStream(os);
}

