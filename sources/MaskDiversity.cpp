#include "SequenceAlignment.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* data = new FileSequenceAlignment(argv[1]);
	ifstream is(argv[2]);
	int n;
	is >> n;
	string* mask = new string[n];
	for (int j=0; j<n; j++)	{
		is >> mask[j];
	}
	data->MaskDiversity(mask,n);

	ofstream os(argv[3]);
	data->ToStream(os);
}

