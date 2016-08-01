#include "SequenceAlignment.h"
#include "Random.h"

int main(int argc, char* argv[])	{

	FileSequenceAlignment* data = new FileSequenceAlignment(argv[1]);
	ifstream is(argv[2]);
	int n1,n2;
	is >> n1 >> n2;
	string* check = new string[n1];
	string* mask = new string[n2];
	for (int j=0; j<n1; j++)	{
		is >> check[j];
	}
	for (int j=0; j<n2; j++)	{
		is >> mask[j];
	}
	data->MaskOutgroup(check,mask,n1,n2);

	ofstream os(argv[3]);
	data->ToStream(os);
}

